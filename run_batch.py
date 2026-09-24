#!/usr/bin/env python3
"""Submit a GLORI batch to Slurm. Run --help for public options."""
import argparse
import hashlib
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import uuid

from batch_common import (ANALYSIS, BLAS_ENV, DEFAULTS, VERSION, BatchLock, load_samples,
                          now, read_json, write_json, limit_math_threads)
from batch_report import collect, read_table, SITE_FIELDS, FDR_FIELDS

ROOT = Path(__file__).resolve().parent


def positive_cpus(value):
    value = int(value)
    if value < 3:
        raise argparse.ArgumentTypeError('至少申请 3 CPU，以容纳排序主线程和额外线程。')
    return value


def positive(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError('必须至少为 1。')
    return value


def scheduler(args):
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True,
                            errors='replace', timeout=60)
    if result.returncode:
        raise RuntimeError(f'{args[0]} 失败: {(result.stderr or result.stdout).strip()}')
    return result.stdout.strip()


def assert_idle(batch):
    """Require all recorded jobs to have left squeue before retries/manual collection."""
    identifiers = [str(item[key]) for item in batch['rounds']
                   for key in ['array_job', 'collect_job'] if item.get(key)]
    if not identifiers:
        return
    # Query user's live jobs without -j: completed IDs can make some squeue versions exit nonzero.
    output = scheduler(['squeue', '--noheader', '--user', str(os.getuid()), '--format=%F'])
    active = set(output.split()) & set(identifiers)
    if active:
        raise ValueError('本批次还有运行、排队或挂起的作业: ' + ', '.join(sorted(active))
                         + '。请等它们结束后再重试或重新汇总。')


def scheduler_states(batch):
    if not shutil.which('sacct'):
        return {}
    mapping = {}
    for item in batch['rounds']:
        if not item.get('array_job'):
            continue
        try:
            output = scheduler(['sacct', '--noheader', '--parsable2', '-j', str(item['array_job']),
                                '--format=JobID%100,State%40,ExitCode'])
        except (RuntimeError, subprocess.TimeoutExpired):
            continue
        for line in output.splitlines():
            fields = [field.strip() for field in line.split('|')]
            if len(fields) < 3:
                continue
            match = re.fullmatch(str(item['array_job']) + r'_(\d+)', fields[0])
            if match and int(match[1]) < len(item['tasks']):
                mapping[item['tasks'][int(match[1])]] = fields[1] + ' / ' + fields[2]
    return mapping


def snapshot(destination):
    destination.mkdir()
    for name in ['run_batch.py', 'batch_common.py', 'batch_worker.py', 'batch_report.py']:
        shutil.copy2(ROOT / name, destination / name)
    shutil.copytree(ROOT / 'core', destination / 'core',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in ['LICENSE', 'THIRD_PARTY_NOTICES.md']:
        shutil.copy2(ROOT / name, destination / name)
    shutil.copytree(ROOT / 'LICENSES', destination / 'LICENSES')
    checksums = {path.relative_to(destination).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in destination.rglob('*') if path.is_file()}
    write_json(destination / 'sha256.json', checksums)


def script_text(python, path, round_path, mode, environment):
    lines = ['#!/bin/bash', 'set -euo pipefail']
    for key, value in environment.items():
        lines.append('export ' + key + '=' + shlex.quote(value))
    lines.append('export SLURM_EXPORT_ENV=ALL')
    command = shlex.join([python, str(path), mode, str(round_path)])
    if mode == '--worker':
        command += ' --index "${SLURM_ARRAY_TASK_ID:?Missing array task index}"'
    lines.append('exec ' + command)
    return '\n'.join(lines) + '\n'


def successful(attempt):
    try:
        if read_json(attempt / 'status.json')['state'] != 'SUCCESS':
            return False
        read_table(attempt / 'all_A_sites.tsv', SITE_FIELDS)
        read_table(attempt / 'significant_sites.tsv', ['Sample'] + FDR_FIELDS)
        return True
    except (OSError, ValueError, KeyError):
        return False


def resource_options(args, previous=None):
    settings = dict(previous or DEFAULTS)
    for key in DEFAULTS:
        value = getattr(args, key, None)
        if value is not None:
            settings[key] = value
    if not re.fullmatch(r'[1-9][0-9]*[MGT]', settings['memory']):
        raise ValueError('--memory 格式例如 16G、32000M。')
    if not re.fullmatch(r'(?:\d+-)?\d+:\d{2}:\d{2}', settings['time']):
        raise ValueError('--time 格式例如 08:00:00 或 1-00:00:00。')
    for key in ['partition', 'account']:
        if settings[key] and not re.fullmatch(r'[A-Za-z0-9_.:,/+-]+', settings[key]):
            raise ValueError(f'--{key} 包含不支持的字符。')
    return settings


def submit(args):
    from batch_worker import dependencies
    samples = load_samples(args.samples)
    if args.validate_inputs:
        print(f'输入检查通过：{len(samples)} 个样品。仅检查路径、参考和 FASTQ 首条记录；未提交作业。')
        return 0
    resources = resource_options(args)
    versions = dependencies(slurm=True)
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if not conda_prefix or Path(sys.prefix).resolve() != Path(conda_prefix).resolve():
        raise ValueError('请先 conda activate glori_amplicon，并使用该环境中的 python。')
    if args.check_only:
        print(f'输入及依赖检查通过：{len(samples)} 个样品；未提交作业。')
        return 0
    if not args.retry_failed and not resources['partition']:
        raise ValueError('Specify --partition with the CPU partition name on your cluster.')
    root = Path(args.out).resolve()
    if '%' in str(root) or '\n' in str(root):
        raise ValueError('输出目录不能含 % 或换行（避免 Slurm 日志路径展开错误）。')
    if args.retry_failed:
        if not (root / 'batch.json').is_file():
            raise ValueError('--retry-failed 需要已存在的批次目录及 batch.json。')
    else:
        if root.exists():
            raise ValueError(f'输出目录已存在，拒绝覆盖: {root}。请换新目录或使用 --retry-failed。')
        root.mkdir(parents=True)
    with BatchLock(root):
        if args.retry_failed:
            batch = read_json(root / 'batch.json')
            if batch['version'] not in ('1.0.0', VERSION) or batch['samples'] != samples:
                raise ValueError('样品表、输入文件或版本与原批次不同。请使用新的 --out 目录。')
            assert_idle(batch)
            if any(item.get('uncertain_submission') for item in batch['rounds']):
                raise ValueError('存在未确认的提交记录。请先按 submit.log 检查 Slurm 作业，避免重复提交。')
            tasks = [sample for sample in samples if not successful(root / batch['latest'][sample['sample']])]
            if not tasks:
                print('全部样品均已成功，无需重试。')
                return 0
            resources = resource_options(args, batch['resources'])
            if not resources['partition']:
                raise ValueError('Specify --partition with the CPU partition name on your cluster.')
            # Stale reports must not masquerade as the current retry's result.
            archive = root / f'report_before_retry_{len(batch["rounds"]):03d}'
            archive.mkdir()
            for filename in ['report.html', 'qc_summary.tsv', 'all_A_sites.tsv', 'significant_sites.tsv']:
                if (root / filename).exists():
                    shutil.copy2(root / filename, archive / filename)
        else:
            tasks = samples
            snapshot(root / 'software')
            batch = dict(version=VERSION, created=now(), samples=samples, latest={}, rounds=[],
                         analysis=ANALYSIS, root=str(root))
        batch['resources'] = resources
        round_number = len(batch['rounds']) + 1
        round_dir = root / f'rounds/{round_number:03d}'
        round_dir.mkdir(parents=True)
        (root / 'logs').mkdir(exist_ok=True)
        task_entries = []
        for sample in tasks:
            entry = dict(sample)
            relative = f'samples/{sample["sample"]}/attempt_{round_number:03d}'
            attempt = root / relative
            attempt.mkdir(parents=True)
            write_json(attempt / 'status.json', dict(state='QUEUED', sample=sample['sample']))
            entry['attempt'] = relative
            batch['latest'][sample['sample']] = relative
            task_entries.append(entry)
        environment = {key: os.environ[key] for key in ['PATH', 'CONDA_PREFIX', 'CONDA_DEFAULT_ENV',
                       'LD_LIBRARY_PATH'] if key in os.environ}
        environment.update(BLAS_ENV)
        environment['PYTHONUNBUFFERED'] = '1'
        round_path = round_dir / 'round.json'
        plan = dict(root=str(root), tasks=task_entries, resources=resources, analysis=ANALYSIS,
                    python=sys.executable, environment=environment, versions=versions, created=now())
        write_json(round_path, plan)
        for mode, name in [('--worker', 'array.sh'), ('--collect', 'collect.sh')]:
            (round_dir / name).write_text(script_text(sys.executable, root / 'software/run_batch.py',
                round_path, mode, environment), encoding='utf-8', newline='\n')
        entry = dict(round=round_number, tasks=[sample['sample'] for sample in tasks],
                     array_job=None, collect_job=None, phase='prepared')
        batch['rounds'].append(entry)
        write_json(root / 'batch.json', batch)
        collect(root, batch, final=False)  # Mark queued jobs honestly, not as failed analysis.
        common = ['--parsable', '--export=ALL', '--no-requeue']
        for key in ['partition', 'account']:
            if resources[key]:
                common.append(f'--{key}={resources[key]}')
        token = uuid.uuid4().hex[:10]
        entry['job_name'] = 'glori_' + token

        def sbatch(options, filename, key):
            command = ['sbatch'] + common + options + [str(round_dir / filename)]
            entry.update(phase=f'submitting_{key}', uncertain_submission=True)
            write_json(root / 'batch.json', batch)
            try:
                output = scheduler(command)
            except RuntimeError:
                # sbatch rejected the submission; a timeout instead stays uncertain.
                entry['uncertain_submission'] = False
                raise
            with (round_dir / 'submit.log').open('a', encoding='utf-8') as log:
                log.write(shlex.join(command) + '\n' + output + '\n')
            match = re.fullmatch(r'(\d+)(?:;[^\s]+)?', output)
            if not match:
                raise ValueError(f'sbatch 返回了无法确认的作业编号: {output!r}。请检查 {entry["job_name"]} 作业。')
            entry[key] = match[1]
            entry['uncertain_submission'] = False
            write_json(root / 'batch.json', batch)

        try:
            sbatch(['--hold', '--job-name=' + entry['job_name'],
                    f'--array=0-{len(tasks)-1}%{resources["max_parallel"]}', '--ntasks=1',
                    f'--cpus-per-task={resources["cpus"]}', f'--mem={resources["memory"]}',
                    f'--time={resources["time"]}', '--chdir=' + str(root),
                    '--output=' + str(root / 'logs/array-%A_%a.log')], 'array.sh', 'array_job')
            sbatch(['--job-name=' + entry['job_name'] + '_report', '--ntasks=1', '--cpus-per-task=1',
                    '--mem=2G', '--time=01:00:00', '--dependency=afterany:' + entry['array_job'],
                    '--chdir=' + str(root), '--output=' + str(root / 'logs/report-%j.log')],
                   'collect.sh', 'collect_job')
            # Release only after both job IDs have been durably saved.
            scheduler(['scontrol', 'release', entry['array_job']])
            entry['phase'] = 'submitted'
            write_json(root / 'batch.json', batch)
        except Exception as error:
            entry.update(phase='submission_failed', error=str(error))
            for key in ['collect_job', 'array_job']:
                if entry.get(key):
                    try:
                        scheduler(['scancel', entry[key]])
                    except Exception as cancel_error:
                        entry['cancel_error'] = str(cancel_error)
            for sample in task_entries:
                write_json(root / sample['attempt'] / 'status.json',
                           dict(state='FAILED', sample=sample['sample'], error='提交失败: ' + str(error)))
            write_json(root / 'batch.json', batch)
            collect(root, batch)
            raise
        print(f'已提交 {len(tasks)} 个样品。样品数组作业: {entry["array_job"]}；汇总作业: {entry["collect_job"]}')
        print(f'分析与汇总分区: {resources["partition"]}；具体节点由 Slurm 分配。')
        print(f'最多并行 {resources["max_parallel"]} 个样品，每样品 {resources["cpus"]} CPU / {resources["memory"]}。')
        print(f'可关闭终端；作业结束后查看 {root / "report.html"}')
    return 0


def main():
    limit_math_threads()
    parser = argparse.ArgumentParser(description='GLORI 多样品自动分析（Slurm）')
    parser.add_argument('--samples', help='三列 Tab 分隔样品表')
    parser.add_argument('--out', help='新建的批次结果目录')
    parser.add_argument('--max-parallel', type=positive, help='最多同时运行的样品数，默认 4')
    parser.add_argument('--cpus', type=positive_cpus, help='每样品 CPU，默认 12')
    parser.add_argument('--memory', help='每样品内存，默认 16G')
    parser.add_argument('--time', help='每样品时限，默认 08:00:00')
    parser.add_argument('--partition', help='Slurm CPU partition for analysis and reports; required for a new batch')
    parser.add_argument('--account', help='Slurm 计费账户；省略时使用默认账户')
    parser.add_argument('--keep-all', action='store_true', default=None, help='保留所有可保留的中间文件')
    parser.add_argument('--retry-failed', action='store_true', help='只重跑失败/中断样品，保留旧尝试')
    parser.add_argument('--check-only', action='store_true', help='只检查输入及依赖，不提交')
    parser.add_argument('--validate-inputs', action='store_true', help='仅检查输入路径/参考/首条 FASTQ；支持 Windows')
    parser.add_argument('--summarize', action='store_true', help='在所有作业结束后重新生成报告')
    parser.add_argument('--worker', help=argparse.SUPPRESS)
    parser.add_argument('--index', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--collect', help=argparse.SUPPRESS)
    parser.add_argument('--version', action='version', version=VERSION)
    args = parser.parse_args()
    if args.worker:
        from batch_worker import work
        if args.index is None or args.index < 0:
            parser.error('worker requires non-negative --index')
        return work(args.worker, args.index)
    if args.collect or args.summarize:
        if args.collect:
            root = Path(read_json(args.collect)['root'])
        elif args.out:
            root = Path(args.out).resolve()
        else:
            parser.error('--summarize 需要 --out。')
        with BatchLock(root):
            batch = read_json(root / 'batch.json')
            if args.summarize:
                assert_idle(batch)
            failures = collect(root, batch, scheduler_states(batch))
        print(f'报告已生成: {root / "report.html"}；未成功样品数: {failures}')
        # Report generation succeeded even if a sample failed; sample failures are in the report.
        return 0
    if not args.samples or not (args.out or args.check_only or args.validate_inputs):
        parser.error('请提供 --samples 和 --out；检查模式可省略 --out。')
    return submit(args)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError, ImportError, subprocess.SubprocessError) as error:
        print('错误: ' + str(error), file=sys.stderr)
        sys.exit(1)
