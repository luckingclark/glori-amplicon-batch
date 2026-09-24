"""One sample per Slurm array task; all paths are passed as argument lists."""
import hashlib
import importlib
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback
import zipfile

from batch_common import (fingerprint, fastq_stats, now, read_json, write_json,
                          limit_math_threads, validate_trim_paths)
from batch_report import export_sites

TOOLS = {'bowtie': ['--version'], 'bowtie-build': ['--version'], 'samtools': ['--version'],
         'trim_galore': ['--version'], 'cutadapt': ['--version'], 'fastqc': ['--version']}


def dependencies(slurm=False):
    # Importing numpy/scipy can initialize OpenBLAS before any analysis is run.
    limit_math_threads()
    if sys.platform != 'linux':
        raise ValueError('实际分析和 Slurm 提交需要 Linux；Windows 可使用 --validate-inputs。')
    versions = dict(python=sys.version, executable=sys.executable)
    for name in ['numpy', 'pandas', 'scipy', 'statsmodels', 'Bio', 'pysam']:
        module = importlib.import_module(name)
        versions[name] = getattr(module, '__version__', 'unknown')
    from statsmodels.stats.multitest import multipletests
    from scipy.stats import binomtest
    for name, arguments in TOOLS.items():
        executable = shutil.which(name)
        if executable is None:
            raise ValueError(f'找不到 {name}，请先激活完整的 Conda 环境。')
        result = subprocess.run([executable] + arguments, capture_output=True, text=True,
                                errors='replace', check=True, timeout=60)
        versions[name] = (result.stdout + result.stderr).strip()
    if slurm:
        for name in ['sbatch', 'squeue', 'scontrol', 'scancel']:
            if not shutil.which(name):
                raise ValueError(f'找不到 Slurm 命令 {name}；请在集群登录节点提交。')
    return versions


def fastqc_summary(directory):
    archives = list(Path(directory).glob('*_fastqc.zip'))
    reports = list(Path(directory).glob('*_fastqc.html'))
    if len(archives) != 1 or len(reports) != 1:
        raise ValueError(f'FastQC 报告未完整生成: {directory}')
    with zipfile.ZipFile(archives[0]) as archive:
        names = [name for name in archive.namelist() if name.endswith('/summary.txt')]
        if len(names) != 1:
            raise ValueError(f'FastQC ZIP 缺少 summary.txt: {archives[0]}')
        rows = [line.split('\t') for line in archive.read(names[0]).decode('utf-8').splitlines()]
    if not rows or any(len(row) != 3 or row[0] not in ['PASS', 'WARN', 'FAIL'] for row in rows):
        raise ValueError('FastQC summary.txt 格式不正确。')
    return '; '.join(f'{row[0]}:{row[1]}' for row in rows if row[0] != 'PASS') or 'PASS'


def cleanup_disposable(attempt, sample):
    """Only remove exact files in this attempt after validated exports. No recursive delete."""
    directory = Path(attempt).resolve() / 'output' / 'tmp'
    names = [f'{sample}_AGchanged_2.fq', f'{sample}_A.bed_sorted', f'{sample}_un_2.fq',
             f'{sample}.pileup', f'{sample}.referbase.mpi']
    # Per-contig MPI files can also be huge; keep formatted counts/CR/candidates.
    names += [p.name for p in directory.glob(f'{sample}.referbase.mpi.*')
              if '.formatted.' not in p.name]
    removed = []
    for name in sorted(set(names)):
        path = directory / name
        if path.is_file():
            if path.is_symlink() or path.resolve().parent != directory.resolve():
                raise ValueError(f'临时文件超出本次工作目录，拒绝清理: {path}')
            path.unlink()
            removed.append(name)
    return removed


def work(round_path, index):
    plan = read_json(round_path)
    sample = plan['tasks'][index]
    root = Path(plan['root'])
    attempt = root / sample['attempt']
    # Array jobs must never run the same attempt twice (e.g. automatic requeue).
    claim = attempt / '.started'
    with claim.open('x') as handle:
        handle.write(now() + '\n')
    limit_math_threads()
    os.environ['PYTHONUNBUFFERED'] = '1'
    started = time.monotonic()
    state = dict(state='RUNNING', sample=sample['sample'], start=now(), qc={'warnings': []},
                 slurm_job=os.environ.get('SLURM_JOB_ID'), stage='初始化')
    commands = []
    logdir = attempt / 'logs'
    logdir.mkdir(exist_ok=True)
    status_path = attempt / 'status.json'

    def save():
        state['seconds'] = round(time.monotonic() - started, 3)
        write_json(status_path, state)

    def run(command, stage, capture=False):
        state['stage'] = stage
        save()
        command = list(map(str, command))
        record = dict(stage=stage, command=command, start=now())
        commands.append(record)
        write_json(logdir / 'commands.json', commands)
        print(f'[{sample["sample"]}] {stage}', flush=True)
        with (logdir / 'pipeline.log').open('a', encoding='utf-8') as log:
            log.write('\n[' + now() + '] ' + stage + '\n' + repr(command) + '\n')
            log.flush()
            result = subprocess.run(command, cwd=attempt, stdout=subprocess.PIPE if capture else log,
                                    stderr=log, text=True, errors='replace')
            record.update(end=now(), returncode=result.returncode)
            write_json(logdir / 'commands.json', commands)
            if result.returncode:
                raise RuntimeError(f'{stage} 失败（退出码 {result.returncode}），详见 logs/pipeline.log。')
            return result.stdout.strip() if capture else None

    old_handlers = {}
    def stopped(signum, frame):
        raise RuntimeError(f'收到停止信号 {signum}，请检查 Slurm 时限、内存和取消记录。')

    for sig in [signal.SIGTERM, signal.SIGINT]:
        old_handlers[sig] = signal.signal(sig, stopped)
    save()
    try:
        validate_trim_paths([sample], attempt)
        write_json(attempt / 'versions.json', dependencies())
        if fingerprint(sample['fastq']) != sample['fastq_fingerprint']:
            raise ValueError('FASTQ 大小或修改时间与提交时不同，请用新批次重新提交。')
        if hashlib.sha256(Path(sample['reference']).read_bytes()).hexdigest() != sample['reference_sha256']:
            raise ValueError('参考序列在提交后发生变化，请用新批次重新提交。')
        cpus = plan['resources']['cpus']
        if int(os.environ.get('SLURM_CPUS_PER_TASK', cpus)) < cpus:
            raise ValueError('Slurm 实际分配 CPU 少于本次设置。')
        for subdir in ['qc/raw', 'qc/clean', 'trimmed', 'ref']:
            (attempt / subdir).mkdir(parents=True, exist_ok=True)
        qc = state['qc']
        state['stage'] = '原始 FASTQ 完整性与计数'
        save()
        qc['raw'] = fastq_stats(sample['fastq'])
        save()
        run(['fastqc', '--threads', '1', '--outdir', attempt / 'qc/raw', sample['fastq']], '原始 FastQC')
        qc['raw_fastqc'] = fastqc_summary(attempt / 'qc/raw')
        run(['trim_galore', '--illumina', '--quality', '20', '--length', '25', '--stringency', '5',
             '--error_rate', '0.1', '--cores', '1', '--dont_gzip', '--output_dir', attempt / 'trimmed',
             sample['fastq']], 'Trim Galore 去接头与低质量末端')
        cleaned = list((attempt / 'trimmed').glob('*_trimmed.fq'))
        if len(cleaned) != 1:
            raise ValueError('没有找到唯一的清洗后 *_trimmed.fq 文件。')
        clean = cleaned[0]
        state['stage'] = '清洗后 FASTQ 完整性与计数'
        save()
        qc['clean'] = fastq_stats(clean)
        if qc['clean']['reads'] > qc['raw']['reads']:
            raise ValueError('清洗后 reads 数大于原始数。')
        qc['retention_pct'] = 100 * qc['clean']['reads'] / qc['raw']['reads']
        run(['fastqc', '--threads', '1', '--outdir', attempt / 'qc/clean', clean], '清洗后 FastQC')
        qc['clean_fastqc'] = fastqc_summary(attempt / 'qc/clean')
        sequences = list(sample['records'].values())
        if len({seq.replace('A', 'G') for seq in sequences}) != len(sequences):
            qc['warnings'].append('多个参考经 A→G 转换后完全相同；多重比对 reads 将被排除。')
        if any('N' in seq for seq in sequences):
            qc['warnings'].append('参考含 N；可能降低比对率。')
        reference = attempt / 'ref/amplicons.fa'
        with reference.open('w', encoding='ascii', newline='\n') as handle:
            for name, sequence in sample['records'].items():
                handle.write(f'>{name}\n{sequence}\n')
        core = Path(__file__).resolve().parent / 'core'
        prefix = sample['sample']
        run([sys.executable, core / 'prepare_amplicon_ref.py', '-f', reference, '-pre', prefix,
             '-o', attempt / 'ref', '-p', cpus], '建立参考和索引')
        run([sys.executable, core / 'run_GLORI_amplicon.py', '-i', core, '-q', clean,
             '-f', attempt / f'ref/{prefix}.AG_conversion.fa', '-f2', reference,
             '-b', attempt / f'ref/{prefix}.baseanno', '-pre', prefix, '-o', attempt / 'output',
             '-T', cpus, '--sort-threads', min(2, max(1, cpus - 1)), '--sort-memory', '256M',
             '-m', '2', '-M', '10000', '--cutoff', '3', '-c', '15', '-C', '5',
             '-r', '0.1', '-p', '0.005', '-adp', '0.005', '-s', '0.8', '-R', '0.8',
             '--method', 'binomial', '--keep-tmp'], 'GLORI 全流程')
        bam = attempt / f'output/{prefix}_r.sorted.bam'
        run(['samtools', 'quickcheck', bam], '检查 BAM 完整性')
        if not Path(str(bam) + '.bai').is_file():
            raise ValueError('BAM 索引未生成。')
        qc['mapped_reads'] = int(run(['samtools', 'view', '-c', '-F', '4', bam], '统计比对 reads', True))
        qc['reverse_reads'] = int(run(['samtools', 'view', '-c', '-f', '16', '-F', '4', bam], '检查反向比对', True))
        if not 0 < qc['mapped_reads'] <= qc['clean']['reads'] or qc['reverse_reads']:
            raise ValueError('比对计数异常、零比对或包含不支持的反向比对。')
        qc['mapping_pct'] = 100 * qc['mapped_reads'] / qc['clean']['reads']
        if qc['mapping_pct'] < 50:
            qc['warnings'].append('比对率低于 50%，请结合参考和接头检查；此值是提示阈值。')
        if qc['retention_pct'] < 50:
            qc['warnings'].append('清洗保留率低于 50%，请查看 Trim Galore 报告。')
        qc.update(export_sites(sample, attempt / 'output', attempt))
        save()
        if not plan['resources']['keep_all']:
            try:
                state['removed_temporary_files'] = cleanup_disposable(attempt, prefix)
            except OSError as error:
                qc['warnings'].append(f'结果有效，但临时文件清理不完整: {error}')
        state.update(state='SUCCESS', stage='完成', end=now())
        save()
        return 0
    except Exception as error:
        state.update(state='FAILED', error=f'{state["stage"]}: {error}', end=now())
        save()
        with (logdir / 'pipeline.log').open('a', encoding='utf-8') as log:
            traceback.print_exc(file=log)
        print(f'[{sample["sample"]}] FAILED: {error}', file=sys.stderr, flush=True)
        return 1
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
