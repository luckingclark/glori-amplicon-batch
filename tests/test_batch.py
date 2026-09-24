"""Local regressions use synthetic data and a mocked scheduler, never real submissions."""
import argparse
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import batch_common as common
import batch_report as report
import batch_worker as worker
import run_batch as batch


def fixture(root, name='S1', records=None):
    root = Path(root)
    reference = root / f'{name}.fa'
    records = records or {'AMP1': 'ACGA'}
    reference.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in records.items()))
    fastq = root / f'{name}.fq.gz'
    with gzip.open(fastq, 'wt') as handle:
        handle.write('@r1\nACGT\n+\nIIII\n')
    sheet = root / f'{name}.tsv'
    sheet.write_text(f'sample\tfastq\treference\n{name}\t{fastq.name}\t{reference.name}\n')
    return sheet, common.load_samples(sheet)[0]


def format_line(name='AMP1', position=1, total=100, ag=100, a=25, used=(80, 80, 20)):
    return '\t'.join([name, str(position), '+', name, name, name, 'NA', 'amplicon',
                     str(total), str(ag), str(a), '0', str(total-ag), str(ag-a)]
                     + [str(c) + ';' + ','.join(map(str, used)) for c in [1,2,3,4,5,6,7,8,9,10,15,20]]) + '\n'


def output_fixture(attempt, sample, lines=None, significant=False):
    output = attempt / 'output'
    (output / 'tmp').mkdir(parents=True)
    name = sample['sample']
    (output / 'tmp' / f'{name}.totalformat.txt').write_text(format_line() if lines is None else lines)
    (output / f'{name}.totalCR.txt').write_text('SA\tA-to-G_ratio\n#Median_AMP1\t0.95\n#AMP1\t0.95\n')
    rows = []
    if significant:
        rows = [dict(zip(report.FDR_FIELDS, ['AMP1', 1, '+', 'AMP1', 0.95, 80, 20, 100, 0.25, 0.0001, 0.0001]))]
    common.write_tsv(output / f'{name}.totalm6A.FDR.csv', report.FDR_FIELDS, rows)
    return output


class InputTests(unittest.TestCase):
    def test_relative_absolute_bom_and_multiple_contigs(self):
        with tempfile.TemporaryDirectory(prefix='glori spaces ') as temp:
            root = Path(temp)
            sheet, sample = fixture(root, records={'A1': 'ACGT', 'A2': 'AACN'})
            self.assertEqual(len(sample['records']), 2)
            sheet.write_text(f'\ufeffsample\tfastq\treference\nS1\t{sample["fastq"]}\t{sample["reference"]}\n', encoding='utf-8')
            self.assertEqual(common.load_samples(sheet)[0], sample)

    def test_duplicate_samples_missing_files_bad_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sheet, sample = fixture(root)
            original = sheet.read_text()
            sheet.write_text(original + 's1\tS1.fq.gz\tS1.fa\n')
            with self.assertRaises(ValueError):
                common.load_samples(sheet)
            sheet.write_text(original.replace('S1.fq.gz', 'missing.fq.gz'))
            with self.assertRaises(ValueError):
                common.load_samples(sheet)
            sheet.write_text(original)
            Path(sample['reference']).write_text('>AMP1\nGGGT\n')
            with self.assertRaises(ValueError):
                common.load_samples(sheet)

    def test_fastq_counts_quality_and_truncation(self):
        with tempfile.TemporaryDirectory() as temp:
            file = Path(temp) / 'data.fq'
            file.write_text('@a\nACGT\n+\nI5!I\n@b\nAA\n+\nII\n')
            stats = common.fastq_stats(file)
            self.assertEqual((stats['reads'], stats['bases'], stats['min_length']), (2, 6, 2))
            self.assertAlmostEqual(stats['q30_pct'], 100 * 4 / 6)
            file.write_text('@a\nACGT\n+\nIII\n')
            with self.assertRaises(ValueError):
                common.fastq_stats(file)


class CountTests(unittest.TestCase):
    def test_ratios_before_after_cutoff_significance_and_absent_is_na(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _, sample = fixture(root)
            output = output_fixture(root, sample, significant=True)
            qc = report.export_sites(sample, output, root)
            rows = report.read_table(root / 'all_A_sites.tsv', report.SITE_FIELDS)
            self.assertEqual(qc['significant_sites'], 1)
            self.assertEqual(rows[0]['Ratio'], '0.25')
            self.assertEqual(rows[0]['Percent'], '25.0')
            self.assertEqual(rows[0]['Gcov'], '60')
            self.assertEqual(rows[0]['FDR_status'], 'significant')
            self.assertEqual(rows[1]['Sites'], '4')
            self.assertEqual(rows[1]['Acov'], 'NA')
            self.assertEqual(rows[1]['Observation_status'], 'no_usable_record')

    def test_zero_denominator_and_empty_significant_table(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _, sample = fixture(root)
            output = output_fixture(root, sample, lines=format_line(a=5, used=(0,0,0)))
            report.export_sites(sample, output, root)
            row = report.read_table(root / 'all_A_sites.tsv', report.SITE_FIELDS)[0]
            self.assertEqual(row['Ratio'], 'NA')
            self.assertEqual(row['Acov'], '0')
            self.assertEqual(row['Observation_status'], 'no_AG_after_cutoff')
            self.assertEqual(report.read_table(root / 'significant_sites.tsv', ['Sample'] + report.FDR_FIELDS), [])

    def test_malformed_counts_and_missing_fdr_are_errors(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _, sample = fixture(root)
            output = output_fixture(root, sample, lines=format_line(used=(100,100,101)))
            with self.assertRaises(ValueError):
                report.export_sites(sample, output, root)
            (output / 'S1.totalm6A.FDR.csv').unlink()
            with self.assertRaises(FileNotFoundError):
                report.export_sites(sample, output, root)

    def test_failed_interrupted_samples_do_not_leak_partial_ratios(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _, sample = fixture(root)
            attempt = root / 'samples/S1/attempt_001'
            attempt.mkdir(parents=True)
            output_fixture(attempt, sample, significant=True)
            report.export_sites(sample, attempt / 'output', attempt)
            common.write_json(attempt / 'status.json', dict(state='RUNNING', qc={}))
            data = dict(samples=[sample], latest={'S1': 'samples/S1/attempt_001'}, resources=common.DEFAULTS)
            failures = report.collect(root, data, {'S1': 'TIMEOUT / 0:15'})
            self.assertEqual(failures, 1)
            rows = report.read_table(root / 'all_A_sites.tsv', report.SITE_FIELDS)
            self.assertTrue(all(row['Ratio'] == 'NA' for row in rows))
            self.assertIn('TIMEOUT', (root / 'qc_summary.tsv').read_text(encoding='utf-8'))
            self.assertEqual(report.read_table(root / 'significant_sites.tsv', ['Sample'] + report.FDR_FIELDS), [])

    def test_missing_success_file_is_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _, sample = fixture(root)
            attempt = root / 'a'
            attempt.mkdir()
            common.write_json(attempt / 'status.json', dict(state='SUCCESS'))
            data = dict(samples=[sample], latest={'S1': 'a'}, resources=common.DEFAULTS)
            self.assertEqual(report.collect(root, data), 1)

    def test_cleanup_keeps_counts_and_sibling_files(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            directory = root / 'output/tmp'
            directory.mkdir(parents=True)
            names = ['S1_AGchanged_2.fq', 'S1.referbase.mpi.AMP1', 'S1.totalformat.txt',
                     'S1.referbase.mpi.formatted.txt.AMP1', 'S2_AGchanged_2.fq']
            for name in names:
                (directory / name).write_text('data')
            worker.cleanup_disposable(root, 'S1')
            self.assertFalse((directory / names[0]).exists())
            self.assertFalse((directory / names[1]).exists())
            self.assertTrue(all((directory / name).exists() for name in names[2:]))

    def test_fastqc_warn_is_not_execution_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'x_fastqc.html').write_text('report')
            with zipfile.ZipFile(root / 'x_fastqc.zip', 'w') as archive:
                archive.writestr('x_fastqc/summary.txt', 'PASS\tBasic Statistics\tx.fq\nFAIL\tSequence Duplication Levels\tx.fq\n')
            self.assertIn('FAIL:Sequence Duplication Levels', worker.fastqc_summary(root))


class SchedulerTests(unittest.TestCase):
    def test_whitespace_paths_fail_before_dependency_check_or_submission(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sheet, _ = fixture(root)
            args = self.arguments(root, sheet)
            args.out = str(root / 'results with spaces')
            with patch.object(worker, 'dependencies') as dependencies, \
                    patch.object(batch, 'scheduler') as scheduler:
                with self.assertRaisesRegex(ValueError, 'without whitespace'):
                    batch.submit(args)
                args.out = str(root / 'results')
                moved = root / 'reads with spaces.fq.gz'
                (root / 'S1.fq.gz').rename(moved)
                sheet.write_text('sample\tfastq\treference\nS1\treads with spaces.fq.gz\tS1.fa\n')
                with self.assertRaisesRegex(ValueError, 'without whitespace'):
                    batch.submit(args)
                dependencies.assert_not_called()
                scheduler.assert_not_called()
                self.assertFalse(Path(args.out).exists())

    def test_new_batch_requires_explicit_cpu_partition(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sheet, _ = fixture(root)
            args = self.arguments(root, sheet)
            args.partition = None
            with patch.object(worker, 'dependencies', return_value={}), \
                    patch.dict(os.environ, {'CONDA_PREFIX': sys.prefix}), \
                    patch.object(batch, 'scheduler') as scheduler:
                with self.assertRaisesRegex(ValueError, '--partition'):
                    batch.submit(args)
                scheduler.assert_not_called()
                self.assertFalse((root / 'results').exists())

    def arguments(self, root, sheet, retry=False):
        return argparse.Namespace(samples=str(sheet), out=str(root / 'results'), retry_failed=retry,
            validate_inputs=False, check_only=False, cpus=None, max_parallel=None, memory=None,
            time=None, partition='test_cpu', account=None, keep_all=None)

    def test_resources_dependency_retry_and_output_safety(self):
        with tempfile.TemporaryDirectory(prefix='batch_') as temp:
            root = Path(temp)
            sheet, sample = fixture(root)
            _, sample2 = fixture(root, name='S2')
            with sheet.open('a') as handle:
                handle.write('S2\tS2.fq.gz\tS2.fa\n')
            calls = []
            def fake_scheduler(command):
                calls.append(command)
                if command[0] == 'sbatch':
                    return str(100 + sum(c[0] == 'sbatch' for c in calls))
                return ''
            with patch.object(worker, 'dependencies', return_value={}), patch.dict(os.environ, {'CONDA_PREFIX': sys.prefix}), \
                    patch.object(batch, 'scheduler', side_effect=fake_scheduler), patch.object(os, 'getuid', return_value=1, create=True):
                self.assertEqual(batch.submit(self.arguments(root, sheet)), 0)
                array = calls[0]
                self.assertIn('--array=0-1%4', array)
                self.assertIn('--cpus-per-task=12', array)
                self.assertIn('--mem=16G', array)
                self.assertIn('--partition=test_cpu', array)
                self.assertIn('--partition=test_cpu', calls[1])
                self.assertIn('--hold', array)
                self.assertIn('--dependency=afterany:101', calls[1])
                self.assertEqual(calls[2], ['scontrol', 'release', '101'])
                result = root / 'results'
                self.assertEqual((result / 'software/LICENSES/GLORI-tools-MIT.txt').read_bytes(),
                                 (ROOT / 'LICENSES/GLORI-tools-MIT.txt').read_bytes())
                state = common.read_json(result / 'batch.json')
                attempt = result / state['latest']['S1']
                output_fixture(attempt, sample)
                report.export_sites(sample, attempt / 'output', attempt)
                common.write_json(attempt / 'status.json', dict(state='SUCCESS'))
                # Old batches used an unspecified partition; retry remains compatible.
                state['version'] = '1.0.0'
                state['resources']['partition'] = None
                common.write_json(result / 'batch.json', state)
                self.assertEqual(batch.submit(self.arguments(root, sheet, retry=True)), 0)
                newstate = common.read_json(result / 'batch.json')
                self.assertEqual(newstate['resources']['partition'], 'test_cpu')
                self.assertEqual(newstate['latest']['S1'], state['latest']['S1'])
                self.assertNotEqual(newstate['latest']['S2'], state['latest']['S2'])
                self.assertEqual(newstate['rounds'][1]['tasks'], ['S2'])
                self.assertTrue((result / state['latest']['S2']).exists())
                with self.assertRaises(ValueError):
                    batch.submit(self.arguments(root, sheet))

    def test_partition_override_and_legacy_retry_default(self):
        args = self.arguments(Path('/example'), Path('/samples.tsv'))
        legacy = dict(common.DEFAULTS, partition=None)
        self.assertEqual(batch.resource_options(args, legacy)['partition'], 'test_cpu')
        args.partition = 'another_cpu_partition'
        self.assertEqual(batch.resource_options(args)['partition'], 'another_cpu_partition')

    def test_math_threads_limited_before_first_scientific_import(self):
        class StopBeforeImport(Exception):
            pass
        def inspect_import(name):
            self.assertEqual(name, 'numpy')
            self.assertTrue(all(os.environ[key] == '1' for key in common.BLAS_ENV))
            raise StopBeforeImport
        with patch.dict(os.environ, {key: '48' for key in common.BLAS_ENV}), \
                patch.object(worker.sys, 'platform', 'linux'), \
                patch.object(worker.importlib, 'import_module', side_effect=inspect_import):
            with self.assertRaises(StopBeforeImport):
                worker.dependencies(slurm=True)

    def test_active_array_child_prevents_retry(self):
        data = dict(rounds=[dict(array_job='101', collect_job='102')])
        with patch.object(batch, 'scheduler', return_value='101\n') as scheduler, \
                patch.object(os, 'getuid', return_value=1, create=True):
            with self.assertRaises(ValueError):
                batch.assert_idle(data)
            self.assertIn('--format=%F', scheduler.call_args.args[0])

    def test_accounting_maps_array_task_states_not_raw_job_ids(self):
        data = dict(rounds=[dict(array_job='101', tasks=['S1', 'S2'])])
        with patch.object(batch.shutil, 'which', return_value='/bin/sacct'), \
                patch.object(batch, 'scheduler', return_value='101_0|COMPLETED|0:0\n101_1|OUT_OF_MEMORY|0:9\n101_1.batch|FAILED|0:9\n') as scheduler:
            states = batch.scheduler_states(data)
            self.assertIn('OUT_OF_MEMORY', states['S2'])
            self.assertIn('--format=JobID%100,State%40,ExitCode', scheduler.call_args.args[0])

    def test_collection_submission_failure_cancels_held_array(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sheet, _ = fixture(root)
            calls = []
            def fake(command):
                calls.append(command)
                if command[0] == 'sbatch':
                    if len(calls) == 1:
                        return '123'
                    raise RuntimeError('invalid partition')
                return ''
            with patch.object(worker, 'dependencies', return_value={}), patch.dict(os.environ, {'CONDA_PREFIX': sys.prefix}), \
                    patch.object(batch, 'scheduler', side_effect=fake):
                with self.assertRaises(RuntimeError):
                    batch.submit(self.arguments(root, sheet))
            self.assertIn(['scancel', '123'], calls)
            self.assertFalse(any(c[:2] == ['scontrol', 'release'] for c in calls))
            self.assertIn('提交失败', (root / 'results/qc_summary.tsv').read_text(encoding='utf-8'))

    def test_bash_script_quotes_paths(self):
        text = batch.script_text('/a env/bin/python', Path('/x space/run_batch.py'), Path('/x space/r.json'),
                                 '--worker', {'PATH': "/path with ' quote"})
        self.assertIn('"${SLURM_ARRAY_TASK_ID:?Missing array task index}"', text)
        bash = shutil.which('bash') if sys.platform == 'linux' else None
        if Path(bash or '').is_file():
            result = subprocess.run([bash, '--noprofile', '--norc', '-n'], input=text, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)


class WorkerTests(unittest.TestCase):
    def test_sample_failure_isolated_and_success_keeps_key_outputs(self):
        """Exercise real wrapper state transitions with fake external-tool outputs."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            tasks = []
            for name in ['BAD', 'GOOD']:
                _, sample = fixture(root, name=name)
                with gzip.open(sample['fastq'], 'wt') as handle:
                    handle.write(''.join(f'@r{i}\nACGA\n+\nIIII\n' for i in range(100)))
                sample['fastq_fingerprint'] = common.fingerprint(sample['fastq'])
                sample['attempt'] = name
                (root / name).mkdir()
                tasks.append(sample)
            plan = root / 'round.json'
            common.write_json(plan, dict(root=str(root), tasks=tasks, resources=common.DEFAULTS))
            def fake_run(command, **kwargs):
                cwd = Path(kwargs['cwd'])
                if command[0] == 'trim_galore':
                    if cwd.name == 'BAD':
                        return subprocess.CompletedProcess(command, 7)
                    with gzip.open(command[-1], 'rt') as source:
                        (cwd / 'trimmed/good_trimmed.fq').write_text(source.read())
                elif command[0] == 'fastqc':
                    output = Path(command[command.index('--outdir')+1])
                    (output / 'x_fastqc.html').write_text('report')
                    with zipfile.ZipFile(output / 'x_fastqc.zip', 'w') as archive:
                        archive.writestr('x_fastqc/summary.txt', 'WARN\tPer base sequence content\tx.fq\n')
                elif command[0] == sys.executable and command[1].endswith('run_GLORI_amplicon.py'):
                    output_fixture(cwd, tasks[1], significant=True)
                    (cwd / 'output/GOOD_r.sorted.bam').write_bytes(b'fake-bam')
                    (cwd / 'output/GOOD_r.sorted.bam.bai').write_bytes(b'fake-index')
                    (cwd / 'output/tmp/GOOD_AGchanged_2.fq').write_text('disposable')
                if command[:2] == ['samtools', 'view']:
                    return subprocess.CompletedProcess(command, 0, stdout='0' if '-f' in command else '100')
                return subprocess.CompletedProcess(command, 0, stdout='')
            with patch.object(worker, 'dependencies', return_value={}), \
                    patch.object(worker.subprocess, 'run', side_effect=fake_run):
                self.assertEqual(worker.work(plan, 0), 1)
                self.assertEqual(worker.work(plan, 1), 0)
            bad, good = common.read_json(root / 'BAD/status.json'), common.read_json(root / 'GOOD/status.json')
            self.assertEqual(bad['state'], 'FAILED')
            self.assertEqual(good['state'], 'SUCCESS')
            self.assertIn('7', bad['error'])
            self.assertEqual(good['qc']['mapping_pct'], 100)
            self.assertFalse((root / 'GOOD/output/tmp/GOOD_AGchanged_2.fq').exists())
            self.assertTrue((root / 'GOOD/output/tmp/GOOD.totalformat.txt').exists())
            data = dict(samples=tasks, latest={'BAD': 'BAD', 'GOOD': 'GOOD'}, resources=common.DEFAULTS)
            self.assertEqual(report.collect(root, data), 1)
            rows = report.read_table(root / 'all_A_sites.tsv', report.SITE_FIELDS)
            self.assertTrue(all(row['Ratio'] == 'NA' for row in rows if row['Sample'] == 'BAD'))
            self.assertTrue(any(row['Ratio'] == '0.25' for row in rows if row['Sample'] == 'GOOD'))

    def test_changed_input_is_rejected_before_external_analysis(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _, sample = fixture(root)
            sample['attempt'] = 'attempt'
            (root / 'attempt').mkdir()
            plan = root / 'round.json'
            common.write_json(plan, dict(root=str(root), tasks=[sample], resources=common.DEFAULTS))
            Path(sample['reference']).write_text('>AMP1\nAAAT\n')
            with patch.object(worker, 'dependencies', return_value={}):
                self.assertEqual(worker.work(plan, 0), 1)
            state = common.read_json(root / 'attempt/status.json')
            self.assertEqual(state['state'], 'FAILED')
            self.assertIn('参考序列', state['error'])


RUN_LINUX = sys.platform == 'linux' and os.environ.get('RUN_GLORI_INTEGRATION') == '1'


@unittest.skipUnless(RUN_LINUX, 'requires Linux, full Conda dependencies and RUN_GLORI_INTEGRATION=1')
class LinuxIntegrationTests(unittest.TestCase):
    def test_worker_and_direct_core_match(self):
        """Real Trim Galore/FastQC/Bowtie/samtools/GLORI on synthetic reads; no Slurm."""
        worker.dependencies()
        with tempfile.TemporaryDirectory(prefix='glori_full_') as temp:
            root = Path(temp)
            rng = random.Random(724)
            sequence = ''.join(rng.choice('ACGT') for _ in range(80))
            target = sequence.index('A')
            reference = root / 'original.fa'
            reference.write_text('>AMP1\n' + sequence + '\n>UNOBSERVED\n' + 'AC' * 30 + '\n')
            fastq = root / 'input.fq'
            with fastq.open('w') as handle:
                for i in range(200):
                    read = list(sequence.replace('A', 'G'))
                    if i < 80:
                        read[target] = 'A'
                    read = ''.join(read) + 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
                    handle.write(f'@r{i}\n{read}\n+\n' + 'I' * len(read) + '\n')
            sheet = root / 'samples.tsv'
            sheet.write_text('sample\tfastq\treference\nS1\tinput.fq\toriginal.fa\n')
            sample = common.load_samples(sheet)[0]
            sample['attempt'] = 'attempt'
            (root / 'attempt').mkdir()
            plan = root / 'round.json'
            resources = dict(common.DEFAULTS, cpus=3, keep_all=True)
            common.write_json(plan, dict(root=str(root), tasks=[sample], resources=resources))
            result = subprocess.run([sys.executable, str(ROOT / 'run_batch.py'), '--worker', str(plan), '--index', '0'],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr + (root / 'attempt/logs/pipeline.log').read_text())
            rows = report.read_table(root / 'attempt/all_A_sites.tsv', report.SITE_FIELDS)
            target_row = next(row for row in rows if row['Chr'] == 'AMP1' and int(row['Sites']) == target + 1)
            self.assertEqual(float(target_row['Ratio']), 0.4)
            self.assertEqual(target_row['FDR_status'], 'significant')
            self.assertTrue(all(row['Ratio'] == 'NA' for row in rows if row['Chr'] == 'UNOBSERVED'))
            clean = next((root / 'attempt/trimmed').glob('*_trimmed.fq'))
            result = subprocess.run([sys.executable, str(ROOT / 'core/run_GLORI_amplicon.py'), '-q', str(clean),
                '-f', str(root / 'attempt/ref/S1.AG_conversion.fa'), '-f2', str(root / 'attempt/ref/amplicons.fa'),
                '-b', str(root / 'attempt/ref/S1.baseanno'), '-pre', 'S1', '-o', str(root / 'direct'),
                '-T', '3', '--keep-tmp'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            for name in ['S1.totalm6A.FDR.csv', 'tmp/S1.totalformat.txt', 'S1.totalCR.txt']:
                self.assertEqual((root / 'direct' / name).read_bytes(), (root / 'attempt/output' / name).read_bytes())


if __name__ == '__main__':
    unittest.main()
