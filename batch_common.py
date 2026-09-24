"""Shared, standard-library-only input validation and durable file helpers."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
from datetime import datetime, timezone

VERSION = '1.0.1'
DEFAULTS = dict(cpus=12, max_parallel=4, memory='16G', time='08:00:00',
                partition=None, account=None, keep_all=False)
BLAS_ENV = {key: '1' for key in ['OPENBLAS_NUM_THREADS', 'OPENBLAS_DEFAULT_NUM_THREADS',
            'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']}
ANALYSIS = dict(mismatch=2, max_depth=10000, cutoff='3', coverage=15, count=5,
                ratio=0.1, pvalue=0.005, fdr=0.005, signal=0.8, var_ratio=0.8,
                method='binomial', sort_threads=2, sort_memory='256M',
                trim_quality=20, trim_length=25, trim_overlap=5, trim_error=0.1)


def limit_math_threads():
    """Call before scientific imports, including login-node dependency checks."""
    os.environ.update(BLAS_ENV)


def validate_trim_paths(samples, output=None):
    """Trim Galore 0.6.10 rejects whitespace in FASTQ paths; fail before submission."""
    paths = [sample['fastq'] for sample in samples]
    if output is not None:
        paths.append(str(Path(output).resolve()))
    for path in paths:
        if any(character.isspace() for character in str(path)):
            raise ValueError('Trim Galore requires FASTQ and output paths without whitespace: '
                             + str(path))


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(temporary, path)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_tsv(path, fields, rows):
    path = Path(path)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with temporary.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def identifier(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value) or '_AG_converted' in value:
        raise ValueError(f'名称不合法: {value!r}；请以字母或数字开头，仅用字母、数字、_、.、-，勿用 _AG_converted。')
    return value


def read_fasta(path):
    records = {}
    name = None
    with Path(path).open(encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                name = identifier(line[1:].split()[0] if line[1:].split() else '')
                if name in {'ALL', 'ELSE', 'Median', 'Mean', 'discard', 'control'} or name.startswith('Median_'):
                    raise ValueError(f'{path}: {name} 是 GLORI 统计保留名，请更换 FASTA 名称。')
                if name in records:
                    raise ValueError(f'{path}: FASTA 名称重复: {name}')
                records[name] = ''
            elif name is None:
                raise ValueError(f'{path}: FASTA 必须以 >名称 开始。')
            else:
                records[name] += line.upper()
    if not records:
        raise ValueError(f'{path}: 空 FASTA。')
    for name, sequence in records.items():
        if not re.fullmatch('[ACGTN]+', sequence) or 'A' not in sequence:
            raise ValueError(f'{path}: {name} 必须是含 A 的原始 DNA 序列，仅允许 A/C/G/T/N。')
    return records


def fastq_records(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as handle:
        number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip('\r\n')
            plus = handle.readline()
            quality = handle.readline().rstrip('\r\n')
            number += 1
            if (not header.startswith('@') or not plus.startswith('+') or not sequence
                    or len(sequence) != len(quality) or re.search('[^ACGTNacgtn]', sequence)
                    or min(quality, default='!') < '!' or max(quality, default='!') > '~'):
                raise ValueError(f'{path}: 第 {number} 条 FASTQ 不完整或格式错误。')
            yield sequence, quality


def fastq_stats(path):
    count = bases = q20 = q30 = 0
    shortest = None
    longest = 0
    # Count qualities in C rather than a Python per-base loop over millions of reads.
    quality_counts = [0] * 94
    from collections import Counter
    histogram = Counter()
    for sequence, quality in fastq_records(path):
        length = len(sequence)
        count += 1
        bases += length
        shortest = length if shortest is None else min(shortest, length)
        longest = max(longest, length)
        histogram.update(quality)
    if not count:
        raise ValueError(f'{path}: 没有可用 reads。')
    for quality, number in histogram.items():
        quality_counts[ord(quality) - 33] += number
    q20, q30 = sum(quality_counts[20:]), sum(quality_counts[30:])
    return dict(reads=count, bases=bases, min_length=shortest, max_length=longest,
                mean_length=bases / count, q20_pct=100 * q20 / bases, q30_pct=100 * q30 / bases)


def fingerprint(path):
    stat = Path(path).stat()
    return dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def load_samples(path):
    path = Path(path).resolve()
    samples = []
    names = set()
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['sample', 'fastq', 'reference']:
            raise ValueError('样品表表头必须为 sample、fastq、reference，以 Tab 分隔。')
        for line, row in enumerate(reader, 2):
            if None in row or any(not v or not v.strip() for v in row.values()):
                raise ValueError(f'样品表第 {line} 行必须有三列非空内容。')
            name = identifier(row['sample'].strip())
            if name.casefold() in names:
                raise ValueError(f'样品名重复（忽略大小写）: {name}')
            names.add(name.casefold())
            sample = dict(sample=name)
            for key in ['fastq', 'reference']:
                raw = row[key].strip()
                resolved = (path.parent / raw).resolve()
                if not resolved.is_file() or resolved.stat().st_size == 0:
                    raise ValueError(f'{name}: {key} 文件不存在或为空: {resolved}')
                sample[key] = str(resolved)
            if not sample['fastq'].lower().endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
                raise ValueError(f'{name}: FASTQ 后缀须为 .fq/.fastq，可附加 .gz。')
            # Only probe one read on login node; full validation runs in the sample job.
            probe = fastq_records(sample['fastq'])
            try:
                next(probe)
            except StopIteration:
                raise ValueError(f'{name}: FASTQ 中没有 reads。')
            finally:
                probe.close()
            sample['records'] = read_fasta(sample['reference'])
            sample['fastq_fingerprint'] = fingerprint(sample['fastq'])
            sample['reference_sha256'] = hashlib.sha256(Path(sample['reference']).read_bytes()).hexdigest()
            samples.append(sample)
    if not samples:
        raise ValueError('样品表没有样品。')
    return samples


class BatchLock:
    """OS lock is released on process exit, including crashes; never delete lock inode."""
    def __init__(self, root):
        self.path = Path(root) / '.batch.lock'

    def __enter__(self):
        self.handle = self.path.open('a+b')
        self.handle.seek(0)
        if self.path.stat().st_size == 0:
            self.handle.write(b'0')
            self.handle.flush()
            self.handle.seek(0)
        try:
            if os.name == 'nt':
                import msvcrt
                msvcrt.locking(self.handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(self.handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self.handle.close()
            raise ValueError('另一进程正在提交或汇总此批次，请稍后再试。')
        return self

    def __exit__(self, *args):
        self.handle.close()
