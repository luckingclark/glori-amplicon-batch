#!/usr/bin/env python
"""
GLORI amplicon deep sequencing analysis pipeline.
Simplified from GLORI-tools v1.0 for targeted amplicon analysis.

Workflow:
  1. A-to-G conversion of reads
  2. Bowtie alignment to A-to-G converted amplicon reference
  3. Reverse reads back to original sequence
  4. Convert SAM to sorted BAM
  5. Pileup
  6. Reference base correction (restore original FASTA base)
  7. Per-amplicon formatting + m6A calling (binomial test vs background CR)
  8. Merge results + compute CR + FDR correction

Input:
  - Cleaned single-end FASTQ (adapter trimmed; no molecular UMI in this workflow)
  - Amplicon reference prepared by prepare_amplicon_ref.py

Usage:
  python run_GLORI_amplicon.py \
      -i /path/to/your/glori-amplicon \
      -q cleaned_reads.fq \
      -f amplicon.AG_conversion.fa \
      -f2 amplicon.fa \
      -b amplicon.baseanno \
      -pre sample1 \
      -o ./output \
      -T 4

Requirements:
  - bowtie >= 1.3.1, samtools >= 1.10
  - pysam, BioPython, pandas, numpy, scipy, statsmodels
"""

import os
import sys
import re
import argparse
import subprocess
import itertools
import glob
import gzip
import multiprocessing
import time
import shutil
import shlex
import importlib
from heapq import merge
from collections import defaultdict

import pandas as pd
from Bio.Seq import reverse_complement
from Bio import SeqIO


def run_command(command, **kwargs):
    """Run an argument list and stop immediately when a subprocess fails."""
    print('  Running: ' + shlex.join([str(arg) for arg in command]), flush=True)
    return subprocess.run([str(arg) for arg in command], check=True, **kwargs)


def safe_identifier(value):
    """Identifiers become filenames and annotation keys in the upstream tools."""
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise argparse.ArgumentTypeError(
            'Use letters, numbers, underscore, dot or hyphen; start with a letter or number.')
    if '_AG_converted' in value:
        raise argparse.ArgumentTypeError('The suffix _AG_converted is reserved.')
    return value


def positive_integer(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('Must be at least 1.')
    return number


def prepare_output_directory(path):
    """Never silently resume a run with cached files from another analysis."""
    if os.path.exists(path) and (not os.path.isdir(path) or os.listdir(path)):
        raise ValueError(f'Output directory must be new or empty: {path}. Choose a new -o directory.')
    os.makedirs(path, exist_ok=True)


def validate_alignment_orientation(sam_path):
    """The inherited restoration algorithm only supports forward alignments."""
    mapped = 0
    with open(sam_path) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            fields = line.rstrip('\n').split('\t')
            flag = int(fields[1])
            if flag & 4:
                continue
            if flag & 16:
                raise ValueError(
                    'Reverse-strand alignment detected. This amplicon workflow only supports '
                    'reads aligned forward to the RNA-oriented reference. Check the reference '
                    'orientation and read end; the pipeline has stopped before restoring bases.')
            mapped += 1
    if mapped == 0:
        raise ValueError('No reads mapped. Check the input reads and amplicon reference.')


def merge_files(paths, destination, skip_comments=False):
    with open(destination, 'w') as output:
        for path in paths:
            with open(path) as source:
                for line in source:
                    if not (skip_comments and line.startswith('#')):
                        output.write(line)


# ============================================================
# Helper functions extracted from mapping_reads.py
# ============================================================

def embedded_numbers(s):
    """Extract numbers from string for natural sorting of chrom names."""
    s2 = s.strip().split('\t')
    pieces = re.split(r'(\d+)', s2[0])
    pieces[1::2] = list(map(int, pieces[1::2]))
    # Exact-name tie breaking keeps duplicate IDs adjacent even when different
    # names have equal numeric keys (e.g. read1 and read01).
    return pieces, s2[0]


def sort_bedfiles(bedfiles, outputfiles):
    """Sort BED file by chromosome name, then position."""
    prx2 = bedfiles[:-4]
    path = glob.escape(prx2) + '_chunk_*.bed'
    chunksize = 5000000
    fid = 1
    lines = []
    with open(bedfiles, 'r') as f_in:
        f_out = open(prx2 + '_chunk_{}.bed'.format(fid), 'w')
        for line_num, line in enumerate(f_in, 1):
            lines.append(line)
            if not line_num % chunksize:
                lines = sorted(lines, key=embedded_numbers)
                f_out.writelines(lines)
                f_out.close()
                lines = []
                fid += 1
                f_out = open(prx2 + '_chunk_{}.bed'.format(fid), 'w')
        if lines:
            lines = sorted(lines, key=embedded_numbers)
            f_out.writelines(lines)
            f_out.close()
            lines = []
        if not f_out.closed:
            f_out.close()

    chunks = []
    for filename in glob.glob(path):
        chunks += [open(filename, 'r')]

    with open(outputfiles, 'w') as f_out:
        f_out.writelines(merge(*chunks, key=embedded_numbers))
    for chunk in chunks:
        chunk.close()
    for filename in glob.glob(path):
        os.remove(filename)
    # Read restoration uses the first FASTQ header token as a unique key.
    # Scan the sorted BED rather than retaining millions of IDs in a set.
    previous_name = None
    with open(outputfiles) as sorted_bed:
        for row in sorted_bed:
            read_name = row.split('\t', 1)[0]
            if read_name == previous_name:
                raise ValueError(
                    f'Duplicate FASTQ read ID: {read_name}. The first header token must be '
                    'unique; do not combine read ends or files with overlapping read names.')
            previous_name = read_name


def change_reads(fastq, changename, output_bed, outputdir, change_fac, step=10000):
    """
    Convert all A -> G in reads, record original A positions in BED file.
    change_fac: e.g. 'AG' means replace A with G.
    """
    fac_t = change_fac[1]   # target base (G)
    fac_q = change_fac[0]   # query base (A)

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    file_change = open(changename, 'w')
    file_bed = open(output_bed, 'w')

    # Support both plain and gzipped FASTQ input
    if fastq.endswith('.gz'):
        fh = gzip.open(fastq, 'rt')
    else:
        fh = open(fastq, 'r')

    try:
        while True:
            fr = list(itertools.islice(fh, step * 4))
            if len(fr) == 0:
                break

            list_change = []
            list_bed = []
            list_fr = [line.strip() for line in fr]

            for x in range(0, len(list_fr), 4):
                if x + 3 >= len(list_fr):
                    raise ValueError('Incomplete FASTQ record in input.')
                reads_name = list_fr[x].split()[0]
                yr = list_fr[x + 1].upper()
                if (not reads_name.startswith('@') or not list_fr[x + 2].startswith('+')
                        or len(yr) != len(list_fr[x + 3])):
                    raise ValueError(f'Malformed FASTQ record: {reads_name}')
                A_sites = [m.start() for m in re.finditer('A', yr)]

                if len(A_sites) >= 1:
                    A_sites2 = '_'.join(map(str, A_sites))
                    list_bed.append([reads_name[1:], A_sites2])
                    list_change += [
                        reads_name,
                        yr.replace(fac_q, fac_t),
                        list_fr[x + 2],
                        list_fr[x + 3]
                    ]
                else:
                    list_bed.append([reads_name[1:], 'NA'])
                    list_change += [reads_name, yr, list_fr[x + 2], list_fr[x + 3]]

            file_change.write('\n'.join(list_change) + '\n')
            list_bed1 = ['\t'.join(map(str, it)) for it in list_bed]
            file_bed.write('\n'.join(list_bed1) + '\n')

        file_change.close()
        file_bed.close()
    finally:
        fh.close()
        file_change.close()
        file_bed.close()

    sort_bedfiles(output_bed, output_bed + '_sorted')
    os.remove(output_bed)


def multi_sub(string, sitesA, repl):
    """Replace bases at specific positions in a string."""
    if sitesA != 'NA':
        A_list = list(map(int, sitesA.split('_')))
        new = list(string)
        for index in A_list:
            new[index] = repl
        return ''.join(new)
    else:
        return string


def reverse_reads(output_sam, output_bed, reverse_fac, flag, step=10000,
                  sort_threads=2, sort_memory='256M'):
    """
    Reverse the A->G conversion in aligned reads back to original sequence.
    Uses the BED file to know which positions to restore.
    """
    validate_alignment_orientation(output_sam)
    sorted_sam = output_sam[:-4] + '_sorted.sam'
    run_command(['samtools', 'sort', '-n', '-O', 'SAM', '-@', sort_threads,
                 '-m', sort_memory, '-o', sorted_sam, output_sam])

    reverse_sam = output_sam[:-4] + '_r.sam'
    f1 = open(sorted_sam, 'r')
    f2 = open(output_bed + '_sorted', 'r')
    file_reverse = open(reverse_sam, 'w')

    index = 0
    index2 = 0
    old_items = 'la'

    while True:
        fr1 = list(itertools.islice(f1, step))
        if len(fr1) == 0:
            break

        list_reverse = []
        list_fr = [line.strip().split('\t') for line in fr1]

        for items in list_fr:
            if not items[0].startswith('@') and int(items[1]) & int(flag):
                continue
            index2 += 1
            reads_A = items[0]

            if reads_A[0] != '@' and reads_A != old_items[0]:
                for row in f2:
                    its = row.strip().split('\t')
                    reads_S = its[0]
                    if reads_A == reads_S:
                        index += 1
                        reversed_seq = multi_sub(items[9], its[1], reverse_fac)
                        items[9] = reversed_seq
                        list_reverse.append(items)
                        old_items = items
                        break
            elif reads_A[0] == '@':
                list_reverse.append(items)
                index += 1
            elif reads_A == old_items[0]:
                index += 1
                list_reverse.append(old_items)

        list_reverse1 = ['\t'.join(map(str, it)) for it in list_reverse]
        if list_reverse1:
            file_reverse.write('\n'.join(list_reverse1) + '\n')

    f1.close()
    f2.close()
    file_reverse.close()
    print(f'  reverse_reads: reversed {index} of {index2} entries')
    if index != index2:
        raise ValueError('Not every mapped read could be restored. Check unique FASTQ read names.')
    os.remove(sorted_sam)
    os.remove(output_sam)
    return reverse_sam


def mapping_bowtie(fastq, reference, threads, mismatch, fqname, outputdir, mul_max=1):
    """
    Map reads to reference using bowtie.
    Returns: (sam_path, unmapped_fastq_path)
    """
    # Bowtie1 can segfault with too many threads on small references; 16 is safe
    bowtie_threads = min(int(threads), 16)
    if int(threads) > 16:
        print(f'  Note: capping bowtie threads at {bowtie_threads} (requested {threads})')

    outputfile = outputdir + fqname + '.sam'
    unmapfastq = outputdir + fqname + '_un_2.fq'

    cmd = ['bowtie', '-k', '1', '-m', mul_max, '-v', mismatch,
           '--best', '--strata', '-p', bowtie_threads, '-x', reference,
           fastq, '-S', outputfile, '--un', unmapfastq]
    with open(outputfile + '.output', 'w') as log:
        run_command(cmd, stderr=log)
    return outputfile, unmapfastq


def sam_to_sorted_bam(samfile, fac, threads, flag='4', sort_memory='256M'):
    """Convert SAM to sorted BAM, index, and clean up SAM."""
    output_bam = samfile[:-4] + fac
    # reverse_reads has already omitted unmapped entries.
    run_command(['samtools', 'sort', '-@', threads, '-m', sort_memory,
                 '-O', 'BAM', '-o', output_bam, samfile])
    run_command(['samtools', 'quickcheck', output_bam])
    run_command(['samtools', 'index', output_bam])
    os.remove(samfile)
    return output_bam


def get_sites(chr, tool_dir, file9, baseanno, outputprefix, options):
    """
    Process a single chromosome/amplicon:
      - split pileup by chromosome
      - format pileup (annotate + compute CR)
      - call m6A sites
    """
    chr2 = chr.split('_AG_converted')[0]
    file_mpi = file9 + '.' + chr2
    file_format = outputprefix + '.referbase.mpi.formatted.txt.' + chr2
    file_CR = outputprefix + '.CR.txt.' + chr2
    file_sites = outputprefix + '.callsites.' + chr2

    # Split pileup by chromosome
    with open(file9) as source, open(file_mpi, 'w') as output:
        for line in source:
            if line.split('\t', 1)[0] == chr:
                output.write(line)

    # Split annotation by chromosome
    if options['baseanno'] != 'None':
        # Keep run-specific annotation inside this output; never reuse stale splits.
        baseanno_chr = outputprefix + '.baseanno.' + chr2
        with open(options['baseanno']) as source, open(baseanno_chr, 'w') as output:
            for line in source:
                if line.split('\t', 1)[0] == chr2:
                    output.write(line)
        run_command([sys.executable, tool_dir + 'm6A_pileup_formatter.py',
                     '--db', baseanno_chr, '-i', file_mpi, '-o', file_format, '--CR', file_CR])
    else:
        run_command([sys.executable, tool_dir + 'm6A_pileup_formatter.py',
                     '-i', file_mpi, '-o', file_format, '--CR', file_CR])

    # Call m6A sites
    run_command([sys.executable, tool_dir + 'm6A_caller.py',
                 '-i', file_format, '-o', file_sites,
                 '-c', options['coverage'], '-C', options['count'],
                 '-r', options['ratio'], '-p', options['pvalue'],
                 '-s', options['signal'], '-R', options['var_ratio'],
                 '--cutoff', options['A_cutoffs'], '--CR', options['background'],
                 '--method', options['method']])

    return chr2


# ============================================================
# Main pipeline
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='GLORI amplicon deep sequencing analysis pipeline')

    # Required
    group_req = parser.add_argument_group('Required')
    group_req.add_argument('-q', '--fastq', required=True,
                           help='Cleaned single-end FASTQ file')
    group_req.add_argument('-f', '--reference', required=True,
                           help='A-to-G converted amplicon FASTA '
                                '(e.g. amplicon.AG_conversion.fa)')
    group_req.add_argument('-f2', '--reference2', required=True,
                           help='Original amplicon FASTA (non-converted)')

    # Optional
    group_opt = parser.add_argument_group('Optional')
    group_opt.add_argument('-i', '--NSdir', default=os.path.dirname(os.path.abspath(__file__)),
                           help='Directory containing pipelines/ (default: this script directory)')
    group_opt.add_argument('-b', '--baseanno', default='None',
                           help='Base-level annotation file (default: None, '
                                'use overall conversion rate)')
    group_opt.add_argument('-pre', '--outname_prefix', default='sample', type=safe_identifier,
                           help='Output file prefix (default: sample)')
    group_opt.add_argument('-o', '--outputdir', default='./output',
                           help='New or empty output directory (default: ./output); no resume')
    group_opt.add_argument('-T', '--Threads', default=4, type=positive_integer,
                           help='Mapping/pileup threads (default: 4); step 1 is single-threaded')
    group_opt.add_argument('--sort-threads', default=2, type=positive_integer,
                           help='samtools sort additional threads, independent of -T (default: 2)')
    group_opt.add_argument('--sort-memory', default='256M',
                           help='Approximate sort memory per thread, e.g. 256M or 1G (default: 256M)')
    group_opt.add_argument('-m', '--mismatch', default=2, type=int, choices=range(4),
                           help='Max mismatches for bowtie (default: 2)')
    group_opt.add_argument('--keep-tmp', action='store_true',
                           help='Preserve intermediate files in tmp/ subdirectory '
                                'for debugging (changed reads, unmapped reads, '
                                'SAM, pileup, per-amplicon files, etc.)')

    # m6A calling filters
    group_filt = parser.add_argument_group('m6A calling filters')
    group_filt.add_argument('-c', '--coverage', default='15',
                            help='A+G coverage cutoff (default: 15)')
    group_filt.add_argument('-C', '--count', default='5',
                            help='A count cutoff (default: 5)')
    group_filt.add_argument('-r', '--ratio', default='0.1',
                            help='Minimum m6A ratio (default: 0.1)')
    group_filt.add_argument('-p', '--pvalue', default='0.005',
                            help='P-value cutoff (default: 0.005)')
    group_filt.add_argument('-adp', '--adjustpvalue', default='0.005',
                            help='FDR-adjusted p-value cutoff (default: 0.005)')
    group_filt.add_argument('-s', '--signal', default='0.8',
                            help='Signal ratio cutoff (default: 0.8)')
    group_filt.add_argument('-R', '--var_ratio', default='0.8',
                            help='AG/Total ratio cutoff (default: 0.8)')
    group_filt.add_argument('--cutoff', dest='A_cutoffs', default='3',
                            choices=[str(n) for n in range(1, 11)] + ['15', '20', 'None'],
                            help='One read A-count cutoff (default: 3); None disables that cutoff')
    group_filt.add_argument('-M', '--max-depth', default=10000, type=positive_integer,
                            help='Max pileup depth per position (default: 10000). '
                                 'Lower this if pileup produces huge files.')
    group_filt.add_argument('--method', default='binomial',
                            choices=['binomial', 'poisson'],
                            help='Statistical method (default: binomial)')

    args = parser.parse_args()
    try:
        for module in ('pysam', 'scipy.stats', 'numpy'):
            importlib.import_module(module)
        from statsmodels.stats.multitest import multipletests
    except ImportError as error:
        parser.error(f'Missing/incompatible Python dependency: {error}. Activate the documented Conda environment.')
    if not re.fullmatch(r'[1-9][0-9]*[KMG]', args.sort_memory):
        parser.error('--sort-memory must be a positive integer followed by K, M or G, e.g. 256M')
    # Resolve paths once so subprocesses can use isolated working directories.
    args.NSdir = os.path.abspath(args.NSdir)
    args.fastq = os.path.abspath(args.fastq)
    args.reference = os.path.abspath(args.reference)
    args.reference2 = os.path.abspath(args.reference2)
    if args.baseanno != 'None':
        args.baseanno = os.path.abspath(args.baseanno)
    args.outputdir = os.path.abspath(args.outputdir)
    required_files = [args.fastq, args.reference, args.reference2,
                      args.reference + '.fai', args.reference2 + '.fai']
    if args.baseanno != 'None':
        required_files.append(args.baseanno)
    for helper in ['pileup_genome_multiprocessing.py', 'get_referbase.py',
                   'm6A_pileup_formatter.py', 'm6A_caller.py', 'm6A_caller_FDRfilter.py']:
        required_files.append(os.path.join(args.NSdir, 'pipelines', helper))
    for path in required_files:
        if not os.path.isfile(path):
            parser.error(f'Required file not found: {path}. Prepare the reference first.')
    index_parts = ['1', '2', '3', '4', 'rev.1', 'rev.2']
    if not any(all(os.path.isfile(args.reference + '.' + part + '.' + extension)
                   for part in index_parts) for extension in ['ebwt', 'ebwtl']):
        parser.error('A complete Bowtie 1 index is missing. Run prepare_amplicon_ref.py first.')
    for executable in ('bowtie', 'samtools'):
        if shutil.which(executable) is None:
            parser.error(f'{executable} is not on PATH. Activate the Conda environment.')
    original = list(SeqIO.parse(args.reference2, 'fasta'))
    converted = list(SeqIO.parse(args.reference, 'fasta'))
    if not original or len(original) != len(converted):
        parser.error('Original and converted reference FASTA records do not match.')
    seen = set()
    for orig, ag in zip(original, converted):
        try:
            safe_identifier(orig.id)
        except argparse.ArgumentTypeError as error:
            parser.error(f'Invalid FASTA identifier {orig.id}: {error}')
        if orig.id in seen:
            parser.error(f'Duplicate FASTA identifier: {orig.id}')
        seen.add(orig.id)
        if (ag.id != orig.id + '_AG_converted'
                or str(ag.seq).upper() != str(orig.seq).upper().replace('A', 'G')):
            parser.error(f'Converted reference does not match the original record: {orig.id}')

    # ========================================================
    # Setup paths
    # ========================================================
    NSdir = args.NSdir.rstrip('/') + '/pipelines/'
    outputdir = args.outputdir.rstrip('/')
    prx = args.outname_prefix
    prepare_output_directory(outputdir)

    genome_ag = args.reference      # A-to-G converted amplicon reference
    genome_orig = args.reference2    # Original amplicon reference
    baseanno = args.baseanno

    # Compile options for get_sites
    opts = {
        'baseanno': baseanno,
        'coverage': args.coverage,
        'count': args.count,
        'ratio': args.ratio,
        'pvalue': args.pvalue,
        'signal': args.signal,
        'var_ratio': args.var_ratio,
        'A_cutoffs': args.A_cutoffs,
        'background': 'overall' if baseanno == 'None' else 'gene',
        'method': args.method,
    }

    outputprefix = f'{outputdir}/{prx}'

    print('=' * 60)
    print(f'GLORI Amplicon Analysis Pipeline')
    print(f'  Input FASTQ     : {args.fastq}')
    print(f'  AG reference    : {genome_ag}')
    print(f'  Original ref    : {genome_orig}')
    print(f'  Annotation      : {baseanno}')
    print(f'  Output prefix   : {outputprefix}')
    print(f'  Threads         : {args.Threads}')
    print(f'  Sort settings   : {args.sort_threads} additional threads; {args.sort_memory} per thread')
    print(f'  Mismatches      : {args.mismatch}')
    print(f'  Background CR   : {opts["background"]}')
    print('=' * 60)

    # ========================================================
    # Step 1: A-to-G conversion of reads
    # ========================================================
    print('\n[Step 1/8] A-to-G conversion of reads...')
    changed_fq = f'{outputdir}/{prx}_AGchanged_2.fq'
    output_bed = f'{outputdir}/{prx}_A.bed'

    change_reads(args.fastq, changed_fq, output_bed, outputdir, 'AG')
    print('  A-to-G conversion done.')

    # ========================================================
    # Step 2: Map changed reads to A-to-G amplicon reference
    # ========================================================
    print('\n[Step 2/8] Mapping changed reads to A-to-G amplicon reference...')
    sam_file, unmap_fq = mapping_bowtie(
        changed_fq, genome_ag,  # bowtie index basename = FASTA path (with .fa)
        args.Threads, args.mismatch, prx, outputdir + '/'
    )
    print(f'  Mapping done. SAM: {sam_file}')

    # ========================================================
    # Step 3: Reverse reads back to original sequence
    # ========================================================
    print('\n[Step 3/8] Reversing reads back to original sequence...')
    reversed_sam = reverse_reads(sam_file, output_bed, 'A', '4',
                                 sort_threads=args.sort_threads, sort_memory=args.sort_memory)
    print(f'  Reverse done. SAM: {reversed_sam}')

    # ========================================================
    # Step 4: Convert SAM to sorted BAM
    # ========================================================
    print('\n[Step 4/8] Converting to sorted BAM...')
    final_bam = sam_to_sorted_bam(reversed_sam, '.sorted.bam', args.sort_threads,
                                '4', sort_memory=args.sort_memory)
    print(f'  BAM: {final_bam}')

    # ========================================================
    # Step 5: Pileup
    # ========================================================
    print('\n[Step 5/8] Running pileup...')
    pileup_file = outputprefix + '.pileup'
    run_command([sys.executable, NSdir + 'pileup_genome_multiprocessing.py',
                 '-P', args.Threads, '-f', genome_ag, '-i', final_bam,
                 '-o', pileup_file, '-m', args.max_depth])
    print('  Pileup done.')

    # ========================================================
    # Step 6: Fix reference base (use original FASTA)
    # ========================================================
    print('\n[Step 6/8] Fixing reference base with original FASTA...')
    mpi_file = outputprefix + '.referbase.mpi'
    run_command([sys.executable, NSdir + 'get_referbase.py', '-input', pileup_file,
                 '-referFa', genome_orig, '-outname_prx', outputprefix])
    print(f'  Reference base corrected. MPI: {mpi_file}')

    # ========================================================
    # Step 7: Per-chromosome formatting + m6A calling
    # ========================================================
    print('\n[Step 7/8] Processing each amplicon...')

    # Get list of chromosomes (amplicons) from the MPI file
    chr_file = outputprefix + '_chrlist'
    with open(mpi_file) as source:
        chr_list = sorted({line.split('\t', 1)[0] for line in source if line.strip()})
    with open(chr_file, 'w') as output:
        output.write('\n'.join(chr_list) + ('\n' if chr_list else ''))

    print(f'  Found {len(chr_list)} amplicons: {chr_list}')

    # Process each amplicon (with multiprocessing if multiple)
    if len(chr_list) == 1:
        get_sites(chr_list[0], NSdir, mpi_file, baseanno, outputprefix, opts)
    elif chr_list:
        with multiprocessing.Pool(min(args.Threads, len(chr_list))) as pool:
            pool.starmap(get_sites, [
                (chrom, NSdir, mpi_file, baseanno, outputprefix, opts) for chrom in chr_list])

    # ========================================================
    # Step 8: Merge results + CR + FDR
    # ========================================================
    print('\n[Step 8/8] Merging results...')

    # Combine formatted files (remove header lines starting with #)
    final_format = outputprefix + '.totalformat.txt'
    formatted_files = [outputprefix + '.referbase.mpi.formatted.txt.' +
                       chrom.split('_AG_converted')[0] for chrom in chr_list]
    merge_files(formatted_files, final_format, skip_comments=True)

    # Combine CR files
    final_CR = outputprefix + '.totalCR.txt'
    pd_CR_list = []
    for chrom in chr_list:
        chr_x = chrom.split('_AG_converted')[0]
        cr_file = outputprefix + '.CR.txt.' + chr_x
        if not os.path.isfile(cr_file):
            raise FileNotFoundError(f'Expected conversion-rate file was not created: {cr_file}')
        # No A+G observations can produce a legitimate empty per-contig CR file.
        # Counts still carry zero/NA observations; do not invent a conversion rate.
        if os.path.getsize(cr_file) == 0:
            continue
        if os.path.exists(cr_file):
            cr1 = pd.read_csv(cr_file, sep='\t', names=[
                'SA', 'Totalcovered_reads', 'Remained A reads',
                'Non-A-to-G ratio', 'Mapped_area'])
            cr2 = cr1[~cr1['SA'].isin([
                '#ALL', '#90%', '#75%', '#50%', '#25%', '#10%',
                '#Median', '#Mean'
            ])][['SA', 'Non-A-to-G ratio']]
            xx_median = cr1[cr1['SA'] == '#Median'][
                'Totalcovered_reads'].values
            if len(xx_median) > 0:
                pd_median = pd.DataFrame({
                    'SA': [f'#Median_{chr_x}'],
                    'Non-A-to-G ratio': xx_median})
                pd_CR_t = pd.concat([pd_median, cr2])
            else:
                pd_CR_t = cr2.copy()
            pd_CR_t['A-to-G_ratio'] = 1 - pd_CR_t['Non-A-to-G ratio']
            pd_CR_list.append(pd_CR_t[['SA', 'A-to-G_ratio']])

    if pd_CR_list:
        pd_CR = pd.concat(pd_CR_list)
        pd_CR.to_csv(final_CR, sep='\t', index=False)
    else:
        pd.DataFrame(columns=['SA', 'A-to-G_ratio']).to_csv(final_CR, sep='\t', index=False)
        print('  No reference positions survived pileup/reference filtering; CR is unavailable.')
    print(f'  Conversion rate file: {final_CR}')

    # Combine m6A site calls
    final_sites1 = outputprefix + '.totalm6A.txt'
    call_files = [outputprefix + '.callsites.' + chrom.split('_AG_converted')[0] +
                  '.' + args.A_cutoffs + '.txt' for chrom in chr_list]
    # The bundled caller writes an empty file when no candidates pass its filters.
    # A missing file is therefore an execution error, not an empty result.
    merge_files(call_files, final_sites1)

    # FDR correction
    final_sites2 = outputprefix + '.totalm6A.FDR'
    run_command([sys.executable, NSdir + 'm6A_caller_FDRfilter.py',
                 '-i', final_sites1, '-o', final_sites2, '-adp', args.adjustpvalue])
    if not os.path.isfile(final_sites2 + '.csv'):
        raise FileNotFoundError('FDR helper returned without creating its result table.')

    # ========================================================
    # Cleanup: move intermediate files to tmp/ or delete
    # ========================================================
    tmp_patterns = [
        f'{outputdir}/{prx}_AGchanged_2.fq',
        f'{outputdir}/{prx}_A.bed_sorted',
        f'{outputdir}/{prx}.sam',
        f'{outputdir}/{prx}.sam.output',
        f'{outputdir}/{prx}_un_2.fq',
        f'{outputdir}/{prx}_r.sam',
        f'{outputprefix}.pileup',
        f'{outputprefix}.referbase.mpi',
        f'{outputprefix}.referbase.mpi.*',
        f'{outputprefix}.referbase.mpi.formatted.txt.*',
        f'{outputprefix}.CR.txt.*',
        f'{outputprefix}.callsites.*',
        f'{outputprefix}.baseanno.*',
        f'{outputprefix}.totalformat.txt',
        f'{outputprefix}.totalm6A.txt',
        f'{chr_file}',
    ]
    if args.keep_tmp:
        tmpdir = os.path.join(outputdir, 'tmp')
        os.makedirs(tmpdir, exist_ok=True)
        for pat in tmp_patterns:
            for path in glob.glob(os.path.join(glob.escape(outputdir), os.path.basename(pat))):
                if os.path.isfile(path):
                    shutil.move(path, tmpdir)
        print(f'\n  Intermediate files kept in: {tmpdir}/')
    else:
        for pat in tmp_patterns:
            for path in glob.glob(os.path.join(glob.escape(outputdir), os.path.basename(pat))):
                if os.path.isfile(path):
                    os.remove(path)

    # ========================================================
    # Summary
    # ========================================================
    print('\n' + '=' * 60)
    print('GLORI Amplicon Analysis Complete!')
    print('=' * 60)
    print(f'\nOutput files:')
    print(f'  Aligned BAM          : {final_bam}')
    print(f'  Conversion rate (CR) : {final_CR}')
    print(f'  m6A sites (FDR)      : {final_sites2}.csv')
    print(f'\nOutput columns in m6A file:')
    print(f'  Chr, Sites, Strand, Gene, CR, AGcov, Acov, Genecov,')
    print(f'  Ratio (m6A level), Pvalue, P_adjust (FDR)')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        print('Analysis stopped. Keep this output for diagnosis and use a new -o directory for a rerun.',
              file=sys.stderr)
        sys.exit(1)
