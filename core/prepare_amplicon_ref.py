#!/usr/bin/env python
"""
Prepare amplicon reference for GLORI amplicon deep sequencing analysis.

Input: FASTA file of amplicon sequences (one record per amplicon)
Output:
  - {prefix}.AG_conversion.fa    : A-to-G converted reference
  - bowtie index files           : for the converted reference
  - {prefix}.baseanno            : base-level annotation (for per-amplicon CR)
  - FASTA .fai indices           : original and converted FASTA

Usage:
  python prepare_amplicon_ref.py -f amplicons.fa -pre my_panel -o ./ref

Requirements: bowtie >= 1.3.1, samtools >= 1.10, biopython
"""

import argparse
import os
import subprocess
import sys
import re
import glob
import shutil
from Bio import SeqIO
from Bio.Seq import Seq


def main():
    parser = argparse.ArgumentParser(
        description='Prepare amplicon reference for GLORI amplicon analysis'
    )
    parser.add_argument('-f', '--fasta', required=True,
                        help='Amplicon FASTA file (one record per amplicon)')
    parser.add_argument('-p', '--threads', default=4, type=int,
                        help='Threads for bowtie-build (default: 4)')
    parser.add_argument('-pre', '--prefix', default='amplicon',
                        help='Output prefix (default: amplicon)')
    parser.add_argument('-o', '--outputdir', default='./',
                        help='Output directory (default: ./)')
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('Threads must be at least 1.')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.prefix):
        parser.error('Prefix must start with a letter or number and use only letters, numbers, _, . or -.')
    args.fasta = os.path.abspath(args.fasta)
    args.outputdir = os.path.abspath(args.outputdir)
    for executable in ('bowtie-build', 'samtools'):
        if shutil.which(executable) is None:
            parser.error(f'{executable} is not on PATH. Activate the Conda environment.')
    records = list(SeqIO.parse(args.fasta, 'fasta'))
    if not records:
        parser.error('Input FASTA contains no records.')
    seen = set()
    for record in records:
        if (not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', record.id)
                or '_AG_converted' in record.id):
            parser.error(f'Invalid/reserved FASTA identifier: {record.id}. Use letters, numbers, _, . or -.')
        if record.id in seen:
            parser.error(f'Duplicate FASTA identifier: {record.id}')
        seen.add(record.id)
        if not record.seq or re.search('[^ACGTN]', str(record.seq).upper()):
            parser.error(f'{record.id}: sequence must be nonempty DNA containing A, C, G, T or N.')
        if 'A' not in str(record.seq).upper():
            parser.error(f'{record.id}: original RNA-oriented reference must contain at least one A.')

    os.makedirs(args.outputdir, exist_ok=True)
    prefix = os.path.join(args.outputdir, args.prefix)
    if glob.glob(glob.escape(prefix) + '.*'):
        parser.error(f'Reference files with prefix {prefix} already exist. Use a new -o or -pre.')

    # ============================================================
    # 1. Create A-to-G converted reference FASTA
    # ============================================================
    ag_fa = prefix + '.AG_conversion.fa'
    print(f'[1/4] Creating A-to-G converted reference: {ag_fa}')
    with open(ag_fa, 'w') as out:
        for original_record in records:
            record = original_record[:]
            original_seq = str(record.seq).upper()
            converted_seq = original_seq.replace('A', 'G')
            record.seq = Seq(converted_seq)
            record.id = record.id + '_AG_converted'
            record.description = ''
            SeqIO.write(record, out, 'fasta')

    # ============================================================
    # 2. Build bowtie index for the A-to-G converted reference
    # ============================================================
    print(f'[2/4] Building bowtie index...')
    subprocess.run(['bowtie-build', '--threads', str(args.threads), '-q', ag_fa, ag_fa], check=True)
    subprocess.run(['samtools', 'faidx', args.fasta], check=True)
    subprocess.run(['samtools', 'faidx', ag_fa], check=True)

    # ============================================================
    # 3. Generate base-level annotation from amplicon FASTA
    #    Format: chr, pos_0, pos_1, dir, gene, trans, name, isoform, biotype
    #    Only annotates positions with A (+ strand) or T (- strand)
    #    because these are the only potential m6A sites.
    # ============================================================
    anno_file = prefix + '.baseanno'
    print(f'[3/4] Generating base annotation: {anno_file}')
    total_sites = 0
    with open(anno_file, 'w') as anno:
        for record in records:
            seq = str(record.seq).upper()
            amp_name = record.id
            for i, base in enumerate(seq):
                pos_0 = i
                pos_1 = i + 1

                if base == 'A':
                    # + strand: A is a potential m6A site
                    anno.write('\t'.join([
                        amp_name, str(pos_0), str(pos_1), '+',
                        amp_name, amp_name,
                        f'{amp_name}_{pos_1}', 'NA', 'amplicon'
                    ]) + '\n')
                    total_sites += 1

                elif base == 'T':
                    # - strand: T in reference = A on opposite strand
                    anno.write('\t'.join([
                        amp_name, str(pos_0), str(pos_1), '-',
                        amp_name, amp_name,
                        f'{amp_name}_{pos_1}', 'NA', 'amplicon'
                    ]) + '\n')
                    total_sites += 1

    # ============================================================
    # 4. Summary
    # ============================================================
    print(f'[4/4] Done!')
    print(f'  A-to-G reference : {ag_fa}')
    print(f'  Bowtie index      : {ag_fa}.*.ebwt')
    print(f'  Base annotation   : {anno_file}')
    print(f'  FASTA indices     : {args.fasta}.fai; {ag_fa}.fai')
    print(f'  Potential m6A sites annotated: {total_sites}')

    # Print amplicon stats
    print(f'\n  Amplicon summary:')
    for record in records:
        seq = str(record.seq).upper()
        a_count = seq.count('A')
        t_count = seq.count('T')
        print(f'    {record.id}: {len(seq)} bp, '
              f'A sites: {a_count} (+), {t_count} (-)')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
