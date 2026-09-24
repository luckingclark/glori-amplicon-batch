# Input requirements

[Home](../README.md) · [User guide](usage.en.md) · [Chinese version](inputs.zh-CN.md)

Prepare your own data; this repository contains no experimental or downloadable synthetic input datasets.

## FASTQ: one file per sample

- Already demultiplexed, standard four-line FASTQ with Phred+33 qualities; `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`.
- The full FASTQ path and chosen results path must not contain spaces, tabs or other whitespace: Trim Galore 0.6.10 does not support them. The runner rejects these paths before dependency checks and job submission.
- Use one read end, normally R1, consistently with the reference orientation. Do not combine R1 and R2 into the input. Supply the original reads; the workflow performs Illumina adapter and quality trimming.
- Sequences may contain A/C/G/T/N and must match quality-string lengths. Read identifiers must be unique within a sample for reliable restoration of the original A bases.
- This workflow assumes no molecular UMIs. i5/i7 sample indexes are not molecular UMIs; there is no UMI extraction or deduplication step.

## Reference: original amplicon FASTA

- Supply the original, **unconverted** target insert sequence, written in the RNA-sense orientation using DNA letters (T, not U). Do not include sequencing adapters or sample indexes.
- One or more records per sample, each with a unique `>record_id` and a nonempty sequence containing A. Allowed sequence letters are A/C/G/T/N; lowercase is normalized. N is accepted but may reduce mapping.
- Record IDs begin with a letter/digit and contain only letters, digits, `_`, `.` or `-`. Do not use `_AG_converted` within IDs; do not use reserved IDs `ALL`, `ELSE`, `Median`, `Mean`, `discard`, `control` or names beginning `Median_`.
- The first whitespace-delimited token after `>` is the ID. Subsequent header descriptions are not used. Output coordinates are 1-based within this reference.
- Use a plain-text, uncompressed FASTA. No genome, GTF/GFF, prebuilt index or base annotation is needed; the pipeline generates its own converted reference, Bowtie index and annotation.
- Multiple references that become identical after A→G conversion cannot be distinguished reliably. The workflow warns and excludes multimapped reads under its unique-mapping settings.

## samples.tsv: three Tab-separated columns

The supplied file is a header-only template. Add your own rows without changing column order:

| Column | Enter |
|---|---|
| `sample` | Unique sample ID, beginning with a letter/digit and containing only letters, digits, `_`, `.` or `-`; no `_AG_converted` |
| `fastq` | Existing FASTQ path, e.g. `/path/to/your/reads/sample_R1.fq.gz` |
| `reference` | Existing original FASTA path, e.g. `/path/to/your/references/amplicons.fa` |

One row represents one sample. Sample names are unique ignoring case. Use literal Tab separators and UTF-8 text; a UTF-8 BOM is tolerated. Do not add extra columns or quote paths as shell strings. Absolute paths are simplest; relative paths are resolved from the sample sheet's directory. A reference can be reused across samples. All files must be accessible on the cluster compute nodes and must remain unchanged during a run and its retries.

Only the code, this completed sheet, the selected FASTQ files and the corresponding FASTA files need to be available on the cluster. Keep personal input sheets and results out of software commits.
