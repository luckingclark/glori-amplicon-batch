# GLORI Amplicon Batch

A command-line workflow for researchers analyzing targeted GLORI amplicon sequencing on a Slurm cluster. Provide a sample sheet, one FASTQ read end per sample, and the corresponding original amplicon references. The workflow runs trimming, QC and GLORI analysis, then combines per-sample QC, A-site methylation estimates and significant calls into tables and an offline report.

This is an implemented, independent adaptation of GLORI-tools, not an official GLORI release. Its scope is single-end amplicon analysis; it does not perform paired-end mapping, molecular-UMI extraction, deduplication or primer trimming.

| Guide | English | Chinese version |
|---|---|---|
| Installation, commands and troubleshooting | [User guide](docs/usage.en.md) | [User guide](docs/usage.zh-CN.md) |
| FASTQ, reference and sample-sheet requirements | [Input requirements](docs/inputs.en.md) | [Input requirements](docs/inputs.zh-CN.md) |

## Before you start

You need a **Linux cluster with Bash, Conda and Slurm**, permission to submit CPU jobs, and storage shared by login and compute nodes. Slurm is supplied by the cluster, not installed by this project. No GPU is required. Windows supports input-only validation, not this analysis workflow.

Prepare your own inputs:

| Input | Required content |
|---|---|
| FASTQ | Demultiplexed, standard four-line Phred+33 reads; `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`. Choose one read end, normally R1, with unique read identifiers. |
| Amplicon FASTA | Uncompressed, original **unconverted** insert sequence in RNA-sense orientation, using DNA letters A/C/G/T/N. Exclude sequencing adapters and sample indexes. Each record must contain A and have a unique ID. |
| Sample sheet | UTF-8 text with exactly three Tab-separated columns: `sample`, `fastq`, `reference`; one row per sample. The supplied `samples.tsv` has only a header. |

Use libraries without molecular UMIs: i5/i7 sample indexes are not molecular UMIs. You do not need R2, a whole-genome reference, prebuilt indexes or an annotation database. Each sample may use its own reference, or samples may share a reference. See [input requirements](docs/inputs.en.md) for identifier rules and reserved names.

FASTQ and output paths, including parent directories, must contain **no whitespace** because Trim Galore 0.6.10 does not support it. Use a new output directory; existing directories are refused for new batches. The output path must also contain no `%` character. No experimental data or downloadable demonstration dataset is included.

## First run: one sample using your own data

### 1. Get the code and create the environment

Download this repository with **Code → Download ZIP**, extract it, and transfer the complete code directory plus your FASTQ and FASTA to shared cluster storage. Keep `core/`, `LICENSES/` and the root files together.

Run these commands in a Bash terminal. At each prompt, enter the actual path without surrounding quotes; the commands quote it where needed. Continue in the same terminal for the following steps.

```bash
read -r -p "Full path to the extracted code directory: " GLORI_CODE
cd "$GLORI_CODE"
ls run_batch.py environment.yml
mkdir -p "$HOME/.conda/envs" "$HOME/.conda/pkgs"
conda config --prepend envs_dirs "$HOME/.conda/envs"
conda config --prepend pkgs_dirs "$HOME/.conda/pkgs"
export CONDA_CHANNEL_PRIORITY=strict
conda env create --file environment.yml
conda activate glori_amplicon
python --version
```

The configuration commands give your personal environment and package-cache directories priority without changing the shared Conda installation. Environment creation needs access to the package channels and is done once; later sessions only need activation. If Conda is unavailable or this environment name already exists, follow the [installation guide](docs/usage.en.md#2-create-your-conda-environment). The Python version should be 3.10.

The checked-in [environment.yml](environment.yml) specifies the dependencies:

| Component | Environment specification |
|---|---|
| Python and numerical packages | Python 3.10; NumPy 1.26.4; pandas 2.2.3; SciPy 1.14.1; statsmodels 0.14.4 |
| Sequence/BAM packages | Biopython >=1.80,<1.86; pysam >=0.22,<0.24 |
| Mapping and alignment tools | Bowtie **1** 1.3.1; samtools >=1.10,<2 |
| Trimming and QC | Trim Galore 0.6.10; Cutadapt >=4,<5; FastQC >=0.12,<0.13 |

Conda resolves their transitive dependencies. The configuration uses conda-forge and bioconda, with defaults excluded. It is an environment specification, not a lockfile pinning every package build.

### 2. Create a minimal sample sheet

This example runs one of **your** samples and names it `sample01` in the results. Enter absolute paths to existing files. It creates or replaces `samples.first-run.tsv` in the code directory; use that filename only for this example.

```bash
read -r -p "Full path to your FASTQ: " GLORI_FASTQ
read -r -p "Full path to your original amplicon FASTA: " GLORI_REFERENCE
printf 'sample\tfastq\treference\nsample01\t%s\t%s\n' "$GLORI_FASTQ" "$GLORI_REFERENCE" > samples.first-run.tsv
cat samples.first-run.tsv
```

`printf` writes actual Tab separators. For multiple samples, prepare additional rows with unique sample IDs and the appropriate FASTQ/reference paths; do not combine R1 and R2 into one input. Relative paths, if used, are resolved from the sample sheet's directory.

### 3. Check, then submit

```bash
python run_batch.py --samples samples.first-run.tsv --check-only
```

This checks file paths, the reference, the first FASTQ record, dependencies and Slurm commands, without submitting a job. Full FASTQ validation and FastQC run inside the compute jobs. On success the command exits with status 0 and prints a confirmation; command-line messages are currently mostly Chinese. If it reports an error, resolve that before submitting.

Normally, this check and submission run on the login node. The script limits numerical-library threads before imports. If login-node restrictions still prevent dependency checks, use the small CPU allocation described in the [guide](docs/usage.en.md#4-check-inputs-and-dependencies), subject to your cluster's submission policy.

```bash
read -r -p "CPU partition name on your cluster: " GLORI_PARTITION
read -r -p "Full path to a NEW results directory: " GLORI_OUT
python run_batch.py --samples samples.first-run.tsv --out "$GLORI_OUT" --partition "$GLORI_PARTITION"
squeue -u "$USER"
```

`run_batch.py` is the user entry point. It calls the reference preparation and analysis scripts in `core/`; you do not run those separately for a batch. Both analysis and report jobs use the partition you specify. Slurm chooses the actual nodes. A successful submission prints **two job IDs**: the sample array and the report job. You may close the terminal after submission succeeds.

Default requests are **at most 4 concurrent samples**, each with **12 CPUs, 16 GB and 8 hours**. With four samples, that is up to 48 CPUs and 64 GB total; the samples need not run on one node. The dependent report job requests 1 CPU, 2 GB and 1 hour. These are resource requests, not measured minimum requirements or runtime guarantees. Adjust them with `--max-parallel`, `--cpus`, `--memory` and `--time` as described in the [guide](docs/usage.en.md#5-submit-all-samples).

### 4. Confirm completion

Wait for the sample array and report job to finish. A report job waiting with `Dependency` is normal. In the same terminal, inspect the final table:

```bash
cat "$GLORI_OUT/qc_summary.tsv"
ls "$GLORI_OUT/report.html" "$GLORI_OUT/all_A_sites.tsv" "$GLORI_OUT/significant_sites.tsv"
```

If you reopened the terminal, set `GLORI_OUT` to the results path again first. A completed, successful example has a `sample01` row with **`Status` equal to `SUCCESS`**, the three files listed above, and usable per-sample QC links in `report.html`. For a batch, check every sample's status. Download the complete results directory to preserve report links, then open `report.html` locally.

**File existence, disappearance from `squeue`, or a successful report job alone does not establish analysis success.** Initial reports are written while jobs are queued, and the report job can complete while listing failed samples. FastQC WARN/FAIL flags are recorded separately from software execution status; inspect them before interpreting the estimates. A successful sample may have no significant sites, in which case `significant_sites.tsv` can contain only its header.

## Read the results

| File | What to use it for |
|---|---|
| `report.html` | Offline overview of sample status, QC, warnings and links; explanatory text is currently Chinese. |
| `qc_summary.tsv` | Read counts, trimming retention, Q30, mapping, reverse reads, significant-site counts and errors. |
| `all_A_sites.tsv` | Every positive-strand A in each sample's original reference, with filtered/unfiltered counts, ratios and missing-observation reasons. |
| `significant_sites.tsv` | Combined significant calls with sample names and the original per-sample FDR values. |
| `samples/` | Separate sample attempts containing cleaned reads, reference/index files, BAM/index, counts, QC and logs. |
| `software/`, `rounds/`, `batch.json` | Software snapshot, submission settings and execution records. |

In `all_A_sites.tsv`, **`Ratio = Acov / (Acov + Gcov)`** after the default A-cutoff filter, and **`Percent = 100 * Ratio`** is the estimated m6A percentage. Counts and ratios ending in `_all` precede that filter but remain subject to alignment, pileup quality and depth limits. `Sites` is a **1-based amplicon coordinate**, not an inferred genome coordinate.

Treat `NA` as missing, not 0% methylation. A numeric ratio does not by itself indicate a significant call, and `not_reported` does not establish absence of methylation. FDR correction is performed on the core's prefiltered candidates within each sample, not across the combined batch. Ratios are not automatically corrected for conversion background or PCR bias.

The batch uses Illumina adapter trimming, Q20 and minimum length 25; unique Bowtie mapping with at most 2 mismatches; pileup depth capped at 10,000; and A-cutoff 3. Complete calling thresholds and column definitions are in the [results guide](docs/usage.en.md#7-read-the-results). Unsupported reverse alignments cause a failure. Raising CPU requests does not parallelize the first A-to-G conversion step.

## More samples, retries and checks

Use the same entry point with your completed multi-sample sheet. Failed samples do not stop other samples. Once all recorded jobs have ended, `--retry-failed` creates new attempts for unsuccessful samples and skips validated successes; keep the original inputs unchanged. `--summarize` only regenerates reports. See [retry instructions](docs/usage.en.md#8-retry-failures-or-regenerate-the-report). Retries retain the original software snapshot; use a new batch after changing analysis code or inputs.

For developers, the regression and real-tool integration checks can be run on Linux in the complete environment:

```bash
RUN_GLORI_INTEGRATION=1 python -m unittest discover -s tests -v
```

[Linux validation of commit `b2140b0`](https://github.com/luckingclark/glori-amplicon-batch/actions/runs/35960259874) passed **22 tests with no skips**, including real trimming, FastQC, Bowtie, samtools and GLORI analysis, and a comparison of batch-worker output with a direct core run. The tests generate small inputs in temporary directories. Slurm submission is tested with a mocked scheduler; this does not establish compatibility with every cluster or biological accuracy. The [workflow configuration](.github/workflows/checks.yml) records how those checks run.

For a reproducible software problem, use the repository's [issue tracker](https://github.com/luckingclark/glori-amplicon-batch/issues) with versions and a minimal error description. Remove private paths, sample identifiers and sequencing data from reports.

## Citation and provenance

Use [CITATION.cff](CITATION.cff) to cite this adaptation, and record the repository URL and the commit actually used. Cite the original methods and software below as applicable; using this program does not establish which wet-lab method generated the data.

1. Liu, C., Sun, H., Yi, Y., et al. (2023). **Absolute quantification of single-base m6A methylation in the mammalian transcriptome using GLORI.** *Nature Biotechnology*, 41, 355–366. https://doi.org/10.1038/s41587-022-01487-9 . First published online in 2022; the volume year is 2023.
2. Sun, H., Lu, B., Zhang, Z., et al. (2025). **Mild and ultrafast GLORI enables absolute quantification of m6A methylome from low-input samples.** *Nature Methods*, 22, 1226–1236. https://doi.org/10.1038/s41592-025-02680-9 . Cite this when the GLORI 3.0 experimental method is used; software use alone does not establish which wet-lab method generated a dataset.
3. Lu, B. (2024). **Codes for “Mild and ultrafast GLORI enables absolute quantification of m6A methylome from low-input samples”.** Zenodo. https://doi.org/10.5281/zenodo.14233421 . This is a separate collection of companion analysis scripts, not a release of this amplicon adaptation.
4. **GLORI-tools**, by Cong Liu and contributors. https://github.com/liucongcas/GLORI-tools . Its upstream software citation points to https://doi.org/10.5281/zenodo.7014168 . This repository is not asserted to be an exact copy of that archived version.
5. **RNA-m5C**, SYSU-zhanglab and Jianheng Liu. https://github.com/SYSU-zhanglab/RNA-m5C . This is an upstream code source acknowledged by GLORI-tools.

Code sources and local changes are described in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Cite analysis dependencies as appropriate for your work.

## License

Original additions, modifications and documentation are covered by the [root MIT license](LICENSE), whose scope does not replace upstream copyright. Derived code retains the complete [GLORI-tools MIT notice](LICENSES/GLORI-tools-MIT.txt) and [RNA-m5C MIT notice](LICENSES/RNA-m5C-MIT.txt). Keep applicable upstream notices when redistributing the code; separately installed dependencies retain their own licenses.

This independent adaptation does not imply upstream collaboration or endorsement. Paper materials and the separately licensed Zenodo companion scripts are not bundled.

## Author

**PKU-Gaolab, Ming-Ao Lu** — authors of this adaptation and its documentation.
