# Installation and batch analysis

[Home](../README.md) · [简体中文](usage.zh-CN.md) · [Input requirements](inputs.en.md)

This guide starts from a new cluster account with Conda available. All commands run in a Linux Bash terminal. Copy only the commands inside code blocks; do not copy terminal prompts. Replace every `/path/to/your/...` and `YOUR_CPU_PARTITION` placeholder with your own value before running it.

## 1. Prepare the code and inputs

Download the repository with GitHub **Code → Download ZIP**, extract it, and place the complete code folder in your cluster account. Keep the Python files, `core/`, `LICENSES/`, licenses and `environment.yml` together. Alternatively, clone the repository if Git is available and your account has access.

Upload only your demultiplexed FASTQ files and their original amplicon FASTA references in addition to the code. One sample uses one read end, normally R1. R2, sample index reads, genome indexes, annotation databases, previous analysis results and local Conda environments are not required. There are no example input datasets in this repository. See the [input requirements](inputs.en.md) before preparing the sample sheet.

Enter the code folder and confirm its contents:

```bash
cd /path/to/your/glori-amplicon-batch
pwd
ls
```

`cd` changes directory. `pwd` prints the current directory. `ls` lists files; you should see `run_batch.py`, `environment.yml`, `samples.tsv` and `core`.

The code, inputs and output directory must be on storage accessible to both login and compute nodes. Do not use a login-node-only temporary directory.

FASTQ and output paths must contain no whitespace, including spaces in parent-directory names, because Trim Galore 0.6.10 does not support them. Input and submission checks reject these paths early.

## 2. Create your Conda environment

If `conda` is not found, first load your cluster's Conda module or follow its administrator's Conda setup instructions. Then create the environment in your own account:

```bash
mkdir -p "$HOME/.conda/envs" "$HOME/.conda/pkgs"
conda config --prepend envs_dirs "$HOME/.conda/envs"
conda config --prepend pkgs_dirs "$HOME/.conda/pkgs"
export CONDA_CHANNEL_PRIORITY=strict
conda env create --file environment.yml
conda activate glori_amplicon
python --version
```

`mkdir -p` creates your personal environment and package-cache directories if needed. The two `conda config` commands set your account's preferred directories and persist across sessions; they do not modify the shared system installation. `export` enables strict channel priority for this shell. `conda env create` reads `environment.yml`, downloads dependencies, and creates the named environment. It can take several minutes and needs access to Conda package servers. `conda activate` selects that environment; `python --version` should report Python 3.10.

The environment includes Python, NumPy, pandas, SciPy, statsmodels, Biopython and pysam, plus Bowtie 1, samtools, Trim Galore, Cutadapt and FastQC. Slurm commands are supplied by the cluster, not by this environment.

Create the environment once. In future sessions, only activate it. If an environment named `glori_amplicon` already exists, check whether it meets this specification; do not delete a working environment merely to repeat this guide. You can create a separate named environment with `conda env create --name YOUR_NEW_ENV --file environment.yml` and activate that name instead. The runner checks that Python belongs to the active environment, not that it has a particular name.

## 3. Fill in samples.tsv

`samples.tsv` initially contains a header only. Open it in a text editor and add one row per sample, following [Input requirements](inputs.en.md). Separate the three fields with **Tab**, not spaces or commas. You may provide different references for different samples or reuse one reference for multiple samples.

```bash
cat samples.tsv
```

`cat` displays the sheet so that you can check its sample names and paths. Relative input paths are resolved from the directory containing `samples.tsv`, not from the output directory. Names must be unique even when case is ignored. Do not enter literal placeholder paths as if they were real files.

## 4. Check inputs and dependencies

Normally, stay on the login node, activate the environment and run:

```bash
conda activate glori_amplicon
python run_batch.py --samples samples.tsv --check-only
```

`--samples` selects your sheet. `--check-only` checks paths, all reference sequences, the first FASTQ record, Python dependencies, tools and Slurm commands without submitting jobs. Full FASTQ integrity scans and FastQC run later inside the sample jobs. The script limits numerical-library threads before scientific imports to avoid unnecessary OpenBLAS allocation.

You do **not normally need to enter a CPU node first**. If the login node still blocks imports or memory allocation, or cluster policy requires it, request a small interactive CPU allocation and repeat the check:

```bash
srun --partition=YOUR_CPU_PARTITION --cpus-per-task=1 --mem=2G --time=00:20:00 --pty bash
conda activate glori_amplicon
cd /path/to/your/glori-amplicon-batch
python run_batch.py --samples samples.tsv --check-only
exit
```

`srun` requests temporary compute resources and may wait in the queue. Replace the partition name with an actual CPU partition. `exit` releases this temporary allocation and returns to the login shell. If necessary, run the submission command in step 5 before `exit`, because submission repeats dependency checks. The batch jobs request their own resources; the interactive allocation does not supply the 48 CPUs used by four concurrent samples. If the administrator prohibits job submission from compute nodes, follow the site's supported submission route.

For input-only checks on a machine without analysis dependencies, including Windows, use `python run_batch.py --samples samples.tsv --validate-inputs`. This does not validate the cluster environment.

## 5. Submit all samples

Choose a new output directory and explicitly specify your CPU partition:

```bash
python run_batch.py --samples samples.tsv --out /path/to/your/results/run01 --partition YOUR_CPU_PARTITION
```

`--out` is the results directory that the script will create. It must not already exist for a new batch. `--partition` is mandatory for a new batch; both the sample array and report job use it. The runner does not silently accept a cluster's default partition. The exact compute node is chosen by Slurm, not by the directory or host from which you submit.

The runner prints the sample-array and report job IDs. It snapshots the software, records inputs and resource settings, submits the sample array with a hold, submits a report job depending on completion of that array, and releases the array only after both job IDs are recorded. You can close the terminal after successful submission.

Default resource requests:

| Resource | Default |
|---|---|
| Concurrent samples | At most 4 |
| CPUs per sample | 12 |
| Memory per sample | 16 GB |
| Time per sample | 8 hours |
| Report job | 1 CPU, 2 GB, 1 hour |

Four active samples request up to 48 CPUs and 64 GB in total. They are not forced onto the same node; available memory, partition policy and other users' jobs affect scheduling. Sorting uses 2 additional threads with a 256 MB per-thread limit independently of Bowtie threads. The A-to-G read conversion step is still single-threaded, so raising `--cpus` does not accelerate that step.

If necessary, choose different resources when submitting a new output directory:

```bash
python run_batch.py --samples samples.tsv --out /path/to/your/results/run02 --partition YOUR_CPU_PARTITION --max-parallel 4 --cpus 8 --memory 16G --time 12:00:00
```

This uses at most 32 CPUs across four samples and allows 12 hours per sample. `--cpus` must be at least 3. Add `--account YOUR_SLURM_ACCOUNT` only if your cluster requires a billing account. Add `--keep-all` at submission to retain the large intermediates that the wrapper would otherwise clean up after successful export; it does not undo cleanup internal to the analysis core.

## 6. Follow progress

```bash
squeue -u "$USER"
cat /path/to/your/results/run01/qc_summary.tsv
ls /path/to/your/results/run01/logs
```

`squeue` shows your active jobs, partition, state and assigned node. A report job waiting with `Dependency` is normal. Sample logs are in `logs/array-<job>_<index>.log`; report logs are in `logs/report-<job>.log`. The first summary marks samples as queued. Summary files are regenerated by the report job after the array ends, rather than updated continuously. During execution, inspect each attempt's `status.json` and `logs/pipeline.log` for live progress.

Each sample automatically runs: complete FASTQ validation/counts → raw FastQC → Illumina adapter and quality trimming → clean FASTQ validation/counts and FastQC → reference preparation/indexing → GLORI conversion, mapping, A restoration, sorting, pileup and calling → BAM/count validation → per-site export. FastQC WARN/FAIL flags are recorded without automatically treating an otherwise successful analysis as a software failure.

## 7. Read the results

Download the **complete results directory** to your computer and open `report.html` in a browser. Keeping the directory structure preserves links to FastQC and logs. The report currently uses Chinese explanatory text; tables have stable English column names.

| File | Contents |
|---|---|
| `report.html` | Offline QC overview, sample state and report links |
| `qc_summary.tsv` | Read counts, retention, Q30, mapping, reverse reads, significant-site counts, warnings and failures |
| `all_A_sites.tsv` | Every A position in each sample's original reference, including missing observations |
| `significant_sites.tsv` | Significant calls from successful samples, with original per-sample FDR values |
| `samples/<sample>/attempt_001/` | Per-sample reads, references, BAM/index, counts, QC, logs and status |
| `software/`, `rounds/`, `batch.json` | Software snapshot and execution/submission records |

Key columns in `all_A_sites.tsv`:

| Column | Meaning |
|---|---|
| `Sample`, `Chr`, `Sites`, `Strand` | Sample, FASTA record ID, 1-based amplicon coordinate, and `+` strand |
| `Acov`, `Gcov`, `AGcov` | A, G and A+G counts after the A-cutoff filter |
| `Ratio`, `Percent` | A/(A+G), and that fraction multiplied by 100 |
| `*_all` | Counts/ratios before the A-cutoff filter, still subject to alignment and pileup filtering/depth cap |
| `Signal` | Fraction of total pileup observations retained by the A-cutoff filter |
| `Background_CR` | Amplicon-level A-to-G conversion rate reported by the core |
| `FDR_status` | `significant` means a reported FDR call; `not_reported` does not distinguish individual filtering reasons |
| `Observation_status` | Observed, missing usable record, no A/G after cutoff, or sample failure information |

Use `Percent` for the estimated m6A percentage after the default A-cutoff filter. It is not automatically corrected for conversion background, PCR bias or other experimental effects. `NA` means no usable estimate; it must not be interpreted as 0% methylation. A ratio can exist without a significant call. Failed samples produce missing estimates rather than partial results presented as valid. Check `Status` and QC before interpreting results.

The significant table is allowed to contain only a header when a sample succeeds but has no significant sites. The core applies FDR to its prefiltered candidates within each sample; the batch report does not perform a second correction across samples. Coordinates refer to your FASTA, not automatically to a genome.

Analysis defaults: Illumina adapter trimming, Q20, minimum length 25, adapter overlap 5, error rate 0.1; Bowtie 1 with at most 2 mismatches and unique mapping; pileup maximum depth 10,000; A-cutoff 3; gene/amplicon background; minimum A+G coverage 15, A count 5, ratio 0.1, signal 0.8, A+G/total 0.8; binomial P threshold 0.005 and FDR threshold 0.005. The pileup cap means deeply sequenced samples' site counts need not equal all mapped reads.

No molecular-UMI deduplication, paired-end mapping or primer trimming is performed. Ensure your references describe the sequenced insert in the expected orientation. Unsupported reverse alignments cause an explicit failure. Biological suitability of the defaults remains a separate experimental decision.

## 8. Retry failures or regenerate the report

Wait until both the sample array and its report job have left `squeue`. Keep the original input files unchanged, inspect the error, correct the resource or environment issue, and retry only failed/interrupted samples:

```bash
python run_batch.py --samples samples.tsv --out /path/to/your/results/run01 --partition YOUR_CPU_PARTITION --retry-failed
```

Successful samples are skipped. Failed samples receive a new attempt directory; previous attempts and a copy of the earlier report remain available. To increase memory or time for the retry, add `--memory 32G --time 16:00:00`. If you change input data, references or analysis software, use a new output directory for a fresh batch. Retries reuse the original software snapshot, so updating files beside the launcher does not update an existing batch's analysis code.

If the automatic report job failed but all sample jobs have ended, regenerate the report:

```bash
python run_batch.py --out /path/to/your/results/run01 --summarize
```

This does not rerun analysis. It refuses to run while recorded jobs remain active. If a submission times out and its job IDs are uncertain, inspect the recorded `rounds/.../submit.log` and your Slurm queue before doing anything that might create duplicate jobs. Do not delete execution records to bypass this protection.

## 9. Software checks

On a Linux machine with the complete environment, developers can run:

```bash
RUN_GLORI_INTEGRATION=1 python -m unittest discover -s tests -v
```

Unit tests cover input validation, missing values, QC exports, failed-sample isolation, submission/retry behavior and resource handling. A real-tool integration test creates tiny temporary reads, runs trimming/QC/analysis and compares the wrapper with a direct core run. Slurm itself is mocked in tests; testing does not validate any particular cluster or establish biological accuracy. See [GitHub Actions](https://github.com/luckingclark/glori-amplicon-batch/actions) for actual Linux results.
