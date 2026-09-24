# GLORI Amplicon Batch

Automated multi-sample, single-end GLORI amplicon analysis on Slurm, with per-sample QC and combined A-site quantification.

**Authors of this adaptation: PKU-Gaolab, Ming-Ao Lu.**

This independent, unofficial adaptation packages the GLORI amplicon analysis core with a batch runner: raw-read validation and FastQC, Trim Galore, reference/index preparation, GLORI calling, output checks, and an offline report. It is not an official GLORI release and does not imply collaboration with or endorsement by the upstream authors.

## Documentation

| Document | English | 简体中文 |
|---|---|---|
| Installation and complete workflow | [User guide](docs/usage.en.md) | [操作指南](docs/usage.zh-CN.md) |
| Input files and sample sheet | [Input requirements](docs/inputs.en.md) | [输入文件规范](docs/inputs.zh-CN.md) |

## Installation and submission

Requires Linux, Bash, Conda, and access to Slurm. Supply your own demultiplexed R1 (or another single read end) and original RNA-oriented amplicon FASTA. This workflow assumes no molecular UMIs; i5/i7 sample indexes are not UMIs.

From the downloaded code directory:

```bash
export CONDA_CHANNEL_PRIORITY=strict
conda env create --file environment.yml
conda activate glori_amplicon
```

Fill the header-only `samples.tsv` with your sample names, FASTQ paths and FASTA paths, following the [input requirements](docs/inputs.en.md). Then:

```bash
python run_batch.py --samples samples.tsv --check-only
python run_batch.py --samples samples.tsv --out /path/to/your/results/run01 --partition YOUR_CPU_PARTITION
```

Replace `YOUR_CPU_PARTITION` and the output path with real values. A CPU partition must be explicitly provided for a new batch; this avoids accidentally using a cluster's GPU default. The same partition is used for analysis and report jobs. Default resources are **up to 4 concurrent samples, 12 CPUs, 16 GB and 8 hours per sample**. Sort uses 2 additional threads and 256 MB per thread independently of mapping; the first A-to-G conversion step remains single-threaded. Submission can run from a login node; if local resource policy blocks dependency imports, use a small CPU interactive allocation as described in the guide.

## Results and scope

- `report.html`: offline overview with links to QC and logs; explanatory text is currently Chinese.
- `qc_summary.tsv`: read counts, retention, Q30, mapping, orientation, status and errors.
- `all_A_sites.tsv`: every positive-strand A in each original reference, including observed counts before/after A-cutoff, ratios, and missing-observation reasons.
- `significant_sites.tsv`: original per-sample GLORI FDR calls combined with sample names.

`Ratio = A/(A+G)` and `Percent = 100 * Ratio`; neither is automatically background/PCR corrected. Missing observations remain `NA`, not 0%. FDR is applied to the core's prefiltered candidates within each sample, never recalculated across the batch. A header-only significant table is valid only when the sample succeeded. Coordinates are 1-based within the amplicon FASTA, not inferred genomic coordinates.

Samples run in isolated directories. A failed sample does not stop the others; `--retry-failed` creates new attempts only for failed/interrupted samples. Successful validated runs keep cleaned reads, BAM/index, references, counts, QC and logs while removing selected large intermediates. New batches refuse existing output directories.

No paired-end mapping, automatic UMI extraction/deduplication, or reverse-read quantification is provided. No experimental data, target sequences, results, or downloadable example datasets are distributed. Tests generate tiny inputs only in temporary directories.

## Validation and support

```bash
RUN_GLORI_INTEGRATION=1 python -m unittest discover -s tests -v
```

The [Linux workflow](https://github.com/luckingclark/glori-amplicon-batch/actions) creates the documented environment and runs unit tests plus a real-tool integration test comparing the wrapper with the core on the same cleaned synthetic reads. Slurm tests use a mocked scheduler; GitHub CI does not prove operation on a particular production cluster or biological accuracy. See the actual Actions result for the tested commit rather than assuming it passed.

Open an issue with versions and a minimal error description, removing private paths, sample identifiers and data. Update both guides when changing commands or inputs.

## Citation and provenance

For this adaptation, cite **PKU-Gaolab, Ming-Ao Lu. GLORI Amplicon Batch**, with the [repository URL](https://github.com/luckingclark/glori-amplicon-batch) and commit or release used. Machine-readable citation information is available in [CITATION.cff](CITATION.cff).

Please also cite the original methods and software below. Cite GLORI 3.0 when that experimental method applies; using this software alone does not establish which wet-lab method generated the data.

1. Liu, C., Sun, H., Yi, Y., et al. (2023). **Absolute quantification of single-base m6A methylation in the mammalian transcriptome using GLORI.** *Nature Biotechnology*, 41, 355–366. https://doi.org/10.1038/s41587-022-01487-9 . First published online in 2022; the volume year is 2023.
2. Sun, H., Lu, B., Zhang, Z., et al. (2025). **Mild and ultrafast GLORI enables absolute quantification of m6A methylome from low-input samples.** *Nature Methods*, 22, 1226–1236. https://doi.org/10.1038/s41592-025-02680-9 . Cite this when the GLORI 3.0 experimental method is used; software use alone does not establish which wet-lab method generated a dataset.
3. Lu, B. (2024). **Codes for “Mild and ultrafast GLORI enables absolute quantification of m6A methylome from low-input samples”.** Zenodo. https://doi.org/10.5281/zenodo.14233421 . This is a separate collection of companion analysis scripts, not a release of this amplicon adaptation.
4. **GLORI-tools**, by Cong Liu and contributors. https://github.com/liucongcas/GLORI-tools . Its upstream software citation points to https://doi.org/10.5281/zenodo.7014168 . This repository is not asserted to be an exact copy of that archived version.
5. **RNA-m5C**, SYSU-zhanglab and Jianheng Liu. https://github.com/SYSU-zhanglab/RNA-m5C . This is an upstream code source acknowledged by GLORI-tools.


Cite analysis dependencies as appropriate for your work. Code sources and the scope of local changes are described in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## License

Original additions, modifications, and documentation by PKU-Gaolab, Ming-Ao Lu are licensed under the [MIT License](LICENSE). Portions derived from GLORI-tools and RNA-m5C retain their original copyright and MIT license notices in [LICENSES/](LICENSES/); see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for their scope. This independent adaptation does not imply collaboration with or endorsement by the upstream authors. Third-party dependency licenses remain applicable. Zenodo companion scripts and paper materials are not bundled.
