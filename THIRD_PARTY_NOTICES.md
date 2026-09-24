# Third-party notices

This independent, unofficial amplicon adaptation is maintained by **PKU-Gaolab, Ming-Ao Lu**. Attribution to upstream authors does not imply collaboration, endorsement, or their involvement in maintaining this repository.

## Code and license scope

| Source or contribution | Scope | License notice |
|---|---|---|
| [GLORI-tools](https://github.com/liucongcas/GLORI-tools), Cong Liu | Adapted helper scripts in `core/pipelines/` and upstream portions incorporated into the amplicon workflow. The two entry points in `core/` are local adaptations, not asserted to be upstream files. | [Original MIT license, Copyright (c) 2022 Cong Liu](LICENSES/GLORI-tools-MIT.txt) |
| [RNA-m5C](https://github.com/SYSU-zhanglab/RNA-m5C), SYSU-zhanglab, Jianheng Liu | An upstream code source acknowledged by GLORI-tools; this does not assert direct copying into every local file. | [Original MIT license, Copyright (c) 2019 SYSU-zhanglab, Jianheng Liu](LICENSES/RNA-m5C-MIT.txt) |
| PKU-Gaolab, Ming-Ao Lu | Original local additions, modifications, documentation, environment configuration and tests; no ownership claim over upstream work. | [MIT license, Copyright (c) 2026 PKU-Gaolab, Ming-Ao Lu](LICENSE) |

Both upstream license files are retained in full. Preserve them and any existing source attribution when redistributing the relevant code. Licenses of separately installed dependencies remain applicable.

## Adaptation and changes

The adaptation uses custom amplicon FASTA and single-base annotation while retaining the upstream A-to-G mapping, A restoration and site-calling approach. Changes include current SciPy/statsmodels interfaces, batched conversion writes, independent sorting memory/thread settings, and error propagation. Input safeguards reject duplicate read identifiers, unsupported reverse alignments and nonempty output directories; zero calls produce a header-only table. All amplicons use the requested pileup depth, including names containing `GL`, which previously triggered a special case. These are adaptation and implementation changes, not a claim to have invented GLORI.

The batch layer adds isolated Slurm jobs, checked QC/trimming, per-site count exports, missing-value handling, retained code snapshots and failure-aware aggregation. A guard in the bundled runner permits a legitimate empty per-contig conversion-rate file without fabricating a rate. The analysis core is derived from the same-maintainer [GLORI Amplicon adaptation](https://github.com/luckingclark/glori-amplicon).

The historical upstream commit of the supplied local copy is unknown; no exact release equivalence is claimed.

## Referenced methods and materials

The [GLORI 3.0 companion code](https://doi.org/10.5281/zenodo.14233421) is separately licensed under CC BY 4.0 and was consulted as a method/workflow reference. Its files are not bundled or relabeled MIT. Paper PDFs, figures and study data are not distributed. Full method and software references are in [README.md](README.md#author-licensing-and-citation).
