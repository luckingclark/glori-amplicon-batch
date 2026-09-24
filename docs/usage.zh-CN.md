# 安装与批量分析

[首页](../README.md) · [English](usage.en.md) · [输入文件规范](inputs.zh-CN.md)

本指南从一个可以使用 Conda 的集群新账户开始。所有命令都在 Linux 的 Bash 终端执行。只复制代码框中的命令，不要复制终端提示符。运行前，把所有 `/path/to/your/...` 和 `YOUR_CPU_PARTITION` 占位符替换成自己的真实路径和分区名。

## 1. 准备代码与输入文件

在 GitHub 点击 **Code → Download ZIP**，解压后把完整代码文件夹放到自己的集群账户下。Python 文件、`core/`、`LICENSES/`、许可证和 `environment.yml` 要保留在一起。如果安装了 Git 且账户有访问权限，也可以克隆仓库。

除了代码，只需要上传已经按样品拆分的 FASTQ，以及相应的原始扩增子 FASTA。每个样品选择一个 read end，通常是 R1。不需要 R2、样本 index reads、全基因组索引、注释数据库、以前的结果或本地 Conda 环境。本仓库不提供示例数据集。填写样品表前，请先看[输入文件规范](inputs.zh-CN.md)。

进入代码目录并查看文件：

```bash
cd /path/to/your/glori-amplicon-batch
pwd
ls
```

`cd` 切换目录；`pwd` 显示当前路径；`ls` 列出文件。应能看到 `run_batch.py`、`environment.yml`、`samples.tsv` 和 `core`。

代码、输入和输出都必须位于登录节点和计算节点共同可访问的存储中。不要使用只有登录节点能访问的临时目录。

## 2. 建立自己的 Conda 环境

如果提示找不到 `conda`，先加载集群提供的 Conda module，或按管理员的说明配置 Conda。然后在自己的账户下创建环境：

```bash
mkdir -p "$HOME/.conda/envs" "$HOME/.conda/pkgs"
conda config --prepend envs_dirs "$HOME/.conda/envs"
conda config --prepend pkgs_dirs "$HOME/.conda/pkgs"
export CONDA_CHANNEL_PRIORITY=strict
conda env create --file environment.yml
conda activate glori_amplicon
python --version
```

`mkdir -p` 创建个人环境目录和安装包缓存目录，已经存在时不会报错。两条 `conda config` 命令把它们设为个人账户的优先目录，此设置在以后登录时仍有效，不会修改共享的软件安装。`export` 为当前终端启用严格的 channel 优先级。`conda env create` 根据 `environment.yml` 下载依赖并创建命名环境，可能需要几分钟，也需要访问 Conda 软件源。`conda activate` 激活环境；`python --version` 应显示 Python 3.10。

环境包含 Python、NumPy、pandas、SciPy、statsmodels、Biopython、pysam，以及 Bowtie 1、samtools、Trim Galore、Cutadapt、FastQC。Slurm 命令由集群提供，不在这个 Conda 环境里安装。

环境只需创建一次。以后重新登录，只需要激活。如果已经有 `glori_amplicon`，先确认它是否符合本环境规格，不要为了重走教程删除正常工作的环境。也可以用 `conda env create --name YOUR_NEW_ENV --file environment.yml` 创建另一个命名环境，随后激活相应名称。脚本检查 Python 是否属于当前激活的环境，不限制环境一定叫什么名字。

## 3. 填写 samples.tsv

`samples.tsv` 初始只有表头。用文本编辑器打开，按[输入文件规范](inputs.zh-CN.md)每个样品添加一行。三列之间必须使用 **Tab**，不能用空格或逗号。不同样品可以使用不同参考，也可以共用一个参考。

```bash
cat samples.tsv
```

`cat` 显示文件内容，方便核对样品名和路径。相对输入路径以 `samples.tsv` 所在目录为起点解析，不以输出目录为起点。样品名忽略大小写后也不能重复。不要把教程中的占位路径当成实际文件路径。

## 4. 检查输入和依赖

通常直接留在登录节点，激活环境并执行：

```bash
conda activate glori_amplicon
python run_batch.py --samples samples.tsv --check-only
```

`--samples` 指定样品表。`--check-only` 检查文件路径、所有参考序列、FASTQ 第一条记录、Python 依赖、软件和 Slurm 命令，但不会提交作业。完整 FASTQ 扫描和 FastQC 放在后面的样品计算作业中。脚本会在导入科学计算包前限制数值库线程数，减少不必要的 OpenBLAS 内存申请。

**一般不需要先进入 CPU 节点。** 如果登录节点仍然限制导入或内存申请，或者集群规定必须这样做，可以申请少量交互式 CPU 资源后再检查：

```bash
srun --partition=YOUR_CPU_PARTITION --cpus-per-task=1 --mem=2G --time=00:20:00 --pty bash
conda activate glori_amplicon
cd /path/to/your/glori-amplicon-batch
python run_batch.py --samples samples.tsv --check-only
exit
```

`srun` 申请临时计算资源，可能需要排队。请替换成真实 CPU 分区名。`exit` 释放这次交互资源并返回登录终端。如果提交时在登录节点也会报同样错误，应在 `exit` 之前执行第 5 步，因为提交本身会重复检查依赖。批量作业会另外申请资源，交互式申请的 1 CPU 不负责提供四个样品所需的 48 CPU。如果管理员禁止从计算节点提交作业，则按集群支持的方式提交。

在没有分析依赖的机器上，包括 Windows，可以用 `python run_batch.py --samples samples.tsv --validate-inputs` 仅检查输入。这不能证明集群环境已经可用。

## 5. 一次提交全部样品

选择一个尚不存在的输出目录，并明确指定 CPU 分区：

```bash
python run_batch.py --samples samples.tsv --out /path/to/your/results/run01 --partition YOUR_CPU_PARTITION
```

`--out` 指定脚本将创建的结果目录。首次提交时，此目录必须尚不存在。新批次必须填写 `--partition`，样品数组作业和汇总作业都会使用这个分区。脚本不会默默使用集群默认分区。实际分配到哪台计算节点由 Slurm 决定，与从哪个目录、哪台主机提交无直接绑定关系。

脚本会输出样品数组作业和报告作业的编号。它先保存软件快照、输入和资源记录，提交暂挂的样品数组，再提交依赖该数组结束的汇总作业，确认两个编号保存后才释放数组开始排队。成功提交后可以关闭终端。

默认申请的资源：

| 资源 | 默认值 |
|---|---|
| 同时运行的样品数 | 最多 4 个 |
| 每个样品 CPU | 12 |
| 每个样品内存 | 16 GB |
| 每个样品时限 | 8 小时 |
| 汇总作业 | 1 CPU、2 GB、1 小时 |

四个样品同时运行时合计最多申请 48 CPU、64 GB 内存，但不强制放在同一节点。实际排队还取决于可用内存、分区规则及其他作业。排序单独使用 2 个额外线程、每线程 256 MB 的内存限额，与 Bowtie 线程设置分开。最初的 reads A→G 转换仍然是单线程，提高 `--cpus` 不能加速这一步。

如有需要，可在提交另一个新目录时调整资源：

```bash
python run_batch.py --samples samples.tsv --out /path/to/your/results/run02 --partition YOUR_CPU_PARTITION --max-parallel 4 --cpus 8 --memory 16G --time 12:00:00
```

此时四个样品最多使用 32 CPU，每样品允许 12 小时。`--cpus` 至少为 3。如果集群要求计费账户，再添加 `--account YOUR_SLURM_ACCOUNT`。提交时添加 `--keep-all` 可保留原本由批量脚本在成功导出后清理的大型中间文件，但不会恢复分析核心内部已清理的文件。

## 6. 查看运行进度

```bash
squeue -u "$USER"
cat /path/to/your/results/run01/qc_summary.tsv
ls /path/to/your/results/run01/logs
```

`squeue` 显示自己的活跃作业、分区、状态及节点。汇总作业显示 `Dependency` 等待属于正常现象。样品日志为 `logs/array-<job>_<index>.log`，汇总日志为 `logs/report-<job>.log`。初始汇总表将样品标为排队状态；汇总文件在数组结束后由报告作业重新生成，并非实时刷新。运行过程中可查看各次尝试的 `status.json` 和 `logs/pipeline.log`。

每个样品自动执行：完整 FASTQ 检查与计数 → 原始 FastQC → Illumina 接头及低质量末端清理 → 清洗后 FASTQ 检查与 FastQC → 建立参考及索引 → GLORI 转换、比对、恢复 A、排序、pileup 和位点检测 → BAM 与计数检查 → 位点比例导出。FastQC 的 WARN/FAIL 会记录下来，不会仅凭这些标记把已正常完成的分析视为软件运行失败。

## 7. 阅读结果

把**整个结果目录**下载到自己的电脑，用浏览器打开 `report.html`。保留目录结构，才能正常打开其中的 FastQC 和日志链接。报告目前使用中文说明，表格列名为固定英文。

| 文件 | 内容 |
|---|---|
| `report.html` | 离线 QC 概览、样品状态和报告链接 |
| `qc_summary.tsv` | reads 数、保留率、Q30、比对率、反向 reads、显著位点数、警告与失败信息 |
| `all_A_sites.tsv` | 每个样品原始参考中的所有 A，包括缺失观测的位置 |
| `significant_sites.tsv` | 成功样品的显著位点，保留各样品原本的 FDR 值 |
| `samples/<sample>/attempt_001/` | 单样品 reads、参考、BAM/索引、计数、QC、日志和状态 |
| `software/`、`rounds/`、`batch.json` | 软件快照及运行、提交记录 |

`all_A_sites.tsv` 的主要列：

| 列名 | 含义 |
|---|---|
| `Sample`、`Chr`、`Sites`、`Strand` | 样品、FASTA 记录名称、从 1 开始的扩增子坐标、`+` 链 |
| `Acov`、`Gcov`、`AGcov` | A-cutoff 过滤后的 A、G、A+G 计数 |
| `Ratio`、`Percent` | A/(A+G)，及其乘以 100 后的百分比 |
| `*_all` | A-cutoff 过滤前的计数与比例，仍然受比对、pileup 过滤及深度上限影响 |
| `Signal` | A-cutoff 过滤后保留的总 pileup 观测比例 |
| `Background_CR` | 核心流程报告的扩增子级 A→G 转换率 |
| `FDR_status` | `significant` 表示通过显著性筛选；`not_reported` 不区分具体未通过原因 |
| `Observation_status` | 已观测、无可用记录、过滤后无 A/G，或样品失败信息 |

默认 A-cutoff 过滤后的 m6A 百分比估计看 `Percent`。它不会自动校正转换背景、PCR 偏差等实验因素。`NA` 表示没有可用估计，不能当成 0% 修饰。有比例值不代表一定显著。失败样品会给出缺失估计，不会把部分结果冒充有效数据。解释结果前先检查 `Status` 和 QC。

如果样品成功完成但没有显著位点，显著位点表只有表头是允许的。核心流程在每个样品的预筛选候选位点中进行 FDR 校正，汇总时不会再跨样品校正。坐标是相对于输入 FASTA 的坐标，不会自动换算为基因组坐标。

分析默认参数：Illumina 接头、Q20、最短长度 25、接头最短重叠 5、错误率 0.1；Bowtie 1 最多 2 个错配且只保留唯一比对；pileup 最大深度 10,000；A-cutoff 为 3；gene/扩增子背景；最低 A+G 覆盖 15、A 数 5、比例 0.1、signal 0.8、A+G/total 0.8；二项检验 P 阈值 0.005，FDR 阈值 0.005。由于 pileup 深度上限，深度很高的样品位点计数不一定等于全部比对 reads 数。

流程不做分子 UMI 去重、双端比对或引物剪切。请确认参考描述的是实际测序 insert，并使用预期方向。遇到不支持的反向比对会明确报错。默认参数在生物学上是否适用于实验，需要结合实验设计判断。

## 8. 重跑失败样品或重新汇总

等样品数组及对应报告作业都从 `squeue` 中消失，保持原始输入不变，查看报错并修正资源或环境问题后，仅重跑失败或中断的样品：

```bash
python run_batch.py --samples samples.tsv --out /path/to/your/results/run01 --partition YOUR_CPU_PARTITION --retry-failed
```

成功样品会跳过；失败样品创建新的尝试目录，原尝试和早前报告副本保留。如果要增加内存或时限，可以添加 `--memory 32G --time 16:00:00`。如果改变了数据、参考或分析软件，应换一个新的输出目录重新提交。重试沿用原批次的软件快照，因此仅更新启动脚本旁边的代码不会更新现有批次使用的分析代码。

如果自动汇总作业失败，但样品作业都已结束，可重新生成报告：

```bash
python run_batch.py --out /path/to/your/results/run01 --summarize
```

这不会重跑分析，且记录中的作业仍活跃时会拒绝执行。如果提交超时导致作业编号不确定，先核对 `rounds/.../submit.log` 与 Slurm 队列，再决定如何处理，避免重复作业。不要通过删除执行记录来绕过保护。

## 9. 软件测试

开发者可在装好完整环境的 Linux 机器上执行：

```bash
RUN_GLORI_INTEGRATION=1 python -m unittest discover -s tests -v
```

单元测试覆盖输入检查、缺失值、QC 导出、失败样品隔离、提交与重试及资源处理。真实工具集成测试会临时生成少量 reads，运行清洗、QC 和分析，并与直接运行核心流程的结果比较。测试中的 Slurm 是模拟的，不代表某个具体集群已验证，也不能证明生物学准确性。实际 Linux 结果见 [GitHub Actions](https://github.com/luckingclark/glori-amplicon-batch/actions)。
