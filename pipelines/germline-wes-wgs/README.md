# germline-wes-wgs

单样本 WES/WGS germline 变异检测流水线。一条命令从原始双端 FASTQ 走到带注释的 VCF，
三个 caller、三个注释器可自由组合，参考数据和外部工具统一由 `deploy/` 管理。

> 本流程是 [NGSpipeline](../../README.md) 仓库中的一个流程，位于 `pipelines/germline-wes-wgs/`。
> 下面所有命令都在这个目录下执行：`cd pipelines/germline-wes-wgs`。

```
FASTQ ──▶ 质控 ──▶ 修剪 ──▶ 比对 ──▶ BAM 精修 ──▶ 变异检测 ──▶ 注释 ──▶ 报告 + IGV
         FastQC   fastp   bwa mem   fixmate/markdup   bcftools │ GATK │ DeepVariant
                                                                  VEP │ ANNOVAR │ Funcotator
```

## 目录结构

```
germline-wes-wgs/
├── bin/germline-wes-wgs           主程序（只做调度，约 90 行）
├── src/                           业务模块，一个阶段一个文件
│   ├── 00_common.sh               日志、run、need、阶段标记
│   ├── 10_cli.sh                  参数定义与解析
│   ├── 20_context.sh              校验、路径派生、目录创建
│   ├── 30_stage_qc.sh             ① FastQC
│   ├── 31_stage_trim.sh           ② fastp
│   ├── 32_stage_align.sh          ③ bwa mem
│   ├── 33_stage_bam.sh            ④ fixmate → sort → markdup
│   ├── 40_call_bcftools.sh        ⑤-A
│   ├── 41_call_gatk.sh            ⑤-B
│   ├── 42_call_deepvariant.sh     ⑤-C
│   ├── 43_stage_call.sh           ⑤ 调度
│   ├── 50_annotate_vep.sh         ⑥-A
│   ├── 51_annotate_annovar.sh     ⑥-B
│   ├── 52_annotate_funcotator.sh  ⑥-C
│   ├── 53_stage_annotate.sh       ⑥ 调度
│   ├── 60_stage_report.sh         ⑦ 汇总 + IGV 批处理
│   └── 70_preflight.sh            预检（工具/索引/命名一致性）
├── config/pipeline.conf           默认值（用户主要改这里）
├── deploy/                        依赖安装与下载，一个来源一个脚本
│   ├── deploy.sh                  总入口（doctor / status / env / 单项 / all）
│   ├── 00_doctor.sh               软件体检：已装了什么、还缺什么
│   ├── lib.sh                     dl / dl_tar / 空间检查
│   ├── 10_ref_genome.sh           参考基因组 + 索引
│   ├── 20_gatk_bundle.sh          dbsnp / hapmap / Mills / 1000G
│   ├── 30_vep_cache.sh            VEP cache（自动匹配版本）
│   ├── 40_annovar.sh              ANNOVAR + humandb
│   ├── 50_deepvariant.sh          DeepVariant sif
│   └── 60_funcotator.sh           Funcotator 数据源
├── reference/                     参考数据（大文件，不进 git）
│   ├── genome/                    fasta + .fai + .dict + bwa index
│   ├── gatk_bundle/               dbsnp_146.hg38 等
│   ├── vep/                       homo_sapiens/<ver>_GRCh38/
│   ├── funcotator/                数据源目录
│   └── bed/                       捕获区间（放进去就自动用）
├── extratools/                    外部工具（非 conda 安装的）
│   ├── deepvariant/               deepvariant.sif
│   ├── annovar/                   annovar + humandb
│   └── bin/                       自己编译的小工具，自动加进 PATH
├── env/environment.yaml           conda 环境定义
├── tests/smoke_test.sh            冒烟测试（不需要生信工具）
├── Makefile
└── README.md
```

## 快速开始

```bash
# 1. 先体检：机器上已经装了什么（不改动任何东西）
bash deploy/deploy.sh doctor

# 2. 只补缺失的软件。已经有的直接复用，全都有就跳过
bash deploy/deploy.sh env
conda activate wes          # 只在第 2 步真的装了东西时才需要

# 3. 装参考数据
bash deploy/deploy.sh status
bash deploy/deploy.sh ref bundle        # 最小可用：bcftools + GATK 路线
bash deploy/deploy.sh vep               # 需要 VEP 注释时再加

# 4. 跑
bin/germline-wes-wgs -1 raw/Case1_R1.fastq.gz -2 raw/Case1_R2.fastq.gz \
                 -c bcftools -a vep -t 16

# 5. 看结果
less results/07_report/summary.txt
igv --batch results/07_report/igv.bat
```

## 已有环境怎么办（不想重装 GATK）

**不用重装。** `deploy.sh env` 先做体检，只装缺的；`bin/germline-wes-wgs` 运行时也会
自己去找，找到就直接用。

探测范围：

```
PATH  →  extratools/bin  →  $CONDA_PREFIX/bin  →  常见环境目录
  /home/biotools/miniforge3/envs/*/bin
  /opt/conda/envs/*/bin
  ~/miniforge3|mambaforge|miniconda3|anaconda3/envs/*/bin
```

所以 GATK 装在 `/home/biotools/miniforge3/envs/gatk/bin` 但从来没 activate 过，
也能被找到。pipeline 会把那个目录临时加进 PATH，并在日志里写明来源：

```
[22:09:47] 找到 gatk: /home/biotools/miniforge3/envs/gatk/bin/gatk
[22:09:47] 已加入 PATH: /home/biotools/miniforge3/envs/gatk/bin
```

`doctor` 的输出示例：

```
工具                 状态     位置 / 用途
--------------------------------------------------------------------------
fastqc               OK      /usr/local/bin/fastqc  (FastQC 原始数据质控)
gatk                 OK      /home/biotools/miniforge3/envs/gatk/bin/gatk  (GATK 路线)
vep                  缺失(可选) -  (VEP 注释)

必需工具: 7/8 就绪
缺失的必需工具: bcftools
建议: bash deploy/deploy.sh env     # 只装缺失的，已有的不动
```

`status` 的软件行用同一套探测，**不要求存在名为 `wes` 的 conda 环境**：

```
软件                 OK         8/8 必需工具就绪（复用本机已装）
```

工具在 `/bioapps/bin`、在别人的 env 里、在 `~/.local/bin`，都算就绪。只有真的
缺必需工具时才会报 `缺失(1)  bcftools → bash deploy/deploy.sh env`。

`env` 的行为：

| 情况 | 动作 |
|---|---|
| 必需工具全都有 | 什么都不装，直接跳过 |
| 缺几个 | 生成 `env/environment.missing.yaml`，只装这几个 |
| `wes` 环境已存在 | 用 `conda install -n wes` 往里补，不重建 |
| 加了 `--with-optional` | 连 multiqc / vep / apptainer 一起装 |

可选工具（multiqc、vep、apptainer、ANNOVAR）默认不装——不用 VEP 注释的样本
没必要下 25 G 的 cache。

## 依赖体积与耗时

| 目标 | 体积 | 下载/构建 | 说明 |
|---|---|---|---|
| `env` | 0-4 G | 0-10 分钟 | 只装缺失的；全都有则跳过 |
| `ref` (broad) | ~8 G | 15 分钟 | 自带索引，不用自己建 |
| `ref --source ensembl` | ~4 G | +40 分钟 | 需要本地 `bwa index` |
| `bundle` | ~1.5 G | 5 分钟 | dbsnp / hapmap / Mills / 1000G |
| `vep` | ~25 G | 20 分钟 | 版本自动匹配本机 vep |
| `funcotator` | ~15 G | 15 分钟 | |
| `deepvariant` | ~5 G | 20 分钟 | singularity 镜像 |
| `annovar` | ~10 G+ | 手动 | 需学术注册，脚本只打印指引 |

全装约 **70 G**。只跑 bcftools + VEP 的话 40 G 够。

## 常用命令

```bash
# 三条路线一起跑，自动出一致性对比
bin/germline-wes-wgs --bam results/03_map/Case1.mdup.bam -c all -a none

# 复用已有 BAM，GATK + Funcotator
bin/germline-wes-wgs --bam Case1.mdup.bam -c gatk -a funcotator

# 只重跑变异检测（改完阈值之后）
bin/germline-wes-wgs --bam Case1.mdup.bam -c all --stages call,report --force

# 上服务器前先看会跑什么
bin/germline-wes-wgs ... --dry-run

# 只查环境和数据库，不跑数据
bin/germline-wes-wgs ... --check

# 看所有阶段
bin/germline-wes-wgs --list-stages
```

参数优先级：**命令行 > 环境变量 > `config/pipeline.conf` > 代码兜底**。
所以临时改阈值不用动文件：

```bash
THREADS=32 MIN_QUAL=20 bin/germline-wes-wgs -1 ... -2 ...
```

## 配置：改一处就够

`config/pipeline.conf` 是唯一的默认值来源。路径全部指向项目内目录，
所以 `deploy` 之后**不需要改任何路径**——pipeline 会自己找到：

```bash
: "${REF:=$(discover_ref)}"                                  # 按索引识别，不看文件名
: "${BUNDLE:=$PROJECT_ROOT/reference/gatk_bundle}"
: "${VEP_CACHE:=$PROJECT_ROOT/reference/vep}"
: "${DV_SIF:=$(first_existing "$PROJECT_ROOT/extratools/deepvariant/deepvariant.sif")}"
```

参考基因组不要求固定文件名，也不需要软链接：`discover_ref` 扫描
`reference/genome/` 下所有 fasta，按有没有 `.fai`/`.bwt`/`.dict` 打分选最优的。
Broad 版、Ensembl 版、你自己重命名的 `ref.fa` 都能直接用。

想用别处的数据，加一行 `--vep-cache /other/path` 或设环境变量即可。

## 扩展

**加一个 caller**：新建 `src/44_call_mycaller.sh` 定义 `call_mycaller()`，
在 `43_stage_call.sh` 的 case 里加一行，在 `20_context.sh` 的 `CALLERS` 展开里加一个名字。
主程序不用动。

**加一个 annotator**：同理，`src/54_annotate_xxx.sh` + `53_stage_annotate.sh` 一个分支。

**加一个阶段**：`src/3x_stage_xxx.sh` 定义 `stage_xxx()`，
在 `bin/germline-wes-wgs` 的 `run_pipeline()` 里加一行，在 `20_context.sh` 的阶段白名单里加名字。

约定：每个 `stage_*` 函数自己调用 `mark_done <name>`，调度层只负责跳过判断。

## 设计上的几个决定

**1. 为什么用 bash 而不是 Snakemake。** 单样本、单机、跑一次，Snakemake 的
DAG 调度和 `.snakemake` 元数据是净负担；bash 每条命令都是明文，答辩时能逐行指着讲。
要扩到几十个样本批量跑时，`bin/germline-wes-wgs` 可以直接被 Snakemake 的 `shell:`
调用，不需要改任何模块。

**2. 为什么拆成多个文件。** 一个 800 行的单文件里，改 fastp 参数要滚到第 300 行。
现在改修剪只开 `31_stage_trim.sh`（30 行），改 GATK 只开 `41_call_gatk.sh`。
`deploy/` 同理：加一个新的数据源只需要新增一个文件，`deploy.sh` 加一行 case。

**3. 为什么每个阶段单独落 `.done` 标记。** 三个 caller 各自标记，
某个路线失败重跑时不会把已经跑完的另外两条重算。注释是 caller × annotator
的笛卡尔积，也各自标记。

**4. 为什么缺数据库时降级而不是崩。** 注释库没配好就跳过那一步并 `WARN`，
不让跑了四小时的比对结果白费。GATK 缺 bundle 时跳过 BQSR，缺 hapmap 时跳过
Tranche——结果仍然可用，只是精度略低。

## 指标怎么看

跑完 `results/07_report/summary.txt` 里有：

| 指标 | 期望 | 异常说明 |
|---|---|---|
| 比对率 (mapped) | > 95% | 低 → 参考版本不对或污染 |
| properly paired | > 90% | 低 → 文库/接头问题 |
| dup rate | 5-20% (WES) | 过高 → 起始量不足 |
| 目标区平均深度 | > 100x (WES) | 低 → 捕获效率差 |
| Ti/Tv (SNP) | 2.0-2.1 | < 1.8 → 假阳性多 |
| het/hom | 1.5-2.0 | 偏离 → 污染或近亲 |
| 三线 SNP 一致率 | > 80% | indel 一致率天然更低 |

## 常见问题

**参考命名不一致。** Ensembl 的 `17` vs GATK bundle 的 `chr17`。
这个坑不报错，但 BQSR 找不到染色体、`isec` 交集为 0。预检会检查并告警。
解决：`deploy/deploy.sh ref --source broad`（默认就是），或 `--rename-chr`。

**`markdup` 的 dup rate 偏低。** 检查 `fixmate -m` 前面有没有 `sort -n`。
fixmate 需要 name-collated 输入，顺序错了不报错但结果不准。

**`bcftools mpileup` 的 AF 不对。** 加 `-d 1000`。默认 250 会削掉高深度位点。

**bundle 是解压的 `.vcf`，要不要压回去？** 不用。`bundle_file` 两种格式都认
（`.vcf` 配 `.idx`，`.vcf.gz` 配 `.tbi`），两种都在时优先 `.vcf.gz`。
未压缩的缺 `.idx` 时会在跑 GATK 前自动 `gatk IndexFeatureFile` 补建。

**VEP 报找不到 cache。** cache 版本必须和 vep 主版本一致。
`deploy/30_vep_cache.sh` 会自动探测，也可以 `--release 110` 指定。

**服务器上 `nohup` 后工具全找不到。** 现在不需要了：pipeline 启动时会自己
探测 `/bioapps/bin`、`$CONDA_PREFIX` 和常见 conda 环境目录，找到就把路径加进去。
只有当工具在某个必须 activate 才生效的私有环境里、且不在探测范围内时，才需要
`conda activate` 之后再 `nohup`。

## 测试

```bash
bash tests/smoke_test.sh      # 57 项，不需要任何生信工具
```

只验证脚手架：参数解析、路径派生、阶段跳过、dry-run 不落标记、deploy 的 dry-run。
真实数据必须上服务器跑。

## 与旧版的关系

`lab5_ngs/` 下的 `run_pipeline.sh` / `lab5_onefile.sh` / `ngs_pipeline.sh`
是早期版本，逻辑已全部并入本项目，保留作对照。新项目只依赖 `bin/germline-wes-wgs`。
