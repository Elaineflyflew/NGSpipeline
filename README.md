# NGSpipeline

NGS 分析流程集合。每个流程是一个**自包含**的子目录：各自的入口、配置、依赖安装、
参考数据目录和测试，互不干扰。加一个新流程不需要改动已有流程的任何文件。

## 流程

| 流程 | 适用数据 | 做什么 | 状态 |
|---|---|---|---|
| [germline-wes-wgs](pipelines/germline-wes-wgs/) | 单样本 WES / WGS | 原始 FASTQ → 质控 → 比对 → 胚系变异检测（bcftools / GATK / DeepVariant 三选或多选）→ 注释（VEP / ANNOVAR / Funcotator）→ 汇总报告 + IGV 批处理 | 可用 |

## 目录约定

每个 `pipelines/<name>/` 保持同一套骨架，看一个流程就能猜出另一个：

```
pipelines/<name>/
├── bin/<name>        入口，只做参数解析与阶段调度
├── src/              业务模块，一个阶段一个文件
├── config/           默认配置（唯一配置源）
├── deploy/           软件体检 + 参考数据下载
├── reference/        参考基因组 / 已知位点 / 注释库（不进 git）
├── extratools/       非 conda 安装的外部工具（不进 git）
├── env/              conda 环境定义
├── tests/            冒烟测试（不依赖生信工具）
├── Makefile
└── README.md         该流程的完整文档
```

这样划分的实际好处：

- **加流程不碰别人**：新流程只在自己的目录里长，根目录只多一行索引。
- **单个流程可整体拷走**：`cp -r pipelines/<name> /somewhere` 就能在别的机器上用。
- **配置只有一处**：路径、线程、阈值都在 `config/`，脚本本身不含硬编码路径。

## 使用

```bash
git clone git@github.com:Elaineflyflew/NGSpipeline.git
cd NGSpipeline/pipelines/germline-wes-wgs

# 1. 先看机器上已经装了什么（不改动任何东西）
bash deploy/deploy.sh doctor

# 2. 只补缺失的软件；全都有就跳过
bash deploy/deploy.sh env

# 3. 下载参考数据（约 40G，不含 VEP cache）
bash deploy/deploy.sh ref bundle

# 4. 跑
bin/germline-wes-wgs -1 raw/Case1_R1.fastq.gz -2 raw/Case1_R2.fastq.gz \
                     -b exome.bed -c bcftools -a vep -t 16
```

各流程的详细用法、参数表、常见坑都在它自己的 `README.md` 里。

## 新增流程

1. 建骨架：`mkdir -p pipelines/<name>/{bin,src,config,deploy,env,tests}`
2. 复用公共层：`src/00_common.sh`（日志 / run / need / 阶段标记）和
   `src/10_cli.sh`（参数解析）在流程之间逻辑一致，直接照 `germline-wes-wgs` 改。
3. 在本文的**流程**表里加一行。
4. `reference/` 与 `extratools/` 由 `.gitignore` 兜底忽略，不用手动排除。

## 关于数据

参考数据与外部工具动辄几十 G，全部不进版本库。根目录 `.gitignore` 用 `**/` 兜底
忽略所有流程下的 `reference/genome|gatk_bundle|vep|funcotator` 和
`extratools/{deepvariant,annovar,bin}`，换机器时用该流程的 `deploy/` 重新拉取。

只克隆代码的话不到 1 MB，几秒完成。
