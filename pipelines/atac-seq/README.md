# atac-seq

单样本 ATAC-seq 开放染色质分析流水线。一条命令从原始双端 FASTQ 走到带自动结论的
汇总报告，15 个步骤各自产出可核查的中间文件，最后一步把散落的指标串成
「这份数据能不能用」的判断。

> 本流程是 [NGSpipeline](../../README.md) 仓库中的一个流程，位于 `pipelines/atac-seq/`。
> 下面所有命令都在这个目录下执行：`cd pipelines/atac-seq`。

![流程结构](docs/pipeline.svg)

## 快速开始

```bash
# 1. 看机器上已经装了什么（不改动任何东西）
bash deploy/deploy.sh doctor

# 2. 看参考数据都指到哪了（本流程默认不下载任何数据）
bash deploy/deploy.sh status

# 3. 跑
bin/atac-seq \
  -1 /public/storage/shaojf/Course/NGS/wkdir4labs/lab9/rawdata/sub20M.SRR13849313_1.fastq.gz \
  -2 /public/storage/shaojf/Course/NGS/wkdir4labs/lab9/rawdata/sub20M.SRR13849313_2.fastq.gz \
  -t 20

# 4. 看结论
less results/10_report/summary.txt
igv --batch results/10_report/igv.bat
```

跑之前想先看会执行什么：加 `--dry-run`，只打印命令不执行。

## 十五步做什么

| 步 | 工具 | 做什么 | 关键点 |
|---|---|---|---|
| ① | FastQC | 原始质量、GC、接头、重复 | Q30 和接头比例决定 ② 的力度 |
| ② | fastp | 去接头 + 质量过滤（Q15，短于 36 bp 丢） | SRA 数据多数已清洗，可跳过 |
| ③ | bowtie2 | `--local --very-fast-local --no-mixed --no-discordant` | `--local` 允许末端错配，ATAC 必须；`--no-mixed/discordant` 只留正常配对 |
| ④ | samtools | `view -q 20` + 去 chrM → sort → index | **不能 `grep -v chrM`**，会连 header 的 `@SQ SN:chrM` 一起删掉 |
| ⑤ | alignmentSieve | `--ATACshift`，正链 +4 / 负链 −5 | Tn5 留下 9 bp 缺口，不校正会让 footprint 偏 4-5 bp |
| ⑥ | samtools | 重排 + 建索引 | alignmentSieve 只保留前 11 列且不再有序，**必须重排** |
| ⑦ | awk | 按 TLEN 分档：NFR / mono / multi | 必须用**已偏移**的 srt2，否则 147 bp 阈值错位 9 bp |
| ⑧ | MACS2 | `callpeak -f BAMPE --nomodel -g hs` × 3（all / NFR / mono） | 不加 `-c`（无 input）；BAMPE 下 `--shift` 会被静默丢弃 |
| ⑨ | bedtools | `intersect -v` 去 ENCODE blacklist | 三份 peak 各自过滤 |
| ⑩ | bedtools + deeptools | FRiP 两种算法交叉验证 + TSS 富集分数 | 两算法差 >5% 说明口径不一致 |
| ⑪ | deeptools | bamCoverage → bigWig（CPM）+ TSS 曲线 + fingerprint | `--ignoreDuplicates --blackListFileName --extendReads` 缺一不可 |
| ⑫ | HOMER + bedtools | annotatePeaks 注释分布 + closest 最近基因 | promoter 比例是"像不像 ATAC"的直接证据 |
| ⑬ | HOMER | `findMotifsGenome -mknown` | 不加 `-mknown` 会做从头 discovery，多跑几小时 |
| ⑭ | 内置 awk | 片段长度分布 + SVG 折线图 | 不依赖 R，避免 ATACseqQC 的版本问题 |
| ⑮ | 内置 | 汇总指标 + PASS/WARN/FAIL + 结论 + IGV 批处理 | 阈值全在 `config/pipeline.conf` |

## 输入

| 项 | 怎么给 | 默认 |
|---|---|---|
| 测序数据 | `-1/-2` 双端 FASTQ，或 `--bam` 已偏移的 BAM | 必填二选一 |
| 物种 | `-g hg38` | 读 `config/genomes/<name>.conf` |
| bowtie2 索引 | `--bt2-index` | 自动探测课程目录下的 `lab9/refs/hg38` |
| 参考 fasta | `--fasta` | `lab9/refs/hg38.fa` |
| blacklist | `--blacklist` | `.../ENCODE/blacklist.v2/hg38-blacklist.v2.bed` |
| 基因注释 | `--gtf` | `lab10/refs/gencode.v43.annotation.gtf` |
| motif 库 | `--motif-db` / `--homer-motifs` | JASPAR2022 meme + HOMER known.motifs |

**默认不下载任何数据。** 所有路径都在 `config/genomes/hg38.conf` 里用
`first_existing` 探测多个候选位置（`/public/workspace` 和 `/public/storage` 两套
课程目录都覆盖），按顺序取第一个存在的。想换物种就新增一个
`config/genomes/mm10.conf`，流程代码一行都不用改。

## 输出

```
results/
├── 00_logs/<sample>.<时间戳>.log     完整日志
├── 01_qc/           *_fastqc.html  fastqc.summary.tsv（Q30 / GC / 读长）
├── 02_clean/        <s>_R1.clean.fq.gz  <s>.fastp.{json,html}  <s>.retention.txt
├── 03_map/          <s>_bowtie2.sam  <s>_bowtie2.log  <s>.align_rate.txt  <s>.mito.tsv
├── 04_bamfiles/     <s>.srt.bam（未偏移）  <s>.srt2.bam（已偏移，下游都用它）
│                    <s>.NFR.bam  <s>.mono.bam  <s>.multi.bam
│                    <s>.fragment_class.tsv（三档占比）
├── 05_peaks/        <s>_peaks.narrowPeak  <s>.filtered.peaks.narrowPeak
│                    <s>_NFR_*.narrowPeak   <s>_mono_*.narrowPeak
│                    peak_counts.tsv（三路线峰数 + 基因组覆盖）
├── 06_deeptools/    <s>.tss_profile.pdf  <s>.tss_heatmap.pdf  <s>.fingerprint.pdf
├── 07_bigwigs/      <s>.bamCoverage.bw
├── 07_annotation/   <s>.homer.anno.txt  <s>.homer.annostats.txt  <s>.bedtools.anno.txt
├── 08_motif/        <s>.homer.motif/knownResults.{txt,html}  <s>.toppeaks.narrowPeak
├── 09_atacqc/       <s>.fragment_length.tsv  <s>.fragment_length.svg
│                    <s>.frip.tsv  <s>.tsse.tsv
├── 10_report/       summary.txt   ← 主产物
│                    metrics.tsv  igv.bat  igv.*.png  multiqc/
└── .status/*.done   断点续跑标记
```

**主产物是 `10_report/summary.txt`**，其余都是它的证据链。报告里每一行指标都
能顺着目录名回溯到原始文件。

## 指标怎么看

`summary.txt` 会对每个指标给 PASS / WARN / FAIL，判定阈值在
`config/pipeline.conf` 的「结论判定阈值」一节。这里说清楚每个数意味着什么：

| 指标 | 期望 | 异常说明 |
|---|---|---|
| 比对率 | > 80% | 低 → 参考版本不对、接头残留、或污染 |
| 线粒体占比 | < 5% 优秀，< 20% 可接受 | > 20% → 文库质量差，核基因组信息严重损失 |
| NFR 比例 | > 40% | 低 → Tn5 消化不足，或细胞核完整性差 |
| 片段分布主峰 | 50-100 bp | 右移 → 消化不足；只有单峰且后面全平 → 过度消化 |
| peak 数 | > 10000 | 过少 → 信号弱或比对有问题 |
| FRiP | > 0.2（ENCODE 底线），> 0.3 好 | 低 → peak 不可靠，或 peak 本身是假的 |
| TSS 富集 | ≥ 5（ENCODE 标准） | 低 → 信号没富集在调控区，可能只是背景 |
| promoter 比例 | 30-50% | 过低 → 不像 ATAC；过高 → 可能混了 RNA-seq |
| FRiP 双算法差 | < 5% | 超过 → BED 坐标口径不一致，或一个用了 `-q` 一个没用 |

三份 peak 的相对关系也是一条结论：**NFR 的峰数应多于 mono**。如果反过来，
说明片段分档或偏移校正有问题，而不是生物学发现。

## 常用命令

```bash
# 只重跑 peak calling（改完阈值之后）
bin/atac-seq --bam results/04_bamfiles/Case1.srt2.bam --stages peak,report --force

# 复用已有 BAM，跳过上游
bin/atac-seq --bam Case1.srt2.bam -t 20

# 先看会跑什么
bin/atac-seq -1 R1.fq.gz -2 R2.fq.gz --dry-run

# 只查工具和参考数据，不跑数据
bin/atac-seq -1 R1.fq.gz -2 R2.fq.gz --check

# 跳过耗时步骤
bin/atac-seq -1 R1.fq.gz -2 R2.fq.gz --no-motif          # 跳过 ⑬
bin/atac-seq -1 R1.fq.gz -2 R2.fq.gz --stages qc,trim,align,bam,shift,fraction,peak,report
```

参数优先级：**命令行 > 环境变量 > `config/pipeline.conf` > 代码兜底**。
临时改阈值不用动文件：

```bash
THREADS=32 MIN_MAPQ=30 bin/atac-seq -1 ... -2 ...
```

## 参数

```
输入
  -1, --r1 FILE          双端 R1          -2, --r2 FILE   双端 R2
      --bam FILE         已偏移的 BAM，跳过 qc/trim/align/bam/shift

核心
  -g, --genome NAME      物种配置（默认 hg38）
  -o, --outdir DIR       输出目录（默认 ./results）
  -s, --sample NAME      样本名（默认从文件名推导）
  -t, --threads N        线程数（默认 20）

阈值
      --min-mapq N       比对质量下限（默认 20）
      --keep-mt          保留线粒体 reads
      --nfr-max N        NFR 上限（默认 147）
      --mono-max N       单核小体上限（默认 294）
      --macs2-qvalue X   MACS2 q 值（默认 0.05）
      --toppeaks N       motif 分析取前 N 个峰（默认 10000）

流程控制
      --stages LIST      qc,trim,align,bam,shift,fraction,peak,metrics,bigwig,annotate,atacqc,report
      --list-stages      打印阶段说明
      --force            忽略完成标记
      --check            只检查工具和数据
      --dry-run          只打印命令
      --no-motif         跳过 motif
      --no-atacqc        跳过片段分布
      --meme             额外跑 MEME（很慢）
      --tobias           额外跑 TOBIAS footprint
```

## 和课程三个脚本的对应

| 课程脚本 | 对应步骤 | 本流程改了什么 |
|---|---|---|
| `refer.sh` | ①-⑪ | 变量抽到 `config/`；补了 FastQC 汇总表；chrM 过滤改成只判断 RNAME 列；补了 FRiP 双算法与 TSSE |
| `PeakCalling.sh` | ⑦-⑭ | 分档加注释说明为什么必须用 srt2；三份 peak 各自过滤 blacklist；注释结果统一到 `07_annotation/` |
| `peak.motif.sh` | ⑬⑭ | motif 默认只跑 HOMER（MEME 5 小时默认关）；ATACseqQC 换成内置实现 |

原脚本里的几个坑已经修掉，并在代码注释里写明原因：

- **`grep -v chrM` 会删 header**：`@SQ SN:chrM LN:16569` 也含 `chrM`，
  过滤后 BAM 缺这条染色体的定义。改成 awk 只判断数据行的第 3 列。
- **`multi_nucleosome.bam` 没有样本名前缀**：多样本跑会互相覆盖。
- **`04.bamstats` 目录名带点**：不影响运行但不利于按前缀匹配。

## 设计上的几个决定

**1. 为什么模仿 snakePipes 的三层配置，而不是一个 config。**
snakePipes 把参数拆成 `defaults.yaml`（通用）+ `organisms/*.yaml`（物种）+ 命令行。
好处是换物种只改物种文件，通用参数不动。这里对应
`config/pipeline.conf` + `config/genomes/<name>.conf` + 命令行，同样的分工。

**2. 为什么用 bash 而不是 Snakemake。**
单样本、单机、跑一次，Snakemake 的 DAG 和元数据是净负担；bash 每条命令都是明文，
答辩时能逐行指着讲。要扩到多样本批量跑，`bin/atac-seq` 可以直接被 Snakemake 的
`shell:` 调用，不需要改任何模块。

**3. 为什么每个阶段单独落 `.done` 标记。**
⑬ motif 要跑 2 小时，⑭ 只要 1 分钟。某一步失败重跑时只补那一步，不会把前面
半小时的比对重算一遍。

**4. 为什么缺工具时降级而不是崩。**
HOMER 没配基因组、MEME 没装、R 包缺失，都只 WARN 跳过对应步骤。跑了四小时的
比对结果不该因为一个注释库没装就作废。只有 6 个必需工具缺失才直接退出。

**5. 为什么 TSS 富集自己算。**
ENCODE 的定义是 TSS ± 50 bp 平均深度 / 两侧 ±(1900-2000) bp 背景深度。
`samtools bedcov` 两次就能算出来，比引入 ATACseqQC 的 R 依赖更可控，
数值口径也和 ENCODE 一致。

## 常见问题

**chrM 比例 30% 以上怎么办。** 这是文库问题，不是分析问题。ATAC 的线粒体占比
反映细胞核完整性——样本冻融过、裂解过度、或起始细胞量太少都会这样。
分析上只能去掉 chrM 继续，但结论里必须说明核基因组信息有损失。

**peak 数很少（几千个）。** 依次查：比对率是否 > 80%、chrM 是否过高、
MACS2 的 `-g` 是否写成了 `hs`（写成 `2.7e9` 也行，但别写成 `1e9`）。

**`alignmentSieve` 之后 BAM 变小很多。** 正常，它只保留前 11 个字段（丢掉 tag）。
但**必须重新 sort**，否则 MACS2 会报 "chromosome blocks not continuous"。

**FRiP 两条算法差很多。** 最常见原因是 `plotEnrichment` 默认不做 MAPQ 过滤，
而 `bedtools coverage` 直接数 BAM 里的所有 read。本流程两边都用同一个 srt2.bam
（已经 `-q 20` 过），所以应该一致；差超过 5% 就检查是不是有人手工换过 BAM。

**HOMER 报找不到 hg38。** 需要先 `perl configureHomer.pl -install hg38`。
没装也不影响主流程，⑫⑬ 会 WARN 跳过，peak 文件仍然有效。

**想加 TOBIAS footprint。** 加 `--tobias`，但要把 `RUN_TOBIAS=1` 对应的
`src/` 逻辑补上（目前只留了配置位，没实现，因为 TOBIAS 的输入要求
未偏移的 srt.bam，和主流程的 srt2 不是同一份，容易搞混）。

## 测试

```bash
bash tests/smoke_test.sh      # 51 项，不需要任何生信工具
```

只验证脚手架：参数解析、路径派生、阶段跳过、命令拼装、dry-run 不落标记。
真实数据必须上服务器跑。
