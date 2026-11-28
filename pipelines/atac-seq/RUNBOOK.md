# 服务器实操手册

> 目标：从登录服务器到拿到 `report.html` 和 IGV 截图，中途断了能续。
> 主入口只有一条：`bin/atac`。配置只改一个文件：`config/atac.conf`。

---

## 0. 上传与体检（5 分钟，省几小时）

```bash
# 本地（results / reference / extratools 都不传，服务器上重新生成）
rsync -avz --exclude results --exclude reference --exclude extratools \
      --exclude '.git' ./ 你的账号@服务器:~/lab9_atac/

# 服务器
cd ~/lab9_atac
df -h .                 # 全程约需 120 GB（含 20 GB SAM）
nproc                   # 拿到真实核数，传给 -t
free -g                 # bowtie2 + sort 峰值约 8–12 GB
```

**先体检再开跑**：

```bash
./bin/atac -1 ... -2 ... --check
```

它逐项打印工具、输入、参考数据、输出空间，缺什么直接指出来，**不碰任何数据**。
不想先给输入文件也行：`./bin/atac --check` 会跳过输入校验，只查工具链。

---

## 1. 环境

**情况 A：能联网、有 conda**

```bash
conda env create -f env/environment.yml
conda activate atac
./bin/atac --check
```

**情况 B：课程镜像已经装好（`/bioapps/bin` 之类）**

不用建环境，只要保证命令在 `PATH` 里。`config/atac.conf` 的 `EXTRA_PATH`
已经把课程服务器的 `/bioapps/bin` 和 `/home/biotools/miniforge3/envs/motif/bin`
写进去了，目录不存在会自动忽略。手工确认：

```bash
for t in fastqc fastp bowtie2 samtools bedtools macs2 alignmentSieve \
         bamCoverage computeMatrix plotProfile plotFingerprint \
         annotatePeaks.pl findMotifsGenome.pl Rscript; do
  printf "%-22s %s\n" "$t" "$(command -v $t || echo MISSING)"
done
```

`annotatePeaks.pl` / `findMotifsGenome.pl` 可能装在单独的 env 里，
把那个 `bin` 加进 `config/atac.conf` 的 `EXTRA_PATH` 即可；`--check` 会把它们标成「可选」。

**情况 C：什么都没装，从零开始**

```bash
./deploy/deploy.sh doctor      # 看缺什么
./deploy/deploy.sh all         # 按依赖顺序全部装（约 7 GB）
```

`deploy.sh` 只负责「把东西放到该放的位置」，不修改你的 shell 配置。
详见 `deploy/deploy.sh -h` 和 `reference/README.md`。

---

## 2. 参考数据（只做一次）

```bash
./deploy/deploy.sh status      # 哪些已经有了
./deploy/deploy.sh all         # 缺什么补什么，已有的自动跳过
```

装完的布局（`config/atac.conf` 默认就指这里，不用改路径）：

```
reference/
├── genome/hg38.fa (+.fai)          deploy.sh genome       ~1.0 G
├── genome/bt2/hg38.*.bt2           deploy.sh bt2          ~4.0 G
├── annotation/gencode.v43.gtf.gz   deploy.sh annotation   ~50 M
└── blacklist/hg38-blacklist.v2.bed deploy.sh blacklist    ~10 K
extratools/
├── bin/{annotatePeaks.pl,findMotifsGenome.pl,igv}
└── homer/                          deploy.sh homer        ~2.0 G
```

**服务器上已经有参考数据**（你说了不用下）→ 三条路任选：

```bash
# 1) 软链接进 reference/，之后完全不用改配置
ln -s /public/ref/hg38.fa        reference/genome/hg38.fa
ln -s /public/ref/bt2/hg38       reference/genome/bt2/hg38

# 2) 命令行直接给
./bin/atac ... -g /public/ref/hg38.fa --bt2 /public/ref/bt2/hg38

# 3) 写进环境变量（优先于配置文件）
export GENOME=/public/ref/hg38.fa BT2=/public/ref/bt2/hg38
```

留空的项会**自动探测**：`reference/genome/` 下按 `.fai` 和文件名打分挑一个 fasta，
`annotation/` 和 `blacklist/` 各取第一个匹配的文件。所以只要把文件丢进对应目录就行。

---

## 3. 开跑

```bash
tmux new -s atac
conda activate atac
cd ~/lab9_atac

./bin/atac \
  -1 /public/storage/shaojf/Course/NGS/wkdir4labs/lab9/rawdata/sub20M.SRR13849313_1.fastq.gz \
  -2 /public/storage/shaojf/Course/NGS/wkdir4labs/lab9/rawdata/sub20M.SRR13849313_2.fastq.gz \
  -t 20 2>&1 | tee logs/run_$(date +%m%d_%H%M).log
```

`Ctrl-b d` 挂起，`tmux attach -t atac` 回来看。
不给 `-o` 时结果落在 `./results/`，日志同时写一份到 `results/00_logs/`。

先看会跑什么（不执行）：

```bash
./bin/atac -1 ... -2 ... --dry-run
./bin/atac --list-stages        # 13 个阶段各自干什么、产出什么
```

---

## 4. 时间预期（20 线程，20 M pairs）

| 步骤 | 阶段名 | 大致耗时 | 说明 |
|---|---|---|---|
| ① FastQC | `qc` | 5–10 min | 两文件各扫一遍 |
| ② fastp | `trim` | 3–5 min | |
| ③ bowtie2 | `align` | 3–5 min | 200 线程时实测 2 m 33 s |
| ④ 过滤 + 排序 | `bam` | 5–8 min | 受磁盘 IO 影响最大 |
| ⑤ alignmentSieve | `bam` | 20–40 s | 实测 18 s |
| ⑥ 重排 | `bam` | 3–5 min | |
| ⑦ 分箱 | `frag` | 5–10 min | 单次遍历写三份 |
| ⑧ MACS2 ×3 | `peaks` | 5–10 min | |
| ⑨ 黑名单 | `peaks` | < 1 min | |
| ⑩ FRiP | `frip` | 5–15 min | bedtools 那一路慢 |
| ⑪ bigWig + TSS | `coverage` | 10–20 min | computeMatrix 是瓶颈 |
| ⑫ 注释 ×3 | `annot` | 10–15 min | |
| ⑬ motif ×3 | `motif` | 30–60 min | 已取 top 1 万；全量要 2 h/路线 |
| ⑭ ATACseqQC | `atacqc` | 10–20 min | 首次要加载 TxDb |
| ⑮ IGV | `igv` | 1 min | 无图形环境时自动降级 |
| 报告 | `report` | < 1 min | |
| **合计** | | **约 2–3 小时** | 主要变量是 ⑬ |

---

## 5. 断点续跑

每个完成的阶段在 `results/.status/` 留一个 `.done` 标记，重跑时自动跳过。

```bash
./bin/atac -1 ... -2 ...                      # 自动跳过已完成，从断点继续
./bin/atac -1 ... -2 ... --stages peaks       # 只跑这一段
./bin/atac -1 ... -2 ... --stages peaks,report --force   # 强制重跑这两段
./bin/atac -1 ... -2 ... --force              # 全部重来
```

失败时会打印 `[FATAL] 函数名() 第 N 行: 命令`（`ERR` trap 抓的，含出错的具体命令），
并提示该从哪个阶段继续。日志在 `results/00_logs/`。

阶段之间有依赖，`--stages` 选的是**从哪开始**而不是「只跑这个」——
比如 `--stages peaks` 会从 ⑧ 一路跑到报告，因为 ⑩⑪⑫ 都读 ⑨ 的 peak。
真要只跑一段，把后面的阶段名也去掉：`--stages qc,trim`。

---

## 6. 跑完要拉回本地的文件

```bash
# 报告和指标（体积小，必拉）
results/99_report/
results/metrics/
results/06_peaks/*.filtered.narrowPeak
results/11_motif/*/knownResults.txt
results/11_motif/*/top_motifs.tsv
results/10_annot/*/annotation_distribution.tsv
results/05_fragbin/*.fragsize_dist.png
results/13_igv/*.png

# 大文件按需（bigWig 可拉到本地 IGV 看）
results/08_bigwig/*.bw
```

打包：

```bash
tar czf lab9_results_$(date +%m%d).tar.gz \
    results/99_report results/metrics results/06_peaks/*.filtered.narrowPeak \
    results/11_motif results/10_annot results/05_fragbin/*.png results/13_igv
```

---

## 7. 报错对照表

| 报错 | 原因 | 处理 |
|---|---|---|
| `Too many open files` | 线程开太多，文件句柄耗尽 | `ulimit -n 65536` 或把 `-t` 降到 8 |
| `samtools index: failed to create index` | 上一步没排好序 / 被 SIGPIPE 截断 | `--stages bam --force` 重跑 |
| `Could not find genome hg38`（HOMER） | 没装 hg38 | `./deploy/deploy.sh homer` |
| `alignmentSieve: command not found` | deeptools 没装 | `conda install -c bioconda deeptools` |
| `Killed`（OOM） | sort 内存超了 | 调小 `SORT_MEM`（如 1G），或降 `-t` |
| `IGV: no display` | 无图形环境 | 自动降级成 `xvfb-run` / 只留 batch 文件 |
| `No space left on device` | SAM 太大 | `KEEP_SAM=0` 重跑 `--stages align,bam --force` |
| `ATACseqQC 步骤 SKIP` | R 包缺失 | `./deploy/deploy.sh atacqc`，不影响主链 |
| `未知阶段: xxx` | `--stages` 拼错 | `./bin/atac --list-stages` 看合法值 |
| `[dry-run] ...` 后停下 | `--dry-run` 只打印不执行 | 去掉 `--dry-run` 重跑 |

---

## 8. 自检清单

跑完对照这几条，全过才算这条链路走通：

- [ ] ③ 比对率 > 95%
- [ ] ③ chrM 占比 < 20%（理想 < 10%）
- [ ] ⑦ 片段分布有三个峰（<100 / ~200 / ~400 bp）
- [ ] ⑦ NFR 占比 30–60%
- [ ] ⑨ 黑名单剔除 < 1%
- [ ] ⑩ FRiP > 0.20，两算法差值 < 0.02
- [ ] ⑪ TSS 富集 > 6（ENCODE），曲线呈 V 形
- [ ] ⑫ promoter 占比 20–35%
- [ ] ⑬ top motif 含 CTCF / NFY / SP1 / ETS 之一
- [ ] ⑮ MYC 区有清晰 peak
- [ ] `99_report/report.html` 能打开，无 SKIP 项（除 ⑭ 外）

`report.html` 里的「三指标合看」段落会自动给这三项下结论
（片段周期 / FRiP / TSS 富集），不用手工比对。
