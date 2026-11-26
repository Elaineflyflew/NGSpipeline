# reference/ —— 参考数据

这个目录放**只读的参考数据**，全部由 `deploy/` 下载，不要手工往里塞文件
（除非你明确知道自己在做什么）。

```
reference/
├── genome/           参考基因组 fasta + 索引（文件名随意，见下）
│   ├── Homo_sapiens_assembly38.fasta        broad 版（chr 前缀）
│   ├── Homo_sapiens_assembly38.fasta.fai
│   ├── Homo_sapiens_assembly38.dict
│   └── Homo_sapiens_assembly38.fasta.{amb,ann,bwt,pac,sa,alt}
├── gatk_bundle/      GATK 已知位点（.vcf 或 .vcf.gz 都可以，见下）
│   ├── dbsnp_146.hg38.vcf        或 .vcf.gz + .tbi
│   ├── hapmap_3.3.hg38.vcf
│   ├── Mills_and_1000G_gold_standard.indels.hg38.vcf
│   └── 1000G_phase1.snps.high_confidence.hg38.vcf
├── vep/
│   └── homo_sapiens/115_GRCh38/             VEP cache（--dir_cache 指到 vep/）
├── funcotator/
│   └── funcotator_dataSources.v1.8.hg38.../ Funcotator 数据源
└── bed/
    └── *.bed                                WES 捕获区间，放进来就自动用
```

## 怎么准备

```bash
bash deploy/deploy.sh status      # 先看缺什么
bash deploy/deploy.sh ref bundle  # 下载
```

## 参考基因组怎么被找到

不需要软链接，也不需要固定文件名。pipeline 扫描 `genome/` 下所有
`*.fa` / `*.fasta` / `*.fna`，按下面的规则打分，取最高分的那个：

| 加分 | 条件 |
|---|---|
| +2 | 有 `.fai` |
| +2 | 有 `.bwt`（bwa 索引） |
| +1 | 有 `.dict` |
| +1 | 名字里带 `assembly38` / `GRCh38` / `hg38` |

所以下面这些都能直接用，不用改任何配置：

```
genome/Homo_sapiens_assembly38.fasta        ← deploy 下载的
genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa   ← 课程服务器那份
genome/ref.fa                                ← 自己重命名的
```

想确认它选中了哪个，跑 `bin/germline-wes-wgs -h` 之外最简单的办法是看启动横幅：

```
=== germline-wes-wgs 2.0.0 | 样本 Case1 ===
参考: /path/to/reference/genome/Homo_sapiens_assembly38.fasta
```

要指定别的文件：`-r /other/ref.fa` 或 `export REF=/other/ref.fa`。

## gatk_bundle 的两种格式

GATK 两种都吃，脚本也都认：

| 格式 | 索引 | 说明 |
|---|---|---|
| `xxx.vcf.gz` | `xxx.vcf.gz.tbi` | `deploy` 下载的官方形态 |
| `xxx.vcf` | `xxx.vcf.idx` | 解压后的形态，GATK 自己建的那种索引 |

判断规则：同名 `.vcf` 或 `.vcf.gz` 任一存在就算就绪，两种都在时优先 `.vcf.gz`。
所以下面这种目录直接能用，`deploy/deploy.sh bundle` 不会重复下载：

```
dbsnp_146.hg38.vcf            dbsnp_146.hg38.vcf.idx
hapmap_3.3.hg38.vcf           hapmap_3.3.hg38.vcf.idx
Mills_and_1000G_gold_standard.indels.hg38.vcf   (+ .idx)
1000G_phase1.snps.high_confidence.hg38.vcf      (+ .idx)
```

**未压缩 VCF 缺 `.idx` 时会自动补建**：预检会提示，真正跑 GATK 时先执行
`gatk IndexFeatureFile -I <file>` 再往下走。`.idx` 为 0 字节也算缺失（半截索引不能用）。

`1000G_omni2.5.hg38.vcf` 之类多出来的文件不影响，pipeline 只按名字取需要的四个。

## 为什么用 broad 版参考

GATK 官方的 `Homo_sapiens_assembly38.fasta`：

- contig 名带 `chr` 前缀，和 `dbsnp_146.hg38`、`hapmap_3.3.hg38` **天然一致**
- 自带 `.fai` / `.dict` / bwa 索引，省掉 40 分钟 `bwa index`
- 带 `.alt` 文件，GATK 做 BQSR 时能正确处理 alt contig

用 Ensembl 版（无 `chr` 前缀）会踩命名不一致的坑，详见 README 常见问题。

## bed/ 的用法

`config/pipeline.conf` 会自动取 `bed/` 下第一个 `.bed` 作为捕获区间。
所以放了 WES 的 bed 之后，命令行不用再加 `-b`。

多个 bed 想指定其中一个，还是用 `-b` 显式传。

## 磁盘

全部约 50 G（不含 ANNOVAR）。这个目录不应该进 git，见 `.gitignore`。
