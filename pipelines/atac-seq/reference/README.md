# reference/ —— 参考数据

**这个流程默认不下载任何数据。** 参考基因组、索引、blacklist、GTF、motif 库
都在 `config/genomes/<name>.conf` 里按候选路径探测，服务器上已有的直接复用。

所以这个目录平时是空的。只有在两种情况下才往里放东西：

1. 服务器上的数据不在 `config/genomes/hg38.conf` 列出的位置，
   你又不想改配置文件 —— 把文件放进来，约定名即可被探测到。
2. 想给这个流程做一份自包含的数据快照（比如拷到另一台机器）。

约定的目录结构：

```
reference/
├── genome/
│   ├── hg38.fa                     参考 fasta
│   ├── hg38.fa.fai                 samtools faidx 生成
│   └── hg38.{1..4}.bt2  hg38.rev.1.bt2  hg38.rev.2.bt2   bowtie2-build 生成
├── blacklist/
│   └── hg38-blacklist.v2.bed       ENCODE blacklist
├── annotation/
│   └── gencode.v43.annotation.gtf  GENCODE 注释
└── motif/
    ├── JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
    └── known.motifs                HOMER knownTFs
```

换物种时不要往这里塞第二个物种的文件，而是新建
`config/genomes/<name>.conf` 指向对应数据 —— 这样同一份代码能同时支持多个物种。
