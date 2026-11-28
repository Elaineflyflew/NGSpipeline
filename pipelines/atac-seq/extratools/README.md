# extratools/ —— 非 conda 安装的外部工具

放那些 conda 装不了、或者你想钉死版本的工具。目录里的可执行文件会被自动加进
`PATH`（`20_context.sh` 初始化时做的），不需要手工 `export`。

```
extratools/
└── bin/            自己编译或下载的单文件工具
```

平时是空的。这个流程用到的 HOMER、IGV、ATACseqQC 通常已经装在系统或 conda 里，
`find_tool` 会自动扫下面这些位置，找到就直接用，找不到才由 `deploy/deploy.sh` 补装：

- `PATH` 本身
- `extratools/bin`（就是这里）
- `$CONDA_PREFIX/bin`
- `/bioapps/bin`（课程镜像的公共工具目录）
- `/home/biotools/miniforge3/envs/*/bin`
- `~/miniforge3|mambaforge|miniconda3|anaconda3/envs/*/bin`

想加自己的搜索目录，设 `TOOL_SEARCH_DIRS`（空格分隔），会插在 `extratools/bin`
之后、系统目录之前。

## 什么时候才需要往这里放东西

- 要用某个特定版本的 `macs2`，但不想动 conda 环境
- 拿到的 IGV 是解压即用的发行包（`igv.sh` 旁边那一堆 jar），不想装到系统里
- 自己改过的 `findMotifsGenome.pl`

放进来之后 `bash deploy/deploy.sh doctor` 就能看到它被识别到了。
