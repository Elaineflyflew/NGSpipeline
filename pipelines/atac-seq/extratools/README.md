# extratools/ —— 非 conda 安装的外部工具

放那些 conda 装不了、或者你想用特定版本的工具。目录里的可执行文件会被
自动加进 PATH（`20_context.sh` 里做的），不需要手工 export。

```
extratools/
└── bin/            自己编译或下载的单文件工具
```

HOMER、MEME、TOBIAS 这些通常装在独立的 conda 环境里，不放在这里 ——
`find_tool` 会自动扫 `/home/biotools/miniforge3/envs/*/bin` 和
`$HOME/{miniforge3,mambaforge,miniconda3,anaconda3}/envs/*/bin`，
找到就把路径临时加进 PATH，所以不 activate 也能用。
