# extratools/ —— 非 conda 安装的外部工具

有些工具不适合放 conda（体积大、依赖冲突、需要许可证），统一放这里。
`bin/wes-pipeline` 启动时会把 `extratools/bin` 加进 PATH，
所以自己编译的小工具丢进 `bin/` 就能直接调用。

```
extratools/
├── deepvariant/
│   └── deepvariant.sif          DeepVariant 容器镜像（约 5 G）
├── annovar/
│   ├── annotate_variation.pl    ANNOVAR 主程序（需学术注册下载）
│   ├── table_annovar.pl
│   └── humandb/                 各数据库
└── bin/                         其他可执行文件，自动加进 PATH
```

## 为什么 DeepVariant 用容器

DeepVariant 依赖特定版本的 TensorFlow/protobuf，和 `gatk4` 装同一个 conda
环境几乎必然冲突。官方推荐容器运行，所以这里放 `.sif`。

```bash
bash deploy/deploy.sh deepvariant          # 默认 1.6.1
bash deploy/deploy.sh deepvariant --version 1.5.0
```

镜像版本要**固定**，不要用 latest：换版本会静默改变 call 结果，报告没法复现。

## 为什么 ANNOVAR 是半手动

ANNOVAR 需要学术注册，脚本不能替你点同意。流程：

1. <https://www.openbioinformatics.org/annovar/annovar_download_form.php> 填表
2. 邮件收到下载链接，下载 `annovar.latest.tar.gz` 放到 `extratools/`
3. `bash deploy/deploy.sh annovar --tarball extratools/annovar.latest.tar.gz`

没有 ANNOVAR 不影响其他步骤，`-a vep` 或 `-a funcotator` 照样能用。

## 容器引擎

`deploy/50_deepvariant.sh` 优先用 `singularity`，其次 `apptainer`，都没有则打印指引。
conda 环境里已经带了 apptainer：

```bash
conda activate wes
apptainer --version
```

如果集群有 module 系统，`module load singularity` 后再跑 deploy 也可以。

## 磁盘

约 25 G（DeepVariant 5 G + ANNOVAR 10 G+ + 余量）。不进 git，见 `.gitignore`。
