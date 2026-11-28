# `_shared/` —— 所有流程共用的基础设施

一个文件：`common.sh`。各流程的 `src/00_common.sh` 在**最前面** source 它：

```bash
_SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../_shared" && pwd)"
source "$_SHARED_DIR/common.sh"
```

路径是按 `BASH_SOURCE` 算的，所以从任何目录调用入口或 deploy 脚本都能找到，
不依赖当前工作目录。

## 收录判据

只有一条：**把这个函数原样搬到另一个流程，一行都不用改。**

放得进来：

| 类别 | 函数 |
|---|---|
| 日志 | `log` `ok` `warn` `die` `section` |
| 路径 | `first_existing` |
| 工具定位 | `find_tool`（PATH → 项目 `extratools/bin` → `/bioapps/bin` → 常见 conda 环境） |
| 前置检查 | `need` `need_file` `need_var` `require_arg` `checkdir` + `MISSING_TOOLS` / `SOFT_MISSING` |
| 命令执行 | `show` `run` `runpipe` `fail_cmd` `try_run` `try_runpipe` |
| 断点状态 | `stage_enabled` `stage_done` `mark_done` `run_stage` |
| 工具兼容 | `fastqc_probe` `pick_python` `pick_rscript` |
| 小工具 | `slug` `stats_field` `fastqc_stem` `fastqc_field` |

放不进来：

- `discover_*` —— 每个流程探测的参考数据结构不同（bwa 的 `.bwt/.dict` 与
  bowtie2 的 `.1.bt2`…`rev.2.bt2` 是两套完全不同的约定）。
- `stage_*` —— 业务逻辑，一个流程一套。
- `bundle_file()` —— GATK bundle 的命名规则，只有变异检测流程关心。
- `metric()` —— 五列指标格式是 atac-seq 定的，别的流程未必照抄。
- `contig_expr()` —— 依赖 `DROP_CONTIGS` 这个 atac-seq 才有的配置项。

## 改这里的风险

`common.sh` 是两边共用的，改错一次两条流程一起坏。所以：

```bash
# 改完必须两边都跑
cd pipelines/germline-wes-wgs && bash tests/smoke_test.sh   # 57 项
cd pipelines/atac-seq        && bash tests/smoke_test.sh    # 51 项
                             && bash tests/run_tests.sh     # 15 项
```

两套测试都不需要生信工具，几秒就能跑完，覆盖 `--dry-run` 的命令拼装和断点跳过，
正好是共用层改动最容易打坏的地方。

## 已知的环境假设

`find_tool` 的搜索列表里写死了一个课程镜像的路径：

- `/home/biotools/miniforge3/envs/*/bin` —— 预装的一批 conda 环境

别的机器上没有这个目录，`[[ -x ... ]]` 判假，不会有任何副作用。

**`/bioapps/bin` 故意没有写在这里。** 它属于"这个流程信任哪些目录"的策略，
放在各流程 `config` 的 `EXTRA_PATH` 里（`20_context.sh` 初始化时加进 `PATH`）。
写进共用层会和 config 重复，而且会让测试没法用桩目录模拟"某个工具没装"——
删掉桩里的 gatk，真的 gatk 还是会从 `/bioapps/bin` 被扫出来。

想加自己的搜索目录，设环境变量 `TOOL_SEARCH_DIRS`（空格分隔），会插在项目内
`extratools/bin` 之后、系统目录之前。
