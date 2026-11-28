#!/usr/bin/env bash
################################################################################
# run_workflow.sh —— Snakemake 入口的包装脚本
#
# 存在的理由只有一个：snakemake 那一串参数（-s / -j / --configfile /
# --rerun-incomplete / --use-conda）每次手敲容易漏，漏一个就是几小时的返工。
# 这里把默认值钉死，剩下的原样透传给 snakemake。
#
# 用法：
#   bash run_workflow.sh                     # 用所有核开跑
#   bash run_workflow.sh -j 32               # 指定核数
#   bash run_workflow.sh -n                  # 只构建 DAG，不执行
#   bash run_workflow.sh --use-conda         # 每个 rule 用 envs/*.yaml 的独立环境
#   bash run_workflow.sh --unlock            # 上次异常中断后解锁
#   bash run_workflow.sh -- dag | dot -Tpdf > dag.pdf
#
# 其余参数（--forcerun / --until / --dry-run / --profile …）直接透传。
################################################################################

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

SNAKEFILE="workflow/Snakefile"
CONFIG="config/config.yaml"

# ---------------------------------------------------------------------------
# 默认核数：优先 nproc，macOS 没有 nproc 就用 sysctl，再不行给 8
# ---------------------------------------------------------------------------
detect_cores() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    elif command -v sysctl >/dev/null 2>&1; then
        sysctl -n hw.ncpu 2>/dev/null || echo 8
    else
        echo 8
    fi
}
CORES="${CORES:-$(detect_cores)}"

die() { printf '\033[31m[ERROR]\033[0m %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 预检：没有 snakemake 就直接说清楚怎么装，别让人对着 "command not found" 猜
# ---------------------------------------------------------------------------
command -v snakemake >/dev/null 2>&1 || die \
    "找不到 snakemake。装一个：
    mamba install -c bioconda -c conda-forge snakemake-minimal>=7.32
  或者不装 snakemake，直接用 bash 主程序：./bin/atac -1 R1 -2 R2"

[[ -f "$SNAKEFILE" ]] || die "找不到 $SNAKEFILE，请在仓库根目录运行本脚本"
[[ -f "$CONFIG"   ]] || die "找不到 $CONFIG"

# 版本下限：wildcard_constraints / resources / 新版 report 行为都要求 >= 7.32
SMK_VER="$(snakemake --version 2>/dev/null | head -1)"
python3 - "$SMK_VER" <<'PY' || die "snakemake 版本过低（$SMK_VER），需要 >= 7.32"
import sys
v = sys.argv[1].strip().split(".")
try:
    major, minor = int(v[0]), int(v[1])
except (IndexError, ValueError):
    sys.exit(1)
sys.exit(0 if (major, minor) >= (7, 32) else 1)
PY

# ---------------------------------------------------------------------------
# 参数分流：-j / --cores 由本脚本接管，其余全部透传
# ---------------------------------------------------------------------------
EXTRA=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--cores)   CORES="${2:?-j 后面要跟数字}"; shift 2 ;;
        -j*)          CORES="${1#-j}"; shift ;;
        --cores=*)    CORES="${1#--cores=}"; shift ;;
        -h|--help)
            sed -n '3,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *)            EXTRA+=("$1"); shift ;;
    esac
done

printf '\033[36m[INFO]\033[0m snakemake %s  |  核数 %s  |  %s\n' "$SMK_VER" "$CORES" "$SNAKEFILE"

exec snakemake \
    -s "$SNAKEFILE" \
    --configfile "$CONFIG" \
    -j "$CORES" \
    --rerun-incomplete \
    --printshellcmds \
    "${EXTRA[@]}"
