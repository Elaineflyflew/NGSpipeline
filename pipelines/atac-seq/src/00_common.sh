#!/usr/bin/env bash
# =============================================================================
#  src/00_common.sh —— 公共基础设施
#
#  约定：本文件只定义函数，source 时不产生副作用（不建目录、不改 PATH）。
#  所有模块都依赖 log / run / need / stage_done，所以它必须最先加载。
#
#  与 germline-wes-wgs 的同名文件保持同一套接口，改一个流程时另一个不用重学。
# =============================================================================

VERSION="1.0.0"
PROG="atac-seq"

# ---------------------------------------------------------------- 日志
log()  { printf '[%s] %s\n'      "$(date '+%m-%d %H:%M:%S')" "$*"; }
warn() { printf '[%s] WARN  %s\n' "$(date '+%m-%d %H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] ERROR %s\n' "$(date '+%m-%d %H:%M:%S')" "$*" >&2; exit 1; }

section() { printf '\n[%s] ===== %s =====\n' "$(date '+%m-%d %H:%M:%S')" "$*"; }

# ---------------------------------------------------------------- 路径工具
first_existing() {
  local p
  for p in "$@"; do
    [[ -e "$p" ]] && { printf '%s' "$p"; return 0; }
  done
  printf ''
}

# 找一个目录下第一个匹配的文件（用于"放进去就自动用"的约定）
first_file() {
  local dir="$1" pat="$2" f
  [[ -d "$dir" ]] || return 0
  for f in "$dir"/$pat; do
    [[ -f "$f" ]] && { printf '%s' "$f"; return 0; }
  done
  printf ''
}

# ---------------------------------------------------------------- 物种配置
# 换物种时要把上一个物种探测出来的路径清掉，否则 ${VAR:=...} 不会覆盖旧值
GENOME_VARS="BT2_INDEX FASTA BLACKLIST GTF MOTIF_DB HOMER_MOTIFS CHROM_SIZES EFFECTIVE_GENOME"

load_genome() {
  local g="$1" conf="$PROJECT_ROOT/config/genomes/$1.conf" v
  [[ -s "$conf" ]] || die "没有这个物种的配置: ${conf}（可用的：$(ls "$PROJECT_ROOT/config/genomes" 2>/dev/null | sed 's/\.conf$//' | tr '\n' ' '))"
  for v in $GENOME_VARS; do unset "$v"; done
  GENOME="$g"
  # shellcheck source=/dev/null
  source "$conf"
  log "物种配置: $conf"
}

# ---------------------------------------------------------------- 工具定位
# 先看 PATH，再到项目自带目录和常见 conda/mamba 环境里找。
# 用户机器上往往已经装好了（可能在某个没激活的 env 里），逼他们重装既浪费
# 空间又可能换版本。找到就用，找不到才提示装。
find_tool() {
  local t="$1" p
  p=$(command -v "$t" 2>/dev/null) && { printf '%s' "$p"; return 0; }

  local home="${HOME:-/nonexistent}"
  for p in "$PROJECT_ROOT"/extratools/bin \
           "${CONDA_PREFIX:-/nonexistent}/bin" \
           /bioapps/bin \
           /home/biotools/miniforge3/envs/*/bin \
           /opt/conda/envs/*/bin \
           "$home"/miniforge3/envs/*/bin \
           "$home"/mambaforge/envs/*/bin \
           "$home"/miniconda3/envs/*/bin \
           "$home"/anaconda3/envs/*/bin; do
    [[ -x "$p/$t" ]] && { printf '%s' "$p/$t"; return 0; }
  done
  return 1
}

# ---------------------------------------------------------------- 前置检查
#   dry-run    只记录"需要什么"，不检查
#   --check    缺了就记一笔继续走，最后统一汇报
#   正常运行   缺了立刻退出
MISSING_TOOLS=0
need() {
  (( ${DRY_RUN:-0} )) && return 0

  local p
  if p=$(find_tool "$1"); then
    case ":$PATH:" in
      *":$(dirname "$p"):"*) ;;
      *) export PATH="$(dirname "$p"):$PATH" ;;
    esac
    return 0
  fi

  if (( ${CHECK_ONLY:-0} )); then
    warn "缺少命令: $1"
    MISSING_TOOLS=$(( MISSING_TOOLS + 1 ))
    return 0
  fi
  die "找不到命令 $1 —— 跑 'deploy/deploy.sh doctor' 看缺什么，再 'deploy/deploy.sh env' 补装"
}

# 工具存在但不报错：用于"有就多做一步，没有就跳过"的可选步骤
have() { (( ${DRY_RUN:-0} )) && return 0; find_tool "$1" >/dev/null 2>&1; }

require_arg() {
  [[ -n "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] 缺少 %s（实际运行会报错）\n' "$2"; return 0; }
  die "$2 未指定：$3"
}

checkfile() {
  [[ -s "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] 文件不存在（实际运行会报错）: %s\n' "$1"; return 0; }
  die "$2 不存在: $1 —— 检查 config/genomes/${GENOME}.conf"
}

checkdir() {
  [[ -d "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] 目录不存在（实际运行会报错）: %s\n' "$1"; return 0; }
  die "$2 不存在: $1"
}

# ---------------------------------------------------------------- 命令执行
show() {
  local a
  for a in "$@"; do
    if [[ "$a" == *" "* || -z "$a" ]]; then printf "'%s' " "$a"; else printf '%s ' "$a"; fi
  done
  printf '\n'
}

run() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] '; show "$@"; return 0; fi
  printf '+ '; show "$@"
  "$@"
}

# 管道用：$1 是完整命令字符串（仅用于本脚本内部拼好的命令）
runpipe() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] %s\n' "$1"; return 0; fi
  printf '+ %s\n' "$1"
  eval "$1"
}

# ---------------------------------------------------------------- 阶段状态
stage_enabled() { [[ ",${STAGES:-}," == *",$1,"* ]]; }
stage_done()    { [[ -f "$STATUS/$1.done" ]] && (( ! ${FORCE:-0} )); }
mark_done()     { (( ${DRY_RUN:-0} )) && return 0; touch "$STATUS/$1.done"; }

run_stage() {
  local name="$1" fn="$2"
  if stage_done "$name"; then
    log "阶段 [$name] 已完成，跳过（--force 可强制重跑）"
    return 0
  fi
  "$fn"
}

# ---------------------------------------------------------------- 文本排版
# 按显示宽度补齐：中文占 2 列，否则表格会歪
disp_width() {
  LC_ALL=C awk -v s="$1" 'BEGIN{
    n = 0
    for (i = 1; i <= length(s); i++) {
      c = substr(s, i, 1)
      n += (c ~ /[\x80-\xff]/) ? 2 : 1
    }
    print n
  }'
}

pad() {
  local s="$1" want="$2" w fill
  w=$(disp_width "$s")
  fill=$(( want - w ))
  (( fill < 0 )) && fill=0
  printf '%s%*s' "$s" "$fill" ""
}

# 一行两列表格
row() {
  local k="$1" v="$2" w="${3:-14}"
  printf '  %s %s\n' "$(pad "$k" "$w")" "$v"
}

# 人类可读大小，跟随软链接（du -h 默认不跟，参考基因组会显示 0）
human_size() {
  local p="$1" s
  [[ -e "$p" ]] || { printf '-'; return 0; }
  s=$(du -Lsh "$p" 2>/dev/null | awk '{print $1}')
  [[ -n "$s" ]] || s=$(ls -ld "$p" 2>/dev/null | awk '{print $5}')
  printf '%s' "${s:--}"
}

# 百分比：两个整数 → 保留一位小数
pct() {
  awk -v a="$1" -v b="$2" 'BEGIN{ printf "%.1f", (b > 0 ? 100 * a / b : 0) }'
}
