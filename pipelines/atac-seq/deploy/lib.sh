#!/usr/bin/env bash
# =============================================================================
#  deploy/lib.sh —— deploy 公共函数
#
#  日志 / die / need / first_existing 直接复用 src/00_common.sh，不在这里
#  重复实现——同一个函数两处维护迟早会不一致。
#  本文件只放 deploy 特有的东西：路径常量、下载、空间检查、表格输出。
# =============================================================================

: "${PROJECT_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PROJECT_ROOT

# shellcheck source=/dev/null
source "$PROJECT_ROOT/src/00_common.sh"

# 路径全部用 := 默认值，允许用环境变量改到别处（多套数据共存时用得上）
: "${REF_DIR:=$PROJECT_ROOT/reference}"
: "${GENOME_DIR:=$REF_DIR/genome}"
: "${BT2_DIR:=$GENOME_DIR/bt2}"
: "${ANN_DIR:=$REF_DIR/annotation}"
: "${BL_DIR:=$REF_DIR/blacklist}"
: "${TOOLS_DIR:=$PROJECT_ROOT/extratools}"
: "${HOMER_DIR:=$TOOLS_DIR/homer}"
: "${IGV_DIR:=$TOOLS_DIR/igv}"
: "${BIN_DIR:=$TOOLS_DIR/bin}"

: "${DRY_RUN:=0}"
: "${FORCE:=0}"

# ---------------------------------------------------------------- 下载
# dl URL OUT —— 下载单个文件，已存在则跳过
dl() {
  local url="$1" out="$2"
  mkdir -p "$(dirname "$out")"

  if [[ -s "$out" ]] && (( ! FORCE )); then
    log "已存在，跳过: $(basename "$out")"
    return 0
  fi
  if (( DRY_RUN )); then
    printf '[dry-run] 下载 %s\n           -> %s\n' "$url" "$out"
    return 0
  fi

  log "下载 $(basename "$out")"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 8 -s 8 -c --auto-file-renaming=false --allow-overwrite=true \
           --console-log-level=warn --summary-interval=0 \
           -d "$(dirname "$out")" -o "$(basename "$out")" "$url" \
      || die "下载失败: $url"
  else
    curl -fL --retry 3 --retry-delay 5 -C - -o "$out" "$url" || die "下载失败: $url"
  fi
  [[ -s "$out" ]] || die "下载后文件为空: $out"
}

# 简单的空间检查：下载前提醒，避免下到一半 No space left
check_space() {
  local need_gb="$1" dir="$2" avail
  mkdir -p "$dir" 2>/dev/null || true
  avail=$(df -g "$dir" 2>/dev/null | awk 'NR==2{print $4}')
  [[ -z "$avail" ]] && return 0
  if (( avail < need_gb )); then
    warn "$dir 剩余 ${avail}G，本步骤大约需要 ${need_gb}G"
  fi
}

# ---------------------------------------------------------------- 表格输出
# 按"显示宽度"补齐：中文占 2 列，printf 的 %-22s 按字节补会错位
_dwidth() {
  local s="$1" b c
  b=$(printf '%s' "$s" | wc -c | tr -d ' ')
  c=${#s}
  printf '%s' $(( c + (b - c) / 2 ))
}
pad() {
  local s="$1" w="$2" n
  n=$(_dwidth "$s")
  printf '%s%*s' "$s" $(( w > n ? w - n : 0 )) ''
}
# row 标签 状态 说明
row() { printf '%s %s %s\n' "$(pad "$1" 20)" "$(pad "$2" 12)" "$3"; }

# 人类可读大小。du 默认不跟随软链接，遇到 ln -s 过来的参考会显示 0，所以用 du -L
human_size() {
  local p="$1" s=""
  [[ -e "$p" ]] || { printf '-'; return 0; }
  s=$(du -Lsh "$p" 2>/dev/null | awk '{print $1}' || true)
  [[ -n "$s" ]] || s=$(ls -ldh "$p" 2>/dev/null | awk '{print $5}')
  printf '%s' "${s:-?}"
}

# 项目里现在用的是哪个基因组 / 哪个 bt2 前缀（和 src/20_context.sh 同一套探测）
find_genome() { discover_genome "$GENOME_DIR"; }
find_bt2()    { discover_bt2 "$BT2_DIR"; }
