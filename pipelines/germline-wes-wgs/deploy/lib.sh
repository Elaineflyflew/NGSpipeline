#!/usr/bin/env bash
# =============================================================================
#  deploy/lib.sh —— deploy 公共函数
#
#  日志 / die / need / first_existing / discover_ref 直接复用 src/00_common.sh，
#  不在这里重复实现——同一个函数两处维护迟早会不一致。
#  本文件只放 deploy 特有的东西：路径常量、下载、空间检查。
# =============================================================================

: "${PROJECT_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PROJECT_ROOT

# shellcheck source=/dev/null
source "$PROJECT_ROOT/src/00_common.sh"

# 路径全部用 := 默认值，允许用环境变量改到别处（测试和多套数据共存都用得上）
: "${REF_DIR:=$PROJECT_ROOT/reference}"
: "${GENOME_DIR:=$REF_DIR/genome}"
: "${BUNDLE_DIR:=$REF_DIR/gatk_bundle}"
: "${VEP_DIR:=$REF_DIR/vep}"
: "${FUNC_DIR:=$REF_DIR/funcotator}"
: "${BED_DIR:=$REF_DIR/bed}"
: "${TOOLS_DIR:=$PROJECT_ROOT/extratools}"
: "${DV_DIR:=$TOOLS_DIR/deepvariant}"
: "${ANNOVAR_DIR:=$TOOLS_DIR/annovar}"

# 让 src/00_common.sh 的 bundle_file() 在 deploy 里也能用（它读 ${BUNDLE}）
: "${BUNDLE:=$BUNDLE_DIR}"

# 判断某个已知位点是否已就绪（两种格式任一存在即可）
have_bundle_vcf() {
  local n="$1"
  [[ -s "$BUNDLE_DIR/$n.vcf" || -s "$BUNDLE_DIR/$n.vcf.gz" ]]
}

: "${DRY_RUN:=0}"
: "${FORCE:=0}"

# gs://bucket/path → https://storage.googleapis.com/bucket/path
gs_to_https() {
  case "$1" in
    gs://*) printf 'https://storage.googleapis.com/%s' "${1#gs://}" ;;
    *)      printf '%s' "$1" ;;
  esac
}

# dl URL OUT —— 下载单个文件，已存在则跳过
dl() {
  local url out
  url=$(gs_to_https "$1"); out="$2"
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

# dl_tar URL DIR —— 把压缩包下到 DIR 下并解压
dl_tar() {
  local url="$1" dest="$2" tarball
  tarball="$dest/$(basename "${url%%\?*}")"
  dl "$url" "$tarball"
  if (( DRY_RUN )); then
    printf '[dry-run] 解压 %s -> %s\n' "$tarball" "$dest"
    return 0
  fi
  log "解压 $(basename "$tarball")"
  case "$tarball" in
    *.tar.gz|*.tgz) tar -xzf "$tarball" -C "$dest" ;;
    *.tar.bz2)      tar -xjf "$tarball" -C "$dest" ;;
    *.zip)          unzip -q -o "$tarball" -d "$dest" ;;
    *)              die "不认识的压缩格式: $tarball" ;;
  esac
}

# 简单的空间检查：下载前提醒，避免下到一半 No space left
check_space() {
  local need_gb="$1" dir="$2" avail
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
# 状态列也按显示宽度补，"缺失(1)" 这种按字节算是 12 列，会顶歪整张表
row() { printf '%s %s %s\n' "$(pad "$1" 20)" "$(pad "$2" 10)" "$3"; }

# 人类可读大小。du 默认不跟随软链接，遇到 ln -s 过来的参考基因组会显示 0，
# 所以用 du -L；du 不可用（某些网络盘）时退回 ls。
human_size() {
  local p="$1" s=""
  [[ -e "$p" ]] || { printf '-'; return 0; }
  s=$(du -Lsh "$p" 2>/dev/null | awk '{print $1}' || true)
  [[ -n "$s" ]] || s=$(ls -ldh "$p" 2>/dev/null | awk '{print $5}')
  printf '%s' "${s:-?}"
}

# 参考基因组：按索引自动识别，不看文件名，也不要求软链接
find_ref() {
  local r
  r=$(discover_ref "$GENOME_DIR")
  [[ -n "$r" ]] && { printf '%s' "$r"; return 0; }
  return 1
}
