#!/usr/bin/env bash
# =============================================================================
#  deploy/30_vep_cache.sh —— VEP cache
#
#  cache 版本必须和本机 vep 的主版本号一致，否则 VEP 直接报
#  "Cannot find cache for homo_sapiens ... GRCh38"。所以这里先探测 vep 版本，
#  探测不到再用 --release 或默认值。
#
#  解压后目录是 reference/vep/homo_sapiens/<release>_GRCh38/，
#  正是 --dir_cache 需要的层级（pipeline 会自动规整路径）。
# =============================================================================
set -euo pipefail

VEP_RELEASE="${VEP_RELEASE:-112}"

detect_vep_release() {
  command -v vep >/dev/null 2>&1 || return 1
  # 不同版本输出格式不同，抓第一个像 "release 112" / "version 112" 的数字
  vep 2>&1 | grep -oiE '(release|version)[[:space:]]*[0-9]+' | head -1 | grep -oE '[0-9]+'
}

deploy_vep() {
  local release=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --release) release="$2"; shift 2 ;;
      -h|--help) echo "用法: $0 [--release N]"; return 0 ;;
      *) die "未知参数: $1" ;;
    esac
  done

  if [[ -z "$release" ]]; then
    if release=$(detect_vep_release) && [[ -n "$release" ]]; then
      log "探测到本机 vep 版本: $release"
    else
      release="$VEP_RELEASE"
      warn "没探测到 vep 版本，使用默认 ${release}（如与实际不符，用 --release 指定）"
    fi
  fi

  mkdir -p "$VEP_DIR"
  check_space 25 "$VEP_DIR"

  # 已经装好就跳过
  if [[ -d "$VEP_DIR/homo_sapiens" ]] && (( ! FORCE )); then
    log "VEP cache 已存在: $VEP_DIR/homo_sapiens"
    ls -1 "$VEP_DIR/homo_sapiens" | sed 's/^/  /'
    return 0
  fi

  local url="https://ftp.ensembl.org/pub/release-${release}/variation/indexed_vep_cache/homo_sapiens_vep_${release}_GRCh38.tar.gz"
  dl_tar "$url" "$VEP_DIR"

  if (( DRY_RUN )); then
    log "（dry-run）cache 将解压到 $VEP_DIR/homo_sapiens"
    return 0
  fi

  [[ -d "$VEP_DIR/homo_sapiens" ]] || die "解压后没有找到 $VEP_DIR/homo_sapiens，检查 tar 包结构"
  log "VEP cache 就绪: $VEP_DIR"

  # HGVS 命名需要参考序列；cache 模式下 VEP 会要求 --fasta
  local ref
  if ref=$(find_ref); then
    log "提示: 需要 HGVS 时给 vep 加 --fasta $ref"
  fi
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
  deploy_vep "$@"
fi
