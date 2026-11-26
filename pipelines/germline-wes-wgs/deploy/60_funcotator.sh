#!/usr/bin/env bash
# =============================================================================
#  deploy/60_funcotator.sh —— GATK Funcotator 数据源
#
#  解压后目录结构：reference/funcotator/<数据源版本>/dbsnp, gencode, ... 
#  pipeline.conf 会自动取第一个子目录作为 --data-sources-path。
#
#  Broad 的 tarball 名字会随版本变。如果默认文件名 404，用下面的命令查：
#      gsutil ls gs://gcp-public-data--broad-references/hg38/v0/funcotator/
#  然后 --url 覆盖。
# =============================================================================
set -euo pipefail

BROAD_FUNC_DIR="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/funcotator"
FUNCOTATOR_TARBALL="${FUNCOTATOR_TARBALL:-funcotator_dataSources.v1.8.hg38.20230908g.tar.gz}"

deploy_funcotator() {
  local url=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --url)  url="$2"; shift 2 ;;
      --file) FUNCOTATOR_TARBALL="$(basename "$2")"; url="$2"; shift 2 ;;
      -h|--help) echo "用法: $0 [--url URL | --file 本地tar.gz]"; return 0 ;;
      *) die "未知参数: $1" ;;
    esac
  done
  [[ -z "$url" ]] && url="$BROAD_FUNC_DIR/$FUNCOTATOR_TARBALL"

  mkdir -p "$FUNC_DIR"
  check_space 15 "$FUNC_DIR"

  # 已经有数据源目录就跳过
  local d
  for d in "$FUNC_DIR"/*/; do
    [[ -d "$d/dbsnp" || -d "$d/gencode" ]] && { log "Funcotator 数据源已存在: ${d%/}"; return 0; }
  done

  if [[ "$url" == /* ]]; then
    log "使用本地文件: $url"
    (( DRY_RUN )) || tar -xzf "$url" -C "$FUNC_DIR"
  else
    dl_tar "$url" "$FUNC_DIR"
  fi

  if (( DRY_RUN )); then
    log "（dry-run）数据源将解压到 $FUNC_DIR/"
    return 0
  fi

  log "Funcotator 数据源就绪: $FUNC_DIR"
  ls -1 "$FUNC_DIR" | sed 's/^/  /'
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
  deploy_funcotator "$@"
fi
