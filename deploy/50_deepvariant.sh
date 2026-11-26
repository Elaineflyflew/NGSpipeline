#!/usr/bin/env bash
# =============================================================================
#  deploy/50_deepvariant.sh —— DeepVariant singularity 镜像
#
#  有 singularity/apptainer：直接从 docker hub 拉成 sif（约 5G）
#  只有 docker：先 pull 再从 daemon 构建
#  都没有：打印安装建议，不报错
#
#  镜像版本要固定，别用 latest —— 换版本会静默改变 call 结果，
#  论文/报告里没法复现。
# =============================================================================
set -euo pipefail

DV_VERSION="${DV_VERSION:-1.6.1}"
DV_IMAGE="docker://google/deepvariant:${DV_VERSION}"

deploy_deepvariant() {
  local version="$DV_VERSION"
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --version) version="$2"; shift 2 ;;
      -h|--help) echo "用法: $0 [--version 1.6.1]"; return 0 ;;
      *) die "未知参数: $1" ;;
    esac
  done
  DV_IMAGE="docker://google/deepvariant:${version}"

  mkdir -p "$DV_DIR"
  local sif="$DV_DIR/deepvariant.sif"

  if [[ -s "$sif" ]] && (( ! FORCE )); then
    log "已存在: $sif ($(du -h "$sif" | cut -f1))"
    return 0
  fi

  check_space 10 "$DV_DIR"

  if (( DRY_RUN )); then
    printf '[dry-run] singularity build %s %s\n' "$sif" "$DV_IMAGE"
    return 0
  fi

  local eng=""
  for e in singularity apptainer; do
    command -v "$e" >/dev/null 2>&1 && { eng="$e"; break; }
  done

  if [[ -n "$eng" ]]; then
    log "用 $eng 构建镜像（约 5G，10-20 分钟）"
    "$eng" build "$sif" "$DV_IMAGE" || die "镜像构建失败，检查网络或换成 --version 指定其他版本"
  elif command -v docker >/dev/null 2>&1; then
    log "没有 singularity，改用 docker 拉取后再转换"
    docker pull "google/deepvariant:${version}" || die "docker pull 失败"
    die "已 pull 到本地 docker，但转换需要 singularity。建议：conda install -c conda-forge singularity"
  else
    cat <<EOF

  既没有 singularity 也没有 docker，无法下载镜像。两种办法：
    1. conda install -c conda-forge apptainer   # 再重跑本脚本
    2. 在别的机器 build 好 deepvariant.sif，scp 到：
         $sif

  没有镜像也能跑 DeepVariant：用 conda 装
         conda install -c bioconda deepvariant
  此时 pipeline 会自动调用本机 run_deepvariant。

EOF
    return 0
  fi

  [[ -s "$sif" ]] || die "构建后文件为空: $sif"
  log "DeepVariant 就绪: $sif"
  log "pipeline.conf 会自动探测到它，不需要改配置"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
  deploy_deepvariant "$@"
fi
