#!/usr/bin/env bash
# =============================================================================
#  deploy/deploy.sh —— 依赖体检与安装
#
#  这个流程默认【不下载任何数据】：参考基因组、索引、blacklist、GTF、motif 库
#  都从 config/genomes/<name>.conf 里探测，服务器上已有的直接复用。
#
#  所以 deploy 只做两件事：
#    doctor  看装了什么、缺什么
#    env     缺的软件用 conda 补（不重建环境，不覆盖已有版本）
#    status  看参考数据都指到哪了
#
#  用法
#    bash deploy/deploy.sh doctor
#    bash deploy/deploy.sh status
#    bash deploy/deploy.sh env
#    bash deploy/deploy.sh all
# =============================================================================
set -Eeuo pipefail

# shellcheck source=lib.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"

usage() {
cat <<EOF
atac-seq deploy

用法: bash deploy/deploy.sh <命令>

  doctor    软件体检：必需/可选工具装了什么、缺什么
  status    参考数据：每个路径指向哪里、大小多少
  env       只补装缺失的必需工具（不重建环境）
  all       env + 打印后续步骤

环境变量
  ENV_NAME   conda 环境名（默认 wes，和课程其他流程共用）
  CONDA_EXE  conda 可执行文件（默认自动找）
EOF
}

cmd_status() {
  section "参考数据"
  printf '  %s %s\n' "$(pad "genome" 12)" "$GENOME_LABEL ($GENOME_NAME)"
  printf '  %s %s\n' "$(pad "config" 12)" "$PROJECT_ROOT/config/genomes/${GENOME}.conf"
  echo
  printf '  %s %s %s\n' "$(pad "item" 14)" "$(pad "size" 10)" "path"
  local pair name path
  for pair in "BT2_INDEX:$BT2_INDEX" "FASTA:$FASTA" "BLACKLIST:$BLACKLIST" \
              "GTF:$GTF" "HOMER_MOTIFS:$HOMER_MOTIFS" "MOTIF_DB:$MOTIF_DB"; do
    name="${pair%%:*}"; path="${pair#*:}"
    if [[ -z "$path" ]]; then
      printf '  %s %s %s\n' "$(pad "$name" 14)" "$(pad "-" 10)" "未找到"
    elif [[ -e "$path" ]]; then
      printf '  %s %s %s\n' "$(pad "$name" 14)" "$(pad "$(human_size "$path")" 10)" "$path"
    else
      printf '  %s %s %s\n' "$(pad "$name" 14)" "$(pad "-" 10)" "路径不存在: $path"
    fi
  done
  echo
  echo "  这些路径全部来自 config/genomes/${GENOME_NAME}.conf，按候选顺序取第一个存在的。"
  echo "  想换位置：export BT2_INDEX=... 或 --bt2-index ..."
}

cmd_env() {
  local missing="" t
  for t in $REQ_TOOLS; do
    find_tool "$t" >/dev/null 2>&1 || missing="$missing $t"
  done

  if [[ -z "$missing" ]]; then
    log "必需工具全部就绪，什么都不用装"
    return 0
  fi

  log "缺少:$missing"

  local pkgs="" p
  for t in $missing; do
    if conda_has "$t"; then
      p=$(pkg_for "$t")
      case " $pkgs " in *" $p "*) ;; *) pkgs="$pkgs $p" ;; esac
    else
      warn "$t 不在 conda 默认渠道里，需要手动装"
    fi
  done

  if [[ -z "$pkgs" ]]; then
    warn "没有可自动安装的包，请手动处理上面列出的工具"
    return 0
  fi

  local conda_exe="${CONDA_EXE:-}"
  [[ -z "$conda_exe" ]] && conda_exe=$(command -v conda 2>/dev/null || true)
  if [[ -z "$conda_exe" ]]; then
    warn "找不到 conda，请手动安装:$pkgs"
    return 0
  fi

  local env_name="${ENV_NAME:-wes}"
  echo
  log "安装到 conda 环境: $env_name"
  log "包:$pkgs"
  run "$conda_exe" install -y -n "$env_name" -c conda-forge -c bioconda $pkgs
}

cmd_doctor() {
  preflight
}

case "${1:-}" in
  doctor) cmd_doctor ;;
  status) cmd_status ;;
  env)    cmd_env ;;
  all)    cmd_env; echo; cmd_status ;;
  -h|--help|"") usage ;;
  *) usage >&2; die "未知命令: $1" ;;
esac
