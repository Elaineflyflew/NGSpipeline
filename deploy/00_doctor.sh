#!/usr/bin/env bash
# =============================================================================
#  deploy/00_doctor.sh —— 软件体检
#
#  先看机器上已经有什么，再决定装什么。理由很实际：
#  多数服务器上 gatk / bcftools / samtools 早就装好了（只是可能不在 PATH 里），
#  再 conda 装一遍既浪费 3 G 空间，又可能装上另一个版本，让结果对不上。
#
#  探测顺序：PATH → extratools/bin → $CONDA_PREFIX → 常见 conda/mamba 环境的 bin
#
#    bash deploy/00_doctor.sh            # 完整表格
#    bash deploy/00_doctor.sh --quiet    # 只输出缺失清单（给 deploy.sh env 用）
#    bash deploy/00_doctor.sh --json     # 机器可读（后续接 CI 用）
# =============================================================================
set -euo pipefail

# 清单格式：组|可执行文件|conda 包|说明|必需(1/0)
# conda 包写 - 表示不走 conda（需手动装或本来就该有）
doctor_rows() {
cat <<'EOF'
质控|fastqc|fastqc|FastQC 原始数据质控|1
质控|fastp|fastp|接头/低质量修剪|1
质控|multiqc|multiqc|MultiQC 报告汇总|0
比对|bwa|bwa|bwa mem 比对|1
比对|samtools|samtools|BAM 排序/去重/统计|1
检测|bcftools|bcftools|bcftools 路线 + 统计|1
检测|bgzip|htslib|bgzip 压缩|1
检测|tabix|htslib|tabix 索引|1
检测|gatk|gatk4|GATK 路线 + Funcotator|1
注释|vep|ensembl-vep|VEP 注释|0
注释|table_annovar.pl|-|ANNOVAR（需学术注册）|0
容器|singularity|apptainer|DeepVariant 镜像|0
容器|apptainer|apptainer|DeepVariant 镜像（替代 singularity）|0
检测|run_deepvariant|-|DeepVariant 本机版（有 sif 就不需要）|0
EOF
}

# 探测结果写进这几个全局变量，deploy_env 直接用
DOCTOR_OK=0
DOCTOR_MISS_REQ=""
DOCTOR_MISS_OPT=""
DOCTOR_MISS_OPT_NOCONDA=""   # 可选工具里 conda 装不了的那些（如 run_deepvariant）
DOCTOR_PKGS=""          # 缺失的"必需"工具对应的 conda 包
DOCTOR_PKGS_OPT=""      # 缺失的"可选"工具对应的 conda 包
DOCTOR_REQ_TOTAL=0      # 必需工具总数
DOCTOR_REQ_MISSING=0    # 缺失的必需工具数（status 用它判断软件是否就绪）

deploy_doctor() {
  local quiet=0 json=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --quiet) quiet=1; shift ;;
      --json)  json=1; shift ;;
      -h|--help)
        echo "用法: $0 [--quiet|--json]"
        echo "  --quiet  只输出缺失的 conda 包（供 deploy.sh env 使用）"
        return 0 ;;
      *) die "未知参数: $1" ;;
    esac
  done

  local grp tool pkg desc req p
  local req_total=0 req_missing=0

  (( quiet || json )) || {
    printf '\n软件体检（先看已有什么，再决定装什么）\n'
    printf '%s %-10s %s\n' "$(pad "工具" 20)" "状态" "位置 / 用途"
    printf '%s\n' "--------------------------------------------------------------------------"
  }

  while IFS='|' read -r grp tool pkg desc req; do
    [[ -z "${tool:-}" ]] && continue
    if p=$(find_tool "$tool"); then
      DOCTOR_OK=$(( DOCTOR_OK + 1 ))
      (( quiet || json )) || row "$tool" "OK" "$p  ($desc)"
    else
      if [[ "$req" == 1 ]]; then
        req_missing=$(( req_missing + 1 ))
        DOCTOR_MISS_REQ="$DOCTOR_MISS_REQ $tool"
        [[ "$pkg" != "-" ]] && DOCTOR_PKGS="$DOCTOR_PKGS $pkg"
        (( quiet || json )) || row "$tool" "缺失" "-  ($desc)"
      else
        DOCTOR_MISS_OPT="$DOCTOR_MISS_OPT $tool"
        if [[ "$pkg" != "-" ]]; then
          DOCTOR_PKGS_OPT="$DOCTOR_PKGS_OPT $pkg"
        else
          DOCTOR_MISS_OPT_NOCONDA="$DOCTOR_MISS_OPT_NOCONDA $tool"
        fi
        (( quiet || json )) || row "$tool" "缺失(可选)" "-  ($desc)"
      fi
    fi
    [[ "$req" == 1 ]] && req_total=$(( req_total + 1 ))
  done < <(doctor_rows)

  # htslib 同时提供 bgzip 和 tabix，去重
  DOCTOR_PKGS=$(printf '%s\n' $DOCTOR_PKGS | awk 'NF && !seen[$0]++' | tr '\n' ' ')
  DOCTOR_PKGS="${DOCTOR_PKGS% }"
  DOCTOR_PKGS_OPT=$(printf '%s\n' $DOCTOR_PKGS_OPT | awk 'NF && !seen[$0]++' | tr '\n' ' ')
  DOCTOR_PKGS_OPT="${DOCTOR_PKGS_OPT% }"

  local have=$(( req_total - req_missing ))
  DOCTOR_REQ_TOTAL=$req_total
  DOCTOR_REQ_MISSING=$req_missing

  if (( json )); then
    printf '{"ok":%d,"missing_required":"%s","missing_optional":"%s","packages":"%s"}\n' \
      "$DOCTOR_OK" "${DOCTOR_MISS_REQ# }" "${DOCTOR_MISS_OPT# }" "$DOCTOR_PKGS"
    return 0
  fi
  (( quiet )) && { printf '%s\n' "$DOCTOR_PKGS"; return 0; }

  printf '\n必需工具: %d/%d 就绪\n' "$have" "$req_total"
  if [[ -z "$DOCTOR_MISS_REQ" ]]; then
    printf '缺失的必需工具: 无 —— 环境已满足，不需要 conda 装任何东西\n'
  else
    printf '缺失的必需工具:%s\n' "$DOCTOR_MISS_REQ"
    printf '建议: bash deploy/deploy.sh env     # 只装缺失的，已有的不动\n'
  fi
  if [[ -n "$DOCTOR_MISS_OPT" ]]; then
    printf '可选工具缺失:%s（不影响主流程）\n' "$DOCTOR_MISS_OPT"
    # 别给空头承诺：run_deepvariant 这类 conda 里没有，--with-optional 也装不上
    if [[ -n "$DOCTOR_PKGS_OPT" ]]; then
      printf '      要一起装: bash deploy/deploy.sh env --with-optional\n'
    fi
    if [[ -n "$DOCTOR_MISS_OPT_NOCONDA" ]]; then
      printf '      需手动安装（conda 里没有）:%s\n' "$DOCTOR_MISS_OPT_NOCONDA"
    fi
  fi
  printf '\n提示: 探测范围包括 PATH、extratools/bin、$CONDA_PREFIX 以及常见的\n'
  printf '      miniforge/mambaforge/miniconda/anaconda 环境目录。\n'
  printf '      已经在别的 env 里装好的工具会被直接复用，不会被复制或重装。\n'
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
  deploy_doctor "$@"
fi
