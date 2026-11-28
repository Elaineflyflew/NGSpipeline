#!/usr/bin/env bash
# =============================================================================
#  deploy/deploy.sh —— 依赖安装总入口
#
#    bash deploy/deploy.sh status         # 先看现在有什么、缺什么、占多大
#    bash deploy/deploy.sh all            # 装齐
#    bash deploy/deploy.sh genome bt2     # 只装指定几项
#    bash deploy/deploy.sh --dry-run all  # 只打印要下载什么，不实际下
#
#  单项：env genome bt2 annotation blacklist homer atacqc igv
# =============================================================================
set -euo pipefail

SELF_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$SELF_DIR/lib.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/00_doctor.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/10_ref_genome.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/15_annotation.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/20_blacklist.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/30_homer.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/40_atacqc.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/50_igv.sh"

usage() {
cat <<EOF
deploy.sh —— 安装 ATAC-seq 流水线的参考数据与外部工具

用法
  bash deploy/deploy.sh [选项] <目标>...

目标
  doctor        软件体检：本机已装了什么、还缺什么（不改动任何东西）
  status        数据体检：参考数据装了哪些、占多大
  all           全部安装
  env           只补装缺失的软件（已有的不重装；全都有则跳过）
  genome        参考基因组 fasta + .fai                     ~1.0G
  bt2           bowtie2 索引                                ~4.5G（建索引约 1 小时）
  annotation    GENCODE v${GENCODE_RELEASE} GTF             ~1.5G
  blacklist     ENCODE hg38 blacklist v2                    ~1M
  homer         HOMER 本体 + hg38 数据                      ~3G
  atacqc        R / Bioconductor ATACseqQC（可选）          ~1G
  igv           IGV 解压到 extratools/（可选）              ~200M

选项
  --dry-run       只打印要下载/安装什么，不实际执行
  --force         已存在的也重新下载/重建
  --with-optional env 时连 multiqc / igv / pygenometracks 一起装
  -h, --help

关于软件
  fastqc / fastp / bowtie2 / samtools / bedtools / macs2 / deeptools 不强制用
  conda 装。deploy 会先探测 PATH、extratools/bin、\$CONDA_PREFIX 以及常见的
  miniforge/mambaforge/miniconda/anaconda 环境目录，已经有的直接复用，只补缺失的。

常用组合
  bash deploy/deploy.sh doctor              # 第一步永远先跑这个
  bash deploy/deploy.sh genome bt2 annotation blacklist
  bash deploy/deploy.sh homer               # 要 ⑫⑬ 才装
  bash deploy/deploy.sh --dry-run all       # 上服务器前先看清楚要下多少

下载源
  基因组     ${UCSC_BASE}
  GTF        ${GENCODE_BASE}/release_${GENCODE_RELEASE}
  黑名单     GitHub Boyle-Lab/Blacklist
  HOMER      ${HOMER_BASE}
  IGV        ${IGV_BASE}
  全部可用环境变量覆盖，见各子脚本顶部的 *_BASE / *_RELEASE / *_VERSION
EOF
}

deploy_env() {
  # 先体检：已有的工具一个都不动，只补缺的
  deploy_doctor --quiet >/dev/null
  local pkgs; pkgs=$(deploy_env_pkgs)

  if [[ -z "${pkgs// /}" ]]; then
    log "必需工具已全部就绪，跳过 conda 安装"
    return 0
  fi

  log "缺失的 conda 包: $pkgs"

  local mgr="conda"
  command -v mamba >/dev/null 2>&1 && mgr="mamba"

  if (( DRY_RUN )); then
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "atac"; then
      printf '[dry-run] %s install -y -n atac -c conda-forge -c bioconda %s\n' "$mgr" "$pkgs"
    else
      printf '[dry-run] %s env create -f env/environment.missing.yaml\n' "$mgr"
      printf '           其中依赖: %s\n' "$pkgs"
    fi
    return 0
  fi

  command -v "$mgr" >/dev/null 2>&1 || die "没有 conda/mamba。请先装 miniforge，或手动安装:$pkgs"

  if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "atac"; then
    log "环境 atac 已存在，只往里补装缺失的包"
    "$mgr" install -y -n atac -c conda-forge -c bioconda $pkgs
  else
    local out="$PROJECT_ROOT/env/environment.missing.yaml"
    {
      printf 'name: atac\nchannels:\n  - conda-forge\n  - bioconda\n  - defaults\ndependencies:\n'
      local p
      for p in $pkgs; do printf '  - %s\n' "$p"; done
    } > "$out"
    log "按缺失清单创建环境 atac: $out"
    "$mgr" env create -f "$out"
  fi

  log "完成。激活后即可使用: conda activate atac"
  log "注意：本机原有环境里的工具不会被搬进来，pipeline 运行时会自动从原位置找到它们"
}

deploy_status() {
  printf '\n%s %s %s\n' "$(pad "项目" 20)" "$(pad "状态" 12)" "路径 / 大小"
  printf '%s\n' "--------------------------------------------------------------------------"

  local fa; fa=$(find_genome)
  if [[ -n "$fa" ]]; then
    row "参考基因组" "OK" "$fa ($(human_size "$fa"))"
    [[ -s "$fa.fai" ]] || row "" "索引不全" "缺 .fai（samtools faidx）"
  else
    row "参考基因组" "缺失" "bash deploy/deploy.sh genome"
  fi

  local bt2; bt2=$(find_bt2)
  if [[ -n "$bt2" ]]; then
    row "bowtie2 索引" "OK" "$bt2 ($(human_size "$BT2_DIR"))"
  else
    row "bowtie2 索引" "缺失" "bash deploy/deploy.sh bt2"
  fi

  local gtf; gtf=$(discover_one "$ANN_DIR" '*.gtf' '*.gtf.gz')
  if [[ -n "$gtf" ]]; then
    row "基因注释 GTF" "OK" "$gtf ($(human_size "$gtf"))"
  else
    row "基因注释 GTF" "缺失" "bash deploy/deploy.sh annotation"
  fi

  local bl; bl=$(discover_one "$BL_DIR" '*.bed' '*.bed.gz')
  if [[ -n "$bl" ]]; then
    row "ENCODE 黑名单" "OK" "$bl ($(wc -l < "$bl" 2>/dev/null || echo 0) 个区间)"
  else
    row "ENCODE 黑名单" "缺失" "bash deploy/deploy.sh blacklist"
  fi

  # HOMER 有两种存在方式：deploy 装在项目内 extratools/homer，或者本来就装在
  # 系统/conda 里（课程镜像是 /home/biotools/miniforge3/envs/motif/share/homer）。
  # 只看项目内那份会把"已经装好、⑫⑬ 一直在用"的 HOMER 报成缺失。
  local hroot="" hbin="" cand=""
  if [[ -d "$HOMER_DIR/data/genomes/hg38" ]]; then
    hroot="$HOMER_DIR"
  else
    hbin=$(find_tool annotatePeaks.pl 2>/dev/null || true)
    if [[ -n "$hbin" ]]; then
      # conda 布局：<prefix>/bin/annotatePeaks.pl → <prefix>/share/homer
      # 经典布局：<homer>/bin/annotatePeaks.pl    → <homer>
      for cand in \
          "$(dirname "$(dirname "$(readlink -f "$hbin")")")/share/homer" \
          "$(dirname "$(dirname "$(readlink -f "$hbin")")")"; do
        [[ -d "$cand/data/genomes/hg38" ]] && { hroot="$cand"; break; }
      done
    fi
  fi
  if [[ -n "$hroot" ]]; then
    row "HOMER + hg38" "OK" "$hroot ($(human_size "$hroot/data/genomes/hg38"))"
  else
    row "HOMER + hg38" "缺失" "bash deploy/deploy.sh homer"
  fi

  if (( DRY_RUN )) || Rscript -e \
      'if (!requireNamespace("ATACseqQC", quietly=TRUE)) quit(status=1)' >/dev/null 2>&1; then
    row "ATACseqQC" "OK" "R 包可用"
  else
    row "ATACseqQC" "可选缺失" "bash deploy/deploy.sh atacqc"
  fi

  if find_tool igv >/dev/null 2>&1 || find_tool pyGenomeTracks >/dev/null 2>&1; then
    row "截图工具" "OK" "$(find_tool igv 2>/dev/null || find_tool pyGenomeTracks)"
  else
    row "截图工具" "可选缺失" "bash deploy/deploy.sh igv"
  fi

  deploy_doctor --quiet >/dev/null 2>&1 || true
  if (( DOCTOR_REQ_MISSING == 0 )); then
    row "软件" "OK" "${DOCTOR_REQ_OK} 个必需工具就绪（复用本机已装）"
  else
    row "软件" "缺失(${DOCTOR_REQ_MISSING})" "${DOCTOR_MISS_REQ# } → bash deploy/deploy.sh env"
  fi

  printf '\n项目根目录: %s\n' "$PROJECT_ROOT"
  printf '看每个工具的路径: bash deploy/deploy.sh doctor\n\n'
}

_run_one() {
  case "$1" in
    doctor)     deploy_doctor; STATUS_ONLY=1 ;;
    status)     deploy_status; STATUS_ONLY=1 ;;
    env)        deploy_env ;;
    genome)     deploy_genome ;;
    bt2)        deploy_bt2 ;;
    annotation) deploy_annotation ;;
    blacklist)  deploy_blacklist ;;
    homer)      deploy_homer ;;
    atacqc)     deploy_atacqc ;;
    igv)        deploy_igv ;;
    *) die "未知目标: $1（可用: doctor status env genome bt2 annotation blacklist homer atacqc igv all）" ;;
  esac
}

main() {
  local -a targets=()
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --dry-run) DRY_RUN=1; shift ;;
      --force)   FORCE=1; shift ;;
      --with-optional) ENV_WITH_OPT=1; shift ;;
      -h|--help) usage; exit 0 ;;
      -*)        die "未知选项: $1" ;;
      *)         targets+=("$1"); shift ;;
    esac
  done
  [[ ${#targets[@]} -gt 0 ]] || { usage; exit 2; }

  local t
  for t in "${targets[@]}"; do
    if [[ "$t" == "all" ]]; then
      deploy_doctor
      # 顺序有讲究：先小后大，网络中断时至少小文件是完整的
      for t in env blacklist annotation genome bt2 homer atacqc igv; do
        _run_one "$t"
      done
      break
    fi
    _run_one "$t"
  done

  (( DRY_RUN )) && log "（dry-run 模式，未实际下载）"
  (( ${STATUS_ONLY:-0} )) || log "完成。查看状态: bash deploy/deploy.sh status"
}

main "$@"
