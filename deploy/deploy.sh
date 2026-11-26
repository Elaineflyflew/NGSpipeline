#!/usr/bin/env bash
# =============================================================================
#  deploy/deploy.sh —— 依赖安装总入口
#
#    bash deploy/deploy.sh status        # 先看现在有什么、缺什么、占多大
#    bash deploy/deploy.sh all           # 装齐（ANNOVAR 需手动，会打印指引）
#    bash deploy/deploy.sh ref bundle    # 只装指定几项
#    bash deploy/deploy.sh --dry-run all # 只打印要下载什么，不实际下
#
#  单项：env ref bundle vep funcotator deepvariant annovar
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
source "$SELF_DIR/20_gatk_bundle.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/30_vep_cache.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/40_annovar.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/50_deepvariant.sh"
# shellcheck source=/dev/null
source "$SELF_DIR/60_funcotator.sh"

usage() {
cat <<EOF
deploy.sh —— 安装 wes-pipeline 的参考数据与外部工具

用法
  bash deploy/deploy.sh [选项] <目标>...

目标
  doctor        软件体检：本机已装了什么、还缺什么（不改动任何东西）
  status        数据体检：参考数据装了哪些、占多大
  all           全部安装（ANNOVAR 需手动下载，只打印指引）
  env           只补装缺失的软件（已有的不重装；全都有则直接跳过）
  ref           参考基因组 + bwa 索引 + .fai + .dict        ~8G
  bundle        GATK 已知位点（dbsnp/hapmap/Mills/1000G）   ~1.5G
  vep           VEP cache（版本自动匹配本机 vep）           ~25G
  funcotator    Funcotator 数据源                           ~15G
  deepvariant   DeepVariant singularity 镜像                ~5G
  annovar       ANNOVAR + humandb（需学术注册，半手动）     ~10G+

选项
  --dry-run       只打印要下载/安装什么，不实际执行
  --force         已存在的也重新下载/重建
  --with-optional env 时连可选工具（multiqc/vep/apptainer）一起装
  -h, --help

关于软件
  gatk / vep / bcftools 这类工具不强制用 conda 装。deploy 会先探测
  PATH、extratools/bin、\$CONDA_PREFIX 以及常见的 miniforge/mambaforge/
  miniconda/anaconda 环境目录，已经有的直接复用，只补缺失的。

常用组合
  bash deploy/deploy.sh doctor                # 第一步永远先跑这个
  bash deploy/deploy.sh ref bundle env        # 只跑 bcftools/GATK 路线
  bash deploy/deploy.sh vep                   # 再加 VEP 注释
  bash deploy/deploy.sh --dry-run all         # 上服务器前先看清楚要下多少

下载源
  参考/资源  https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
  VEP cache  https://ftp.ensembl.org/pub/release-<N>/variation/indexed_vep_cache
  全部可用环境变量覆盖，见各子脚本顶部的 *_BASE / *_RELEASE / *_VERSION
EOF
}

deploy_env() {
  # 先体检：已有的工具一个都不动，只补缺的
  deploy_doctor --quiet >/dev/null

  local pkgs="$DOCTOR_PKGS"
  # 可选工具（multiqc / vep / apptainer ...）默认不装，--with-optional 才带上
  if (( ${ENV_WITH_OPT:-0} )); then
    pkgs="$pkgs $DOCTOR_PKGS_OPT"
    pkgs=$(printf '%s\n' $pkgs | awk 'NF && !seen[$0]++' | tr '\n' ' ')
    pkgs="${pkgs% }"
  fi

  if [[ -z "$pkgs" ]]; then
    log "必需工具已全部就绪，跳过 conda 安装"
    log "（已找到 ${DOCTOR_OK} 个工具，都在原位置复用，不会被复制或覆盖）"
    return 0
  fi

  log "缺失的 conda 包: $pkgs"

  # vep --fork 需要这个 perl 模块，装 vep 时一起带上
  case " $pkgs " in
    *" ensembl-vep "*) pkgs="$pkgs perl-parallel-forkmanager" ;;
  esac

  local mgr="conda"
  command -v mamba >/dev/null 2>&1 && mgr="mamba"

  if (( DRY_RUN )); then
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "wes"; then
      printf '[dry-run] %s install -y -n wes -c conda-forge -c bioconda %s\n' "$mgr" "$pkgs"
    else
      printf '[dry-run] %s env create -f env/environment.missing.yaml\n' "$mgr"
      printf '           其中依赖: %s\n' "$pkgs"
    fi
    return 0
  fi

  command -v "$mgr" >/dev/null 2>&1 || die "没有 conda/mamba。请先装 miniforge，或手动安装:$pkgs"

  if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "wes"; then
    log "环境 wes 已存在，只往里补装缺失的包"
    "$mgr" install -y -n wes -c conda-forge -c bioconda $pkgs
  else
    # 只把缺的写进一个临时 yaml，避免装一堆重复的
    local out="$PROJECT_ROOT/env/environment.missing.yaml"
    {
      printf 'name: wes\nchannels:\n  - conda-forge\n  - bioconda\n  - defaults\ndependencies:\n'
      local p
      for p in $pkgs; do printf '  - %s\n' "$p"; done
    } > "$out"
    log "按缺失清单创建环境 wes: $out"
    "$mgr" env create -f "$out"
  fi

  log "完成。激活后即可使用: conda activate wes"
  log "注意：本机原有环境里的工具不会被搬进来，pipeline 运行时会自动从原位置找到它们"
}

deploy_status() {
  printf '\n%s %-10s %s\n' "$(pad "项目" 20)" "状态" "路径 / 大小"
  printf '%s\n' "--------------------------------------------------------------------------"

  local ref
  if ref=$(find_ref); then
    row "参考基因组" "OK" "$ref ($(human_size "$ref"))"
    local miss=""
    [[ -s "$ref.fai" ]] || miss="$miss .fai"
    [[ -s "${ref%.fa}.dict" || -s "${ref%.fasta}.dict" || -s "${ref%.fna}.dict" ]] || miss="$miss .dict"
    [[ -s "$ref.bwt" ]] || miss="$miss .bwt"
    [[ -n "$miss" ]] && row "" "索引不全" "缺:${miss}"
  else
    row "参考基因组" "缺失" "bash deploy/deploy.sh ref"
  fi

  # .vcf 和 .vcf.gz 都算，按清单逐个判断，不靠 glob 数量
  local n=0 nm fmt=""
  for nm in $GATK_BUNDLE_FILES; do have_bundle_vcf "$nm" && n=$(( n + 1 )); done
  if (( n >= 4 )); then
    [[ -n "$(ls "$BUNDLE_DIR"/*.vcf 2>/dev/null || true)" ]] && fmt="未压缩 .vcf" || fmt=".vcf.gz"
    row "GATK bundle" "OK" "$(human_size "$BUNDLE_DIR")  (${fmt})"
  else
    row "GATK bundle" "缺失(${n}/4)" "bash deploy/deploy.sh bundle"
  fi

  if [[ -d "$VEP_DIR/homo_sapiens" ]]; then
    row "VEP cache" "OK" "$(ls -1 "$VEP_DIR/homo_sapiens" | head -1) ($(human_size "$VEP_DIR"))"
  else
    row "VEP cache" "缺失" "bash deploy/deploy.sh vep"
  fi

  if [[ -n "$(ls -A "$FUNC_DIR" 2>/dev/null)" ]]; then
    row "Funcotator" "OK" "$(human_size "$FUNC_DIR")"
  else
    row "Funcotator" "缺失" "bash deploy/deploy.sh funcotator"
  fi

  if [[ -s "$DV_DIR/deepvariant.sif" ]]; then
    row "DeepVariant" "OK" "$DV_DIR/deepvariant.sif ($(human_size "$DV_DIR/deepvariant.sif"))"
  elif command -v run_deepvariant >/dev/null 2>&1; then
    row "DeepVariant" "OK(本机)" "$(command -v run_deepvariant)"
  else
    row "DeepVariant" "缺失" "bash deploy/deploy.sh deepvariant"
  fi

  if [[ -f "$ANNOVAR_DIR/annotate_variation.pl" ]]; then
    row "ANNOVAR" "OK" "$(human_size "$ANNOVAR_DIR")"
  else
    row "ANNOVAR" "缺失" "bash deploy/deploy.sh annovar（需注册）"
  fi

  # 软件：和 doctor 用同一套探测。有工具就行，不要求存在名为 wes 的 conda 环境——
  # 用 /bioapps/bin 或别的 env 里的 gatk 完全没问题，强行要求 wes 只会误导。
  deploy_doctor --quiet >/dev/null 2>&1 || true
  if (( DOCTOR_REQ_MISSING == 0 )); then
    row "软件" "OK" "${DOCTOR_REQ_TOTAL}/${DOCTOR_REQ_TOTAL} 必需工具就绪（复用本机已装）"
  else
    row "软件" "缺失(${DOCTOR_REQ_MISSING})" "${DOCTOR_MISS_REQ# } → bash deploy/deploy.sh env"
  fi

  printf '\n项目根目录: %s\n' "$PROJECT_ROOT"
  printf '看每个工具的路径: bash deploy/deploy.sh doctor\n\n'
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

  # all 的展开顺序有讲究：先小后大，网络中断时至少小文件是完整的
  local t
  for t in "${targets[@]}"; do
    case "$t" in
      all)
        # 先体检再装：已经有 GATK 的机器不该被塞一份新的
        deploy_doctor
        ENV_WITH_OPT=1   # all = 全都要，包括可选工具
        for t in env ref bundle funcotator vep deepvariant annovar; do
          _run_one "$t"
        done
        break ;;
      *) _run_one "$t" ;;
    esac
  done

  (( DRY_RUN )) && log "（dry-run 模式，未实际下载）"
  (( ${STATUS_ONLY:-0} )) || log "完成。查看状态: bash deploy/deploy.sh status"
}

_run_one() {
  case "$1" in
    doctor)     deploy_doctor; STATUS_ONLY=1 ;;
    status)     deploy_status; STATUS_ONLY=1 ;;
    env)        deploy_env ;;
    ref)        deploy_ref ;;
    bundle)     deploy_bundle ;;
    vep)        deploy_vep ;;
    funcotator) deploy_funcotator ;;
    deepvariant) deploy_deepvariant ;;
    annovar)    deploy_annovar ;;
    *) die "未知目标: $1（可用: doctor status env ref bundle vep funcotator deepvariant annovar all）" ;;
  esac
}

main "$@"
