#!/usr/bin/env bash
# =============================================================================
#  deploy/00_doctor.sh —— 软件体检
#
#  只探测、不安装、不改 PATH。回答两个问题：
#    1. 这条链要的工具，本机有几个？
#    2. 每个工具在哪（PATH 里 / extratools/bin / 某个没激活的 conda env）？
#
#  探测逻辑和 src/00_common.sh 的 find_tool 是同一个，所以"doctor 说有"
#  和"跑的时候能找到"一定一致。
# =============================================================================

# 必需 / 可选工具清单。key 是命令名，value 是用途（打印用）
DOCTOR_REQ="fastqc fastp bowtie2 samtools bedtools macs2 bamCoverage computeMatrix plotProfile plotFingerprint"
DOCTOR_OPT="multiqc annotatePeaks.pl findMotifsGenome.pl Rscript igv pyGenomeTracks unzip aria2c"

deploy_doctor() {
  local quiet=0
  [[ "${1:-}" == "--quiet" ]] && quiet=1

  local ok_n=0 ok_req=0 miss_n=0 miss_list="" p

  (( quiet )) || {
    printf '\n%s %s %s\n' "$(pad "工具" 20)" "$(pad "状态" 12)" "路径"
    printf '%s\n' "--------------------------------------------------------------------------"
  }

  for p in $DOCTOR_REQ $DOCTOR_OPT; do
    local is_req=0
    case " $DOCTOR_REQ " in *" $p "*) is_req=1 ;; esac

    if command -v "$p" >/dev/null 2>&1; then
      (( quiet )) || row "$p" "OK" "$(command -v "$p")"
      ok_n=$(( ok_n + 1 )); (( is_req )) && ok_req=$(( ok_req + 1 ))
    elif find_tool "$p" >/dev/null 2>&1; then
      # 找到但不在 PATH 里——运行时 need() 会临时加进去
      (( quiet )) || row "$p" "OK(未激活)" "$(find_tool "$p")"
      ok_n=$(( ok_n + 1 )); (( is_req )) && ok_req=$(( ok_req + 1 ))
    elif (( is_req )); then
      (( quiet )) || row "$p" "缺失" "→ deploy/deploy.sh env"
      miss_n=$(( miss_n + 1 ))
      miss_list="$miss_list $p"
    else
      (( quiet )) || row "$p" "可选缺失" "用到再装"
    fi
  done

  DOCTOR_OK=$ok_n
  DOCTOR_REQ_OK=$ok_req
  DOCTOR_REQ_MISSING=$miss_n
  DOCTOR_MISS_REQ="$miss_list"

  (( quiet )) || {
    printf '\n必需工具 %d 个就绪，缺 %d 个；可选工具另找到 %d 个\n' \
      "$ok_req" "$miss_n" "$(( ok_n - ok_req ))"
    if (( miss_n > 0 )); then
      printf '补装: bash deploy/deploy.sh env\n'
    else
      printf '无需安装。注意：HOMER 的 annotatePeaks.pl / findMotifsGenome.pl\n'
      printf '要额外 configureHomer.pl -install hg38，见 bash deploy/deploy.sh homer\n'
    fi
  }
  return 0
}

# 按缺失清单生成 conda 包名。samtools/bedtools 这类 bioconda 包名和命令名一致；
# deeptools 提供 bamCoverage/computeMatrix/plotProfile/plotFingerprint 四个命令。
deploy_env_pkgs() {
  local pkgs="" p
  for p in $DOCTOR_MISS_REQ; do
    case "$p" in
      bamCoverage|computeMatrix|plotProfile|plotFingerprint) pkgs="$pkgs deeptools" ;;
      *) pkgs="$pkgs $p" ;;
    esac
  done
  # 去重
  printf '%s\n' $pkgs | awk 'NF && !seen[$0]++' | tr '\n' ' '
}
