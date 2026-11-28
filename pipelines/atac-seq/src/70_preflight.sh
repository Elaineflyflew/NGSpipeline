#!/usr/bin/env bash
# =============================================================================
#  src/70_preflight.sh —— 预检
#
#  在真正开跑之前回答三个问题：
#    1. 需要的软件装了没有（缺了直接说缺哪个，不要跑到一半才崩）
#    2. 参考数据齐不齐
#    3. 参考数据的染色体命名和数据本身一致吗（chr17 vs 17，这个坑不报错）
# =============================================================================

# 必需工具：缺任何一个都跑不下去
REQ_TOOLS="fastqc fastp bowtie2 samtools macs2 bedtools"
# 可选工具：缺了只是少做一步，主流程照样出结果
OPT_TOOLS="multiqc alignmentSieve bamCoverage computeMatrix plotProfile plotHeatmap plotEnrichment plotFingerprint annotatePeaks.pl findMotifsGenome.pl meme-chip Rscript TOBIAS"

tool_version() {
  case "$1" in
    fastqc)    fastqc --version 2>&1 | head -1 | sed 's/^FastQC //' ;;
    fastp)     fastp --version 2>&1 | head -1 | sed 's/^fastp //' ;;
    bowtie2)   bowtie2 --version 2>&1 | head -1 | awk '{print $3}' ;;
    samtools)  samtools --version 2>&1 | head -1 | awk '{print $2}' ;;
    macs2)     macs2 --version 2>&1 | head -1 | sed 's/^macs2 //' ;;
    bedtools)  bedtools --version 2>&1 | head -1 | sed 's/^bedtools //' ;;
    *)         "$1" --version 2>&1 | head -1 ;;
  esac
}

preflight() {
  section "预检"

  local t miss=0
  printf '\n  必需工具\n'
  printf '  %s %s %s\n' "$(pad "tool" 16)" "$(pad "status" 10)" "path / version"
  for t in $REQ_TOOLS; do
    local p
    if p=$(find_tool "$t"); then
      printf '  %s %s %s  %s\n' "$(pad "$t" 16)" "$(pad "OK" 10)" "$p" "$(tool_version "$t" 2>/dev/null)"
    else
      printf '  %s %s %s\n' "$(pad "$t" 16)" "$(pad "MISS" 10)" "-"
      miss=$(( miss + 1 ))
    fi
  done
  printf '\n  必需工具: %d/%d 就绪\n' "$(( 6 - miss ))" 6

  printf '\n  可选工具\n'
  local omiss=""
  for t in $OPT_TOOLS; do
    if find_tool "$t" >/dev/null 2>&1; then
      printf '  %s %s\n' "$(pad "$t" 20)" "OK"
    else
      omiss="$omiss $t"
    fi
  done
  [[ -n "$omiss" ]] && printf '  缺失:%s\n  缺这些只会少做对应步骤，主流程照常出结果\n' "$omiss"

  printf '\n  参考数据\n'
  local f
  for f in "$BT2_INDEX.1.bt2:BT2_INDEX" "$FASTA:FASTA" "$BLACKLIST:BLACKLIST" \
           "$GTF:GTF" "$HOMER_MOTIFS:HOMER_MOTIFS" "$MOTIF_DB:MOTIF_DB"; do
    local path="${f%%:*}" name="${f#*:}"
    if [[ -z "$path" ]]; then
      printf '  %s %s\n' "$(pad "$name" 20)" "未找到（相关步骤会跳过）"
    elif [[ -s "$path" ]]; then
      printf '  %s %s  %s\n' "$(pad "$name" 20)" "OK" "$path"
    else
      printf '  %s %s  %s\n' "$(pad "$name" 20)" "MISS" "$path"
    fi
  done

  # ---- 染色体命名一致性 ----
  # Ensembl 的 17 和 UCSC 的 chr17 混用不会报错，但 peak 会全部落在错误的位置
  if [[ -s "$GTF" && -s "$BLACKLIST" ]]; then
    local g b
    g=$(awk -F'\t' 'NR==1{print ($1 ~ /^chr/ ? "chr" : "nochr")}' "$GTF")
    b=$(awk -F'\t' 'NR==1{print ($1 ~ /^chr/ ? "chr" : "nochr")}' "$BLACKLIST")
    if [[ "$g" != "$b" ]]; then
      warn "命名不一致：GTF 是 ${g}，blacklist 是 ${b} —— peak 注释和过滤都会出错"
    else
      printf '\n  染色体命名: %s（GTF 与 blacklist 一致）\n' "$g"
    fi
  fi

  if (( CHECK_ONLY && miss > 0 )); then
    printf '\n  缺 %d 个必需工具，跑 deploy/deploy.sh env 补装\n\n' "$miss"
    return 1
  fi
  printf '\n'
  return 0
}
