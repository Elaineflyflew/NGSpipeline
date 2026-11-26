#!/usr/bin/env bash
# =============================================================================
#  src/60_stage_report.sh —— ⑦ 报告
#
#  产出两样东西：
#    07_report/summary.txt   比对率/dup/深度 + 各 caller 指标 + 路线一致性 + 注释清单
#    07_report/igv.bat       可直接 igv --batch 执行，出 TP53 总览和候选位点截图
#
#  一致性对比前必须先统一口径：PASS → 拆多等位 → 左对齐 → 位点集合。
#  不做归一化，三个 caller 的交集数字毫无意义（多等位写法和左对齐位置都不同）。
# =============================================================================

# 统一口径的位点集合；放在子 shell 里兜住失败，某个 VCF 异常不该让报告挂掉
vcf_sites() {
  ( bcftools view -f 'PASS,.' "$1" 2>/dev/null \
    | bcftools norm -f "$REF" -m -any 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' 2>/dev/null \
    | sort -u ) || true
}

report_summary() {
  local R="$REP" f c
  {
    echo "=========== 流程信息 ==========="
    echo "样本      : $SAMPLE"
    echo "参考      : $REF"
    echo "参考命名  : $(ref_style)"
    echo "BED       : ${BED:-<未指定，全基因组>}"
    echo "线程      : $THREADS"
    echo "caller    : ${CALLERS[*]}"
    echo "annotator : $ANNOT_DESC"
    echo "输出      : $OUTDIR"
    echo

    echo "=========== 比对 ==========="
    [[ -s "$STAT/${SAMPLE}.flagstat.txt" ]] && head -3 "$STAT/${SAMPLE}.flagstat.txt"
    if [[ -s "$STAT/${SAMPLE}.markdup.txt" ]]; then
      echo "--- 重复标记 ---"
      grep -Ei 'READ|DUPLICATE' "$STAT/${SAMPLE}.markdup.txt" || true
    fi
    if [[ -s "$STAT/${SAMPLE}.target_depth.txt" ]]; then
      echo "--- 目标区深度 ---"
      cat "$STAT/${SAMPLE}.target_depth.txt"
    fi
    echo

    echo "=========== 变异检测 ==========="
    printf '%-12s %10s %10s %10s %8s %8s\n' "caller" "total" "PASS" "SNV" "indel" "TiTv"
    for c in "${CALLERS[@]}"; do
      f=$(final_vcf "$c")
      [[ -s "$f" ]] || continue
      local tot pass snv ind titv
      tot=$(bcftools view -H "$f" | wc -l | tr -d ' ')
      pass=$(bcftools view -H -f PASS "$f" | wc -l | tr -d ' ')
      snv=$(bcftools view -H -v snps "$f" | wc -l | tr -d ' ')
      ind=$(bcftools view -H -v indels "$f" | wc -l | tr -d ' ')
      titv=$(bcftools stats "$f" | awk -F'\t' '$1=="TSTV"{printf "%.3f", $5}')
      printf '%-12s %10s %10s %10s %8s %8s\n' "$c" "$tot" "$pass" "$snv" "$ind" "${titv:--}"
    done
    echo

    echo "=========== FILTER 分布 ==========="
    for c in "${CALLERS[@]}"; do
      f=$(final_vcf "$c"); [[ -s "$f" ]] || continue
      echo "--- $c ---"
      bcftools query -f '%FILTER\n' "$f" | sort | uniq -c | sort -rn
    done
    echo

    if (( ${#CALLERS[@]} > 1 )); then
      echo "=========== 路线一致性（归一化后位点集合）==========="
      for c in "${CALLERS[@]}"; do
        f=$(final_vcf "$c"); [[ -s "$f" ]] && vcf_sites "$f" > "$R/${c}.sites.tsv"
      done
      for c in "${CALLERS[@]}"; do
        [[ -s "$R/${c}.sites.tsv" ]] && printf '%-12s %8s\n' "$c" "$(wc -l < "$R/${c}.sites.tsv" | tr -d ' ')"
      done
      local i j A B sh oa ob
      for ((i=0; i<${#CALLERS[@]}; i++)); do
        for ((j=i+1; j<${#CALLERS[@]}; j++)); do
          A="${CALLERS[i]}"; B="${CALLERS[j]}"
          [[ -s "$R/$A.sites.tsv" && -s "$R/$B.sites.tsv" ]] || continue
          sh=$(comm -12 "$R/$A.sites.tsv" "$R/$B.sites.tsv" | wc -l | tr -d ' ')
          oa=$(comm -23 "$R/$A.sites.tsv" "$R/$B.sites.tsv" | wc -l | tr -d ' ')
          ob=$(comm -13 "$R/$A.sites.tsv" "$R/$B.sites.tsv" | wc -l | tr -d ' ')
          printf '%-10s ∩ %-10s 共享 %6s | 仅%s %6s | 仅%s %6s\n' "$A" "$B" "$sh" "$A" "$oa" "$B" "$ob"
        done
      done
      echo
    fi

    echo "=========== 注释产物 ==========="
    find "$ANN" -maxdepth 1 -type f \
        \( -name "*.vcf.gz" -o -name "*.maf" -o -name "*multianno.txt" \) \
        -exec basename {} \; 2>/dev/null | sort || ls -1 "$ANN"
    echo

    echo "=========== IGV ==========="
    echo "批处理脚本: $R/igv.bat"
    echo "用法: igv --batch $R/igv.bat"
  } > "$R/summary.txt"
}

report_igv() {
  local R="$REP" c f
  {
    echo "# 先确认 IGV 加载的基因组版本与 BAM header 一致"
    echo new
    echo "genome $REF"
    echo "load $MDUP"
    local vcfs=()
    for c in "${CALLERS[@]}"; do
      f=$(final_vcf "$c"); [[ -s "$f" ]] && vcfs+=("$f")
    done
    (( ${#vcfs[@]} )) && echo "load $(IFS=,; echo "${vcfs[*]}")"
    echo "sort base"
    echo "maxPanelHeight 700"
    echo "# TP53 总览"
    echo "goto $IGV_REGION"
    echo "snapshot $R/igv.TP53.overview.png"

    # 候选位点：从第一个 caller 的 PASS 里按 QUAL 挑，纯合/杂合/indel 各最多 2 个
    f=$(final_vcf "${CALLERS[0]}")
    if [[ -s "$f" ]]; then
      bcftools query -i 'FILTER="PASS"' \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%DP\t%AD]\n' "$f" \
        | awk -v max="$IGV_MAX" '
            $7 >= 10 {
              split($6, gt, "[|/]"); split($8, ad, ",")
              af = (ad[1]+ad[2]) > 0 ? ad[2]/(ad[1]+ad[2]) : 0
              type = (length($3) != length($4)) ? "indel" :
                     (gt[1]==gt[2] && af>=0.85) ? "hom" :
                     (gt[1]!=gt[2] && af>=0.20 && af<=0.80) ? "het" : ""
              if (type != "") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"type
            }' \
        | sort -k5,5gr \
        | awk -v max="$IGV_MAX" '{ if (++n[$6] <= 2 && c < max) { print; c++ } }' \
        > "$R/igv.candidates.tsv"

      local chrom pos ra aa qual type
      while IFS=$'\t' read -r chrom pos ra aa qual type; do
        echo "# $type $chrom:$pos $ra>$aa QUAL=$qual"
        echo "goto ${chrom}:$(( pos > 60 ? pos - 60 : 1 ))-$(( pos + 60 ))"
        echo "snapshot $R/igv.${type}.${chrom}.${pos}.png"
      done < "$R/igv.candidates.tsv"
    fi
  } > "$R/igv.bat"
}

stage_report() {
  section "⑦ 生成报告"
  if (( DRY_RUN )); then
    log "（dry-run：报告依赖真实 VCF，跳过）"
    return 0
  fi

  # FastQC 报告多了就用 MultiQC 合成一份，省得逐个开 html
  if command -v multiqc >/dev/null 2>&1 && [[ -n "$(ls -A "$QC" 2>/dev/null)" ]]; then
    ( cd "$OUTDIR" && multiqc -f -o "$REP/multiqc" 01_qc 02_clean >/dev/null 2>&1 ) \
      && log "MultiQC: $REP/multiqc/multiqc_report.html" \
      || warn "MultiQC 运行失败，忽略"
  fi

  report_summary
  report_igv
  mark_done report
  log "报告: $REP/summary.txt"
}
