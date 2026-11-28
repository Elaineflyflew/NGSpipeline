#!/usr/bin/env bash
# =============================================================================
#  src/50_stage_annot.sh —— ⑫ peak 基因组注释
#
#  annotatePeaks.pl -annStats 会把每类注释的数量和百分比写进一个统计文件，
#  直接读它比自己去数 annotated_peaks.txt 的某一列更稳（HOMER 的列顺序
#  在不同版本里变过）。
#
#  判读要点：ATAC 的 peak 应该以 distal intergenic + intron 为主（合计
#  50–70%），promoter 20–35%。promoter 超过 60% 通常是 TSS 附近高信号
#  把 peak 全吸过去了，需要回头看 ⑪ 的 TSS 图和 ⑦ 的片段分布。
# =============================================================================

stage_annot() {
  section "⑫ HOMER peak 注释"
  need annotatePeaks.pl
  need_file "$PEAK/${SAMPLE}.all.filtered.narrowPeak" "⑨ 的 peak（先跑 --stages peaks）"

  local outdir="$ANN/all"
  mkdir -p "$outdir"
  local peaks="$PEAK/${SAMPLE}.all.filtered.narrowPeak"
  local stats="$outdir/annostats.txt"

  runpipe "annotatePeaks.pl '$peaks' '$HOMER_GENOME' -annStats '$stats' > '$outdir/annotated_peaks.txt'"

  run "$PYTHON" "$PROJECT_ROOT/src/py/annot_dist.py" \
      --annstats "$stats" \
      --out "$outdir/annotation_distribution.tsv" \
      --png "$outdir/annotation_distribution.png"

  # ---- ⑫b 第二条注释路线：ChIPpeakAnno（TSS.human.GRCh38 + knownGene）
  # 和 HOMER 用不同的注释库，两条路线的 promoter 占比接近才说明 peak 位置可信。
  # 装不上就跳过，HOMER 的结果已经够用。
  if (( DRY_RUN )) || Rscript -e \
      'if (!requireNamespace("ChIPpeakAnno", quietly=TRUE)) quit(status=1)' >/dev/null 2>&1; then
    local ca="$ANN/chippeakanno"
    mkdir -p "$ca"
    run Rscript "$PROJECT_ROOT/src/py/chippeakanno.R" \
        --peaks "$peaks" --outdir "$ca" --sample "$SAMPLE"
  else
    warn "R 里没有 ChIPpeakAnno，跳过第二条注释路线（HOMER 的结果不受影响）"
  fi

  (( DRY_RUN )) && { mark_done annot; return 0; }

  # ⑫b 的数也进汇总表（加 cpa_ 前缀，和 HOMER 的键区分开）
  local ctsv="$ANN/chippeakanno/${SAMPLE}.chippeakanno_metrics.tsv"
  if [[ -s "$ctsv" ]]; then
    local k v
    while IFS=$'\t' read -r k v; do
      [[ "$k" == "key" || -z "$k" ]] && continue
      metric "12_annotate" "cpa_$k" "$v" "" ""
    done < "$ctsv"
  fi

  local name count pct
  local n=0
  while IFS=$'\t' read -r name count pct; do
    [[ "$name" == "annotation" || -z "$name" ]] && continue
    # HOMER 最后会写一行 Total，它不算一类注释。
    # 不用 ${name,,}：macOS 自带 bash 3.2 不支持大小写转换展开
    case "$name" in
      Total|total|TOTAL)
        metric "12_annotate" "peaks_annotated" "$count" "peaks" ""
        continue ;;
    esac
    n=$(( n + 1 ))
    metric "12_annotate" "$(slug "$name")_pct"   "$pct"   "%"     ""
    metric "12_annotate" "$(slug "$name")_count" "$count" "peaks" ""
  done < "$outdir/annotation_distribution.tsv"

  ok "注释完成：${n} 类（产物 $outdir/）"
  mark_done annot
}
