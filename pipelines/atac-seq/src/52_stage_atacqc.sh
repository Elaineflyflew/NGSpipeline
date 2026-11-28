#!/usr/bin/env bash
# =============================================================================
#  src/52_stage_atacqc.sh —— ⑭ ATACseqQC 交叉验证
#
#  用**未偏移**的 srt.bam：TSSE 和片段分布都基于真实片段长度，而
#  alignmentSieve 挪了坐标却不会同步改 TLEN，两者会自相矛盾。
#
#  这一步的价值是"换个实现再算一遍"：⑪ 用 deepTools 算 TSS 富集，
#  这里用 Bioconductor 的 TSSEscore，两个数接近说明前面的矩阵没算错。
#  装不上 ATACseqQC 不该拦住主链——⑦⑪ 已经给出同口径的数值。
# =============================================================================

stage_atacqc() {
  section "⑭ ATACseqQC 片段分布 / TSSE"
  need Rscript

  # 服务器上常有两套 R（系统 R 与 conda R），ATACseqQC 和 BSgenome 未必装在
  # 同一套里。分两档找：优先"核心包 + BSgenome 都齐"的（能跑 CTCF 足迹），
  # 退一步只要有核心包的（足迹那步会自动跳过）。都不行就降级，不拦主链。
  # 想强行指定某一套 R，设环境变量 RSCRIPT=/usr/bin/Rscript。
  local rbin=""
  if (( ! DRY_RUN )); then
    local core=(ATACseqQC TxDb.Hsapiens.UCSC.hg38.knownGene Rsamtools GenomicAlignments)
    if [[ -n "${RSCRIPT:-}" ]]; then
      rbin="$RSCRIPT"
    else
      rbin=$(pick_rscript "${core[@]}" BSgenome.Hsapiens.UCSC.hg38) \
        || rbin=$(pick_rscript "${core[@]}") || rbin=""
    fi
    if [[ -z "$rbin" || ! -x "$rbin" ]]; then
      warn "没有一套 R 同时具备 ATACseqQC + TxDb.hg38，跳过（bash deploy/deploy.sh atacqc 可补）"
      metric "14_atacqc" "atacqc_available" "0" "" ""
      mark_done atacqc
      return 0
    fi
    local rver
    rver=$("$rbin" -e 'cat(R.version.string)' 2>/dev/null || true)
    log "用 $rbin（${rver:-R ?}）"
  else
    rbin="Rscript"
  fi

  need_file "$BAM_SRT" "④ 的 srt.bam（先跑 --stages bam）"
  # 足迹用 ④ 的 srt2.bam（= alignmentSieve 输出再排序 + 建索引），
  # 不要用 shift.bam：那份没排序，samtools index 会直接报
  # "Unsorted positions on sequence #1"，factorFootprints 就没法扫。
  run "$rbin" "$PROJECT_ROOT/src/py/atacqc.R" \
      --bam "$BAM_SRT" --outdir "$ATACQC" --sample "$SAMPLE" \
      --region "$ATACQC_REGION" --demo-region "$ATACQC_DEMO_REGION" \
      --footprint "${ATACQC_FOOTPRINT:-1}" \
      --footprint-bam "$BAM_SRT2" \
      --peaks "$PEAK/${SAMPLE}.all.filtered.narrowPeak"

  (( DRY_RUN )) && { mark_done atacqc; return 0; }

  local tsv="$ATACQC/${SAMPLE}.atacqc_metrics.tsv"
  if [[ -s "$tsv" ]]; then
    local key val
    while IFS=$'\t' read -r key val; do
      [[ "$key" == "key" || -z "$key" ]] && continue
      case "$key" in
        tss_enrichment_score_r) metric "14_atacqc" "$key" "$val" "fold" ">6" ;;
        *)                      metric "14_atacqc" "$key" "$val" "" "" ;;
      esac
    done < "$tsv"
  fi

  ok "ATACseqQC 完成（产物 $ATACQC/）"
  mark_done atacqc
}
