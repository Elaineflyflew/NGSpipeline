#!/usr/bin/env bash
# =============================================================================
#  src/60_stage_report.sh —— ★ 汇总
#
#  三件事：
#    1. metrics/*.tsv → summary.tsv（按阈值打 PASS / WARN）
#    2. summary.tsv + 所有图 → 自包含 HTML 报告
#    3. MultiQC 把 fastqc / fastp / samtools stats 聚成一页（可选）
#
#  报告回答的不是"跑完了吗"，而是"跑对了吗"：每步产物在哪、指标多少、
#  和阈值比是什么结论。最后一段「三指标合看」是全篇唯一敢说
#  "这批数据能不能用"的地方。
# =============================================================================

stage_report() {
  section "★ 指标汇总 + 报告"

  run "$PYTHON" "$PROJECT_ROOT/src/py/collect_metrics.py" \
      --metrics "$METRIC" --out "$REP/summary.tsv"

  run "$PYTHON" "$PROJECT_ROOT/src/py/make_report.py" \
      --summary "$REP/summary.tsv" --sample "$SAMPLE" --outdir "$REP"

  if find_tool multiqc >/dev/null 2>&1; then
    need multiqc
    run multiqc -f -o "$REP" -n multiqc_report.html "$OUTDIR"
  else
    warn "没有 multiqc，跳过聚合页（report.html 不受影响）"
  fi

  (( DRY_RUN )) && { mark_done report; return 0; }

  local np nw
  np=$(awk -F'\t' 'NR>1 && $6=="PASS"' "$REP/summary.tsv" | wc -l | tr -d ' ')
  nw=$(awk -F'\t' 'NR>1 && $6=="WARN"' "$REP/summary.tsv" | wc -l | tr -d ' ')
  ok "指标 PASS ${np} / WARN ${nw}"
  ok "报告: $REP/report.html"
  mark_done report
}
