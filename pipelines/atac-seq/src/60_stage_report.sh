#!/usr/bin/env bash
# =============================================================================
#  src/60_stage_report.sh —— ⑮ 汇总报告 + 自动结论 + IGV 批处理
#
#  这一步是整个流程的"收口"：把前面 14 步散落在各目录里的产物和数字
#  汇成一份能直接贴进实验记录的东西，并且对每个关键指标给出 PASS/WARN/FAIL。
#
#  设计原则：
#  · 每个指标后面必须跟一句"这个数意味着什么"，光给数字没有用
#  · 判定阈值全部来自 config/pipeline.conf，改阈值不用改代码
#  · 结论分三层：单指标判定 → 分类小结 → 总体结论（能用/勉强能用/不能用）
# =============================================================================

# 从 key-value 文件里取值
kv() {
  local f="$1" k="$2"
  [[ -s "$f" ]] || { printf ''; return 0; }
  awk -F'\t' -v k="$k" '$1 == k { print $2; exit }' "$f"
}

# 判定：$1 值 $2 阈值 $3 方向（ge = 越大越好，le = 越小越好）
judge() {
  local v="$1" th="$2" dir="$3"
  [[ -z "$v" ]] && { printf 'N/A'; return 0; }
  awk -v v="$v" -v t="$th" -v d="$dir" 'BEGIN{
    ok = (d == "ge") ? (v >= t) : (v <= t)
    printf "%s", (ok ? "PASS" : "WARN")
  }'
}

# 带标记的一行：指标名 / 值 / 判定 / 说明
line3() {
  local name="$1" val="$2" mark="$3" note="$4"
  printf '  %s %s %s %s\n' \
    "$(pad "$name" 20)" "$(pad "$val" 14)" "$(pad "$mark" 6)" "$note"
}

report_metrics_tsv() {
  local out="$REP/metrics.tsv"
  {
    printf 'metric\tvalue\tthreshold\tverdict\tmeaning\n'
    printf 'alignment_rate_pct\t%s\t%s\t%s\t比对率，低说明参考版本不对或污染\n' \
      "${RATE:--}" "$MAP_RATE_MIN" "$(judge "$RATE" "$MAP_RATE_MIN" ge)"
    printf 'mito_pct\t%s\t%s\t%s\t线粒体 reads 占比，高说明文库质量差\n' \
      "${MITO:--}" "$MITO_MAX_PCT" "$(judge "$MITO" "$MITO_MAX_PCT" le)"
    printf 'nfr_pct\t%s\t40\t%s\t裸 DNA 片段占比，低说明消化不足或核完整性差\n' \
      "${NFR_PCT:--}" "$(judge "$NFR_PCT" 40 ge)"
    printf 'peaks\t%s\t%s\t%s\tpeak 数，过少说明信号弱或比对有问题\n' \
      "${PEAK_N:--}" "$PEAKS_MIN" "$(judge "$PEAK_N" "$PEAKS_MIN" ge)"
    printf 'frip\t%s\t%s\t%s\t峰内 reads 比例，ENCODE 底线 0.2\n' \
      "${FRIP:--}" "$FRIP_MIN" "$(judge "$FRIP" "$FRIP_MIN" ge)"
    printf 'tss_enrichment\t%s\t%s\t%s\tTSS 富集倍数，ENCODE 要求 >=5\n' \
      "${TSSE:--}" "$TSS_ENRICH_MIN" "$(judge "$TSSE" "$TSS_ENRICH_MIN" ge)"
    printf 'promoter_pct\t%s\t30\t%s\t峰落在启动子的比例，过低不像 ATAC\n' \
      "${PROM:--}" "$(judge "$PROM" 30 ge)"
  } > "$out"
}

report_summary() {
  local R="$REP"

  {
    echo "================================================================"
    echo " ATAC-seq 分析报告"
    echo "================================================================"
    echo
    row "sample"   "$SAMPLE"
    row "genome"   "$GENOME_LABEL ($GENOME_NAME)"
    row "index"    "$BT2_INDEX"
    row "threads"  "$THREADS"
    row "output"   "$OUTDIR"
    echo

    echo "----------------------------------------------------------------"
    echo " 一、上游处理"
    echo "----------------------------------------------------------------"
    local r1n
    r1n=$(kv "$CLEAN/${SAMPLE}.retention.txt" reads_before)
    line3 "1. FastQC"       "${r1n:-?} reads" "-"    "原始质量见 01_qc/*_fastqc.html"
    line3 "2. fastp"        "保留 $(kv "$CLEAN/${SAMPLE}.retention.txt" retention_pct)%" "-" "低于 80% 说明接头或低质量碱基多"
    line3 "3. bowtie2"      "${RATE:--}%" "$(judge "$RATE" "$MAP_RATE_MIN" ge)" "低于 ${MAP_RATE_MIN}% 查参考版本 / 污染"
    line3 "4. dedup + chrM" "${MITO:--}%" "$(judge "$MITO" "$MITO_MAX_PCT" le)" ">20% 判文库不合格，>5% 已经偏高"
    line3 "5. Tn5 shift"    "$([[ -s "$SRT2" ]] && echo done || echo skip)" "-" "正链 +4 / 负链 -5"
    echo

    echo "----------------------------------------------------------------"
    echo " 二、片段特征（⑦ / ⑭）"
    echo "----------------------------------------------------------------"
    line3 "NFR <${NFR_MAX}bp"    "${NFR_PCT:--}%"  "$(judge "$NFR_PCT" 40 ge)" "裸 DNA，开放区信号最强"
    line3 "mono 147-${MONO_MAX}" "${MONO_PCT:--}%" "-" "147-${MONO_MAX} bp"
    line3 "multi >${MONO_MAX}"   "${MULTI_PCT:--}%" "-" ">${MONO_MAX} bp，多为背景"
    line3 "fragment mode"      "${FRAG_MODE:-?} bp" "$(awk -v v="${FRAG_MODE:-0}" 'BEGIN{printf "%s", (v>=50 && v<=100) ? "PASS" : "WARN"}')" "预期 50-100 bp；右移=消化不足"
    echo "  分布图: $AQC/${SAMPLE}.fragment_length.svg"
    echo

    echo "----------------------------------------------------------------"
    echo " 三、peak calling（⑧⑨）"
    echo "----------------------------------------------------------------"
    if [[ -s "$PEAK/peak_counts.tsv" ]]; then
      printf '  %-10s %-12s %-12s %s\n' "set" "raw" "filtered" "genome_pct"
      awk -F'\t' 'NR > 1 { printf "  %-10s %-12s %-12s %s%%\n", $1, $2, $3, $4 }' \
        "$PEAK/peak_counts.tsv"
    fi
    echo "  结论: NFR 峰数应多于 mono；反过来说明片段分档或偏移有问题"
    echo

    echo "----------------------------------------------------------------"
    echo " 四、信号质量（⑩⑪）"
    echo "----------------------------------------------------------------"
    line3 "FRiP (bedtools)"  "${FRIP:--}"  "$(judge "$FRIP" "$FRIP_MIN" ge)" "ENCODE 底线 ${FRIP_MIN}，${FRIP_GOOD} 以上算好"
    line3 "FRiP (deeptools)" "${FRIP2:--}" "-" "两条算法相差应 <5%"
    line3 "TSS enrichment"   "${TSSE:--}"  "$(judge "$TSSE" "$TSS_ENRICH_MIN" ge)" "ENCODE 要求 >=${TSS_ENRICH_MIN}"
    echo "  图: $DEEP/${SAMPLE}.tss_profile.pdf  $DEEP/${SAMPLE}.fingerprint.pdf"
    echo

    echo "----------------------------------------------------------------"
    echo " 五、注释与 motif（⑫⑬）"
    echo "----------------------------------------------------------------"
    line3 "promoter pct"    "${PROM:--}%" "$(judge "$PROM" 30 ge)" "好的 ATAC 有 30-50% 落在 promoter"
    if [[ -s "$MOTIF/${SAMPLE}.homer.motif/knownResults.txt" ]]; then
      echo "  top motif:"
      awk -F'\t' 'NR > 1 && NR <= 6 { printf "    %-28s %s\n", $1, $4 }' \
        "$MOTIF/${SAMPLE}.homer.motif/knownResults.txt"
      echo "  （完整结果: $MOTIF/${SAMPLE}.homer.motif/knownResults.html）"
    else
      echo "  motif 结果未生成"
    fi
    echo
  } > "$R/summary.part1.txt"
}

# 结论段：把上面所有判定汇总成一句人话
report_conclusion() {
  local R="$REP"
  local fails=0 warns=0

  {
    echo "----------------------------------------------------------------"
    echo " 六、结论"
    echo "----------------------------------------------------------------"

    # 逐条判定，同时累计 fail / warn
    verdict_line() {
      local label="$1" mark="$2" text="$3"
      case "$mark" in
        PASS) printf '  [PASS] %s —— %s\n' "$label" "$text" ;;
        WARN) printf '  [WARN] %s —— %s\n' "$label" "$text"; warns=$(( warns + 1 )) ;;
        FAIL) printf '  [FAIL] %s —— %s\n' "$label" "$text"; fails=$(( fails + 1 )) ;;
        *)    printf '  [ -- ] %s —— %s\n' "$label" "$text" ;;
      esac
    }

    # 线粒体三档：<5 优秀 / <20 可接受 / >20 不合格
    local mmito
    mmito=$(awk -v v="${MITO:-0}" -v g="$MITO_GOOD_PCT" -v b="$MITO_MAX_PCT" \
      'BEGIN{ if (v > b) print "FAIL"; else if (v > g) print "WARN"; else print "PASS" }')
    case "$mmito" in
      PASS) verdict_line "线粒体占比" PASS "占比 ${MITO}%，低于 ${MITO_GOOD_PCT}%，文库干净" ;;
      WARN) verdict_line "线粒体占比" WARN "占比 ${MITO}%，偏高但仍在可接受范围，结论需谨慎" ;;
      FAIL) verdict_line "线粒体占比" FAIL "占比 ${MITO}%，超过 ${MITO_MAX_PCT}%，核基因组信息严重损失" ;;
    esac

    verdict_line "比对率" "$(judge "$RATE" "$MAP_RATE_MIN" ge)" "实际 ${RATE:--}%，阈值 ${MAP_RATE_MIN}%"

    local mnfr
    mnfr=$(judge "$NFR_PCT" 40 ge)
    verdict_line "NFR 比例" "$mnfr" "实际 ${NFR_PCT:--}%，40% 以上说明染色质完整性好"

    verdict_line "peak 数" "$(judge "$PEAK_N" "$PEAKS_MIN" ge)" "实际 ${PEAK_N:--}，阈值 ${PEAKS_MIN}"

    local mfrip
    mfrip=$(awk -v v="${FRIP:-0}" -v a="$FRIP_MIN" -v b="$FRIP_GOOD" \
      'BEGIN{ if (v >= b) print "PASS"; else if (v >= a) print "WARN"; else print "FAIL" }')
    case "$mfrip" in
      PASS) verdict_line "FRiP" PASS "${FRIP}，超过 ${FRIP_GOOD}，信号高度富集" ;;
      WARN) verdict_line "FRiP" WARN "${FRIP}，达到 ENCODE 底线 ${FRIP_MIN} 但不算优秀" ;;
      FAIL) verdict_line "FRiP" FAIL "${FRIP}，低于 ${FRIP_MIN}，peak 可能不可靠" ;;
    esac

    verdict_line "TSS 富集" "$(judge "$TSSE" "$TSS_ENRICH_MIN" ge)" "实际 ${TSSE:--}，阈值 ${TSS_ENRICH_MIN}"

    # 一致性：两条 FRiP 算法
    if [[ -n "$FRIP" && -n "$FRIP2" ]]; then
      local diff
      diff=$(awk -v a="$FRIP" -v b="$FRIP2" 'BEGIN{ d = a - b; if (d < 0) d = -d; printf "%.1f", d * 100 }')
      if awk -v d="$diff" 'BEGIN{ exit !(d <= 5) }'; then
        printf '  [PASS] FRiP 交叉验证 —— 两种算法相差 %s 个百分点，口径一致\n' "$diff"
      else
        printf '  [WARN] FRiP 交叉验证 —— 两种算法相差 %s 个百分点，检查 BED 坐标口径\n' "$diff"
        warns=$(( warns + 1 ))
      fi
    fi

    echo
    echo "  总体判定:"
    if (( fails == 0 && warns == 0 )); then
      echo "    数据质量合格，peak 可直接用于下游的差异可及性分析和 motif 解释。"
    elif (( fails == 0 )); then
      echo "    数据可用，但有 ${warns} 项处于临界值。peak 可以用，下结论时建议说明局限。"
    elif (( fails == 1 )); then
      echo "    有 ${fails} 项硬指标不合格、${warns} 项临界。建议先排查该指标再往下游走，"
      echo "    否则基于 peak 的生物学解释不可靠。"
    else
      echo "    有 ${fails} 项硬指标不合格。这份数据不建议直接用于生物学结论，"
      echo "    优先检查：文库起始量、Tn5 消化时间、测序是否污染。"
    fi
    echo
    echo "  下一步可以做的:"
    echo "    · 差异可及性：多样本时把 peaks 合并成 consensus 再 count reads（CSAW / DiffBind）"
    echo "    · 足迹分析：TOBIAS BINDetect 用校正后的信号扫 motif（--tobias 开启）"
    echo "    · 可视化：igv --batch $R/igv.bat"
    echo
  } > "$R/summary.part2.txt"
}

report_igv() {
  local R="$REP"
  {
    echo "# 先确认 IGV 加载的基因组版本与 BAM header 一致"
    echo "new"
    if [[ -n "$FASTA" && -s "$FASTA" ]]; then
      echo "# 用本地 fasta 作为基因组，避免 IGV 联网下载"
      echo "genome $FASTA"
    else
      echo "genome $GENOME_NAME"
    fi
    echo "load $SRT2"
    [[ -s "$BW" ]] && echo "load $BW"
    [[ -s "$PEAKS" ]] && echo "load $PEAKS"
    echo "sort base"
    echo "maxPanelHeight 700"
    echo "# $IGV_GENE 区域"
    echo "goto $IGV_REGION"
    echo "snapshot $R/igv.${IGV_GENE}.overview.png"
    echo "# 换一个 peak 密集区看细节（取 q 值最高的那个峰）"
    if [[ -s "$PEAKS" ]]; then
      local top
      top=$(sort -k9,9nr "$PEAKS" | head -1 | awk -F'\t' '{ mid = int(($2 + $3) / 2); print $1 ":" (mid - 500) "-" (mid + 500) }')
      echo "goto $top"
      echo "snapshot $R/igv.topPeak.png"
    fi
  } > "$R/igv.bat"
}

stage_report() {
  section "⑮ 生成报告"
  if (( DRY_RUN )); then
    log "（dry-run：报告依赖真实结果，跳过）"
    return 0
  fi

  # ---- 收集所有指标 ----
  RATE=$(kv "$MAP/${SAMPLE}.align_rate.txt" overall_alignment_rate)
  MITO=$(kv "$MAP/${SAMPLE}.mito.tsv" mito_pct)
  NFR_PCT=$(awk -F'\t' '$1=="nfr"{print $3}' "$BAM/${SAMPLE}.fragment_class.tsv" 2>/dev/null)
  MONO_PCT=$(awk -F'\t' '$1=="mono"{print $3}' "$BAM/${SAMPLE}.fragment_class.tsv" 2>/dev/null)
  MULTI_PCT=$(awk -F'\t' '$1=="multi"{print $3}' "$BAM/${SAMPLE}.fragment_class.tsv" 2>/dev/null)
  FRAG_MODE=$(kv "$AQC/${SAMPLE}.fragment_class_extra.tsv" fragment_mode_bp)
  PEAK_N=$(awk -F'\t' '$1=="all"{print $3}' "$PEAK/peak_counts.tsv" 2>/dev/null)
  FRIP=$(kv "$AQC/${SAMPLE}.frip.tsv" frip_bedtools)
  FRIP2=$(kv "$AQC/${SAMPLE}.frip.tsv" frip_deeptools)
  TSSE=$(kv "$AQC/${SAMPLE}.tsse.tsv" tss_enrichment)

  # promoter 比例：HOMER annostats 的第一段
  if [[ -s "$ANN/${SAMPLE}.homer.annostats.txt" ]]; then
    PROM=$(awk -F'\t' '$1 ~ /promoter-TSS/ { print $3; exit }' "$ANN/${SAMPLE}.homer.annostats.txt")
  fi

  # ---- MultiQC：把 FastQC / fastp / bowtie2 / samtools 的报告合成一份 ----
  if have multiqc; then
    ( cd "$OUTDIR" && multiqc -f -o "$REP/multiqc" 01_qc 02_clean 03_map >/dev/null 2>&1 ) \
      && log "MultiQC: $REP/multiqc/multiqc_report.html" \
      || warn "MultiQC 运行失败，忽略"
  fi

  report_summary
  report_conclusion
  report_metrics_tsv
  report_igv

  cat "$REP/summary.part1.txt" "$REP/summary.part2.txt" > "$REP/summary.txt"
  rm -f "$REP/summary.part1.txt" "$REP/summary.part2.txt"

  mark_done report
  log "报告: $REP/summary.txt"
}
