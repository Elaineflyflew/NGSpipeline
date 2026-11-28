#!/usr/bin/env bash
# =============================================================================
#  src/52_stage_annotate.sh —— ⑫ peak 注释 + ⑬ motif 富集
#
#  ⑫ 两套注释并行做，互相补充：
#    · HOMER annotatePeaks 给出 peak 相对基因的位置分布（promoter / intron / intergenic），
#      这是判断"数据像不像 ATAC"的直接证据 —— 好的 ATAC 有 30-50% 落在 promoter
#    · bedtools closest 给出每个 peak 最近的基因名，方便挑候选做后续验证
#
#  ⑬ HOMER findMotifsGenome 用 -mknown 只比对已知 motif 库。
#     不加 -mknown 会同时做从头 discovery，多跑几小时且结果更难解释。
#
#  两步都依赖 HOMER 的基因组数据（configureHomer.pl -install hg38）。
#  没装就 WARN 跳过，不影响主流程 —— peak 文件本身已经是有价值的产物。
# =============================================================================

stage_annotate() {
  section "⑫ peak 注释"
  checkfile "$PEAKS" "⑨ 过滤后的 peak"

  # ---- HOMER annotatePeaks ----
  # 这里不用 run：run 会把 "+ 命令" 这行也写进重定向目标，污染注释结果
  if have annotatePeaks.pl; then
    if (( DRY_RUN )); then
      printf '[dry-run] annotatePeaks.pl %s %s -annStats %s\n' \
        "$PEAKS" "$GENOME_NAME" "$ANN/${SAMPLE}.homer.annostats.txt"
    elif ! annotatePeaks.pl "$PEAKS" "$GENOME_NAME" \
          -annStats "$ANN/${SAMPLE}.homer.annostats.txt" \
          > "$ANN/${SAMPLE}.homer.anno.txt" 2> "$ANN/${SAMPLE}.homer.anno.log"; then
      warn "HOMER annotatePeaks 失败（多半是没装 $GENOME_NAME 的基因组数据），跳过"
    fi
  else
    warn "没有 annotatePeaks.pl，跳过 HOMER 注释"
  fi

  # ---- bedtools closest：最近基因 ----
  # GTF 必须先排序，否则 closest 的结果是错的（bedtools 要求 B 文件有序）
  if have bedtools; then
    local gtf_sorted="$ANN/${SAMPLE}.genes.sorted.gtf"
    runpipe "awk -F'\t' '\$3==\"gene\"' '$GTF' | sort -k1,1 -k4,4n > '$gtf_sorted'"
    runpipe "bedtools closest -D b -a '$PEAKS' -b '$gtf_sorted' > '$ANN/${SAMPLE}.bedtools.anno.txt'"
  fi

  # ---- ⑬ motif 富集 ----
  if (( ! RUN_MOTIF )); then
    log "⑬ 跳过 motif（RUN_MOTIF=0）"
    mark_done annotate
    return 0
  fi

  section "⑬ motif 富集（HOMER -mknown）"
  if ! have findMotifsGenome.pl; then
    warn "没有 findMotifsGenome.pl，跳过 motif"
    mark_done annotate
    return 0
  fi

  # top peaks：按第 9 列（-log10 q-value）排序取前 N 个
  runpipe "sort -k9,9nr '$PEAKS' | head -$TOPPEAKS_N > '$TOPPEAKS'"

  if [[ -z "$HOMER_MOTIFS" || ! -s "$HOMER_MOTIFS" ]]; then
    warn "找不到 HOMER known.motifs，跳过 motif（设 --homer-motifs 指定）"
    mark_done annotate
    return 0
  fi

  # -size given：用 peak 自身的实际宽度而不是默认的 200 bp 窗口，
  # ATAC 的峰宽差异很大，固定 200 会把宽峰截断
  if ! run findMotifsGenome.pl "$PEAKS" "$GENOME_NAME" \
        "$MOTIF/${SAMPLE}.homer.motif" \
        -mknown "$HOMER_MOTIFS" \
        -size "${MOTIF_SIZE:-given}" \
        -p "$THREADS"; then
    warn "HOMER motif 失败，跳过"
  fi

  # ---- 可选：MEME（很慢，默认关）----
  if (( RUN_MEME )) && have meme-chip; then
    section "⑬-b MEME motif"
    checkfile "$FASTA" "参考 fasta"
    if [[ -s "$MOTIF_DB" ]]; then
      run bedtools getfasta -fi "$FASTA" -bed "$PEAKS" \
          -fo "$MOTIF/${SAMPLE}.peaks.fasta"
      run meme-chip -oc "$MOTIF/${SAMPLE}.meme" \
          -db "$MOTIF_DB" \
          -p "$THREADS" \
          "$MOTIF/${SAMPLE}.peaks.fasta"
    else
      warn "找不到 JASPAR meme 库，跳过 MEME"
    fi
  fi

  mark_done annotate
}
