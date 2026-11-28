#!/usr/bin/env bash
# =============================================================================
#  src/33_stage_bam.sh —— ④ 质量过滤 + 去线粒体 → sort → index
#
#  两个坑：
#  1. 不能直接 `grep -v chrM` 过滤 SAM 文本 —— header 里的
#     `@SQ SN:chrM LN:16569` 也含 chrM，会被一起删掉，BAM 就丢了这条染色体的
#     定义。必须只对数据行判断第 3 列（RNAME）。
#  2. 用 -q 20 而不是 -q 30：ATAC 的 read 大量落在重复序列和低复杂度区，
#     阈值太高会连真正的开放区一起削掉。20 是 ENCODE 的常用值。
#
#  线粒体占比是 ATAC-seq 第一个硬指标：文库质量差时 chrM 能占到 50% 以上，
#  这些 read 不含任何核基因组信息，必须去掉再算下游指标。
# =============================================================================

stage_bam() {
  section "④ 过滤 + 去线粒体 + 排序"
  need samtools

  if [[ -n "$BAM_IN" ]]; then
    warn "已提供 --bam，跳过 ④（直接使用 ${BAM_IN}）"
    mark_done bam
    return 0
  fi

  checkfile "$SAM" "bowtie2 输出的 SAM"

  # awk 一次扫描同时做三件事：保留 header、丢线粒体、统计占比
  local awkprog
  awkprog='
    BEGIN { n = split(mt, a, " "); for (i = 1; i <= n; i++) m[a[i]] = 1 }
    /^@/ { print; next }
    {
      tot++
      if ($3 in m) { mito++; if (drop) next }
      print
    }
    END { printf "q20_reads\t%d\nmito_reads\t%d\n", tot, mito > stat }
  '

  runpipe "samtools view -@ $THREADS -h -q $MIN_MAPQ '$SAM' \
      | awk -v mt='$MT_NAMES' -v drop=$DROP_MT -v stat='$MAP/${SAMPLE}.mito.tsv' '$awkprog' \
      | samtools view -@ $THREADS -b -o '$Q20' -"

  run samtools sort -@ "$THREADS" -m 1G -T "$TMPD/sort.$SAMPLE" -o "$SRT" "$Q20"
  run samtools index -@ "$THREADS" "$SRT"

  if (( ! DRY_RUN )); then
    samtools flagstat -@ "$THREADS" "$SRT" > "$MAP/${SAMPLE}.flagstat.txt"
    local tot mito
    tot=$(awk -F'\t' '$1=="q20_reads"{print $2}' "$MAP/${SAMPLE}.mito.tsv")
    mito=$(awk -F'\t' '$1=="mito_reads"{print $2}' "$MAP/${SAMPLE}.mito.tsv")
    printf 'mito_pct\t%s\n' "$(pct "${mito:-0}" "${tot:-0}")" >> "$MAP/${SAMPLE}.mito.tsv"
    log "MAPQ>=${MIN_MAPQ}: ${tot:-0} 条，其中线粒体 ${mito:-0} 条（$(pct "${mito:-0}" "${tot:-0}")%）"
  fi

  (( KEEP_TMP )) || { (( DRY_RUN )) || rm -f "$Q20"; }
  mark_done bam
}
