#!/usr/bin/env bash
# =============================================================================
#  src/32_stage_align.sh —— ③ 比对
#
#  两个要点：
#  1. -R 打 read group 是 GATK 的硬性要求，SM 决定最终 VCF 的样本名
#  2. 不落中间 SAM，直接管道成 BAM，省几十 G 磁盘
# =============================================================================

stage_align() {
  section "③ 比对 bwa mem"
  need bwa; need samtools

  local c1="$CLEAN/${SAMPLE}_R1.clean.fq.gz"
  local c2="$CLEAN/${SAMPLE}_R2.clean.fq.gz"
  local sam="$MAP/${SAMPLE}.q${MIN_MAPQ}.bam"
  local rg
  rg=$(printf '@RG\\tID:%s.1\\tSM:%s\\tPL:ILLUMINA\\tLB:lib_%s\\tPU:unit1' \
        "$SAMPLE" "$SAMPLE" "$SAMPLE")

  (( DRY_RUN )) || [[ -s "$c1" && -s "$c2" ]] || die "找不到修剪后的 fastq（先跑 --stages trim）"

  # 2 是 bwa 的 stderr（进度），单独存日志，别混进 BAM
  runpipe "bwa mem -t $THREADS -M -R '$rg' '$REF' '$c1' '$c2' 2> '$LOGS/${SAMPLE}.bwa.log' | samtools view -@ $THREADS -q $MIN_MAPQ -bo '$sam' -"

  (( DRY_RUN )) || samtools quickcheck -v "$sam" || die "bwa 输出 BAM 不完整: $sam"
  log "比对完成: $sam"
  mark_done align
}
