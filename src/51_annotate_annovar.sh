#!/usr/bin/env bash
# =============================================================================
#  src/51_annotate_annovar.sh —— ANNOVAR 注释
#
#  protocol/operation 必须一一对应：g=gene, r=region, f=filter。
#  改 --annovar-protocols 时记得同步 --annovar-ops，数量不一致 table_annovar 会报错。
# =============================================================================

annotate_annovar() {
  local vcf="$1" tag="$2"
  need table_annovar.pl
  require_arg "$ANNOVAR_DB" "--annovar-db" "ANNOVAR humandb 目录"
  checkdir "$ANNOVAR_DB" "ANNOVAR humandb"

  local np no
  np=$(awk -F, '{print NF}' <<< "$ANNOVAR_PROTOCOLS")
  no=$(awk -F, '{print NF}' <<< "$ANNOVAR_OPS")
  [[ "$np" == "$no" ]] || die "--annovar-protocols 有 ${np} 项，--annovar-ops 有 ${no} 项，必须一致"

  run table_annovar.pl "$vcf" "$ANNOVAR_DB" \
      -buildver hg38 -out "$ANN/${SAMPLE}.${tag}.annovar" -remove \
      -protocol "$ANNOVAR_PROTOCOLS" -operation "$ANNOVAR_OPS" \
      -nastring . -vcfinput --thread "$THREADS"

  log "ANNOVAR 完成: $ANN/${SAMPLE}.${tag}.annovar.hg38_multianno.txt"
}
