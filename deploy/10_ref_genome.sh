#!/usr/bin/env bash
# =============================================================================
#  deploy/10_ref_genome.sh —— 参考基因组 + 索引
#
#  两种来源，用 --source 选：
#    broad   （默认，推荐）GATK 官方 hg38 组装，chr 前缀，自带 bwa 索引和 .dict，
#            和 gatk_bundle 里的 dbsnp/hapmap 命名完全一致，下载即用。
#    ensembl Ensembl 主组装，无 chr 前缀，需要本地建索引（bwa index 约 40 分钟）；
#            加 --rename-chr 可转成 chr 前缀以匹配 GATK bundle。
#
#  用法
#    bash deploy/10_ref_genome.sh                 # 默认 broad
#    bash deploy/10_ref_genome.sh --source ensembl --rename-chr
# =============================================================================
set -euo pipefail

BROAD_BASE="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-112}"

deploy_ref() {
  local source="broad" rename_chr=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --source)     source="$2"; shift 2 ;;
      --rename-chr) rename_chr=1; shift ;;
      -h|--help)
        echo "用法: $0 [--source broad|ensembl] [--rename-chr]"
        return 0 ;;
      *) die "未知参数: $1" ;;
    esac
  done

  mkdir -p "$GENOME_DIR"
  check_space 30 "$GENOME_DIR"

  case "$source" in
    broad)   _ref_broad ;;
    ensembl) _ref_ensembl "$rename_chr" ;;
    *) die "--source 只能是 broad 或 ensembl" ;;
  esac

  # dry-run 下文件当然不存在，跳过校验
  (( DRY_RUN )) && { log "（dry-run）参考基因组将放在 $GENOME_DIR"; return 0; }

  local ref
  ref=$(find_ref) || die "参考基因组下载后仍找不到，检查 $GENOME_DIR"
  log "参考基因组就绪: $ref"
  log "pipeline.conf 会自动探测到这个文件，不需要改配置"
}

_ref_broad() {
  log "来源: GATK hg38（Broad）—— 自带索引，命名 chr*，与 gatk_bundle 一致"
  # 主 fasta 约 3G，索引约 5G
  local f
  for f in Homo_sapiens_assembly38.fasta \
           Homo_sapiens_assembly38.fasta.fai \
           Homo_sapiens_assembly38.dict \
           Homo_sapiens_assembly38.fasta.amb \
           Homo_sapiens_assembly38.fasta.ann \
           Homo_sapiens_assembly38.fasta.bwt \
           Homo_sapiens_assembly38.fasta.pac \
           Homo_sapiens_assembly38.fasta.sa \
           Homo_sapiens_assembly38.fasta.alt; do
    dl "$BROAD_BASE/$f" "$GENOME_DIR/$f"
  done

  # 不建软链接：pipeline 按 .fai/.bwt/.dict 自动识别，文件名叫什么都能找到。
  # 软链接反而有坑——rsync 不带 -l 时会被展开成实体文件，或者指向失效。
}

_ref_ensembl() {
  local rename_chr="$1"
  local base="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna"
  local gz="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
  local fa="${gz%.gz}"

  log "来源: Ensembl release ${ENSEMBL_RELEASE}"
  dl "$base/$gz" "$GENOME_DIR/$gz"
  if (( ! DRY_RUN )) && [[ ! -s "$GENOME_DIR/$fa" ]]; then
    log "解压（约 3G）"
    gzip -dc "$GENOME_DIR/$gz" > "$GENOME_DIR/$fa"
  fi

  local ref="$GENOME_DIR/$fa"
  (( DRY_RUN )) && { printf '[dry-run] 建 faidx / bwa index / dict\n'; return 0; }

  # 改名要在建索引之前做，否则索引里的名字和 fasta 不一致
  if (( rename_chr )); then
    log "把 contig 名改成 chr 前缀（对齐 GATK bundle）"
    local map="$GENOME_DIR/rename_chrs.tsv"
    awk '{print $1"\tchr"$1}' "$ref.fai" 2>/dev/null > "$map" || true
    if [[ ! -s "$map" ]]; then
      # 还没建 fai，先直接从 fasta 头取
      grep '^>' "$ref" | sed 's/^>//; s/ .*//' | awk '{print $1"\tchr"$1}' > "$map"
    fi
    awk -v mapf="$map" '
      BEGIN { while ((getline l < mapf) > 0) { split(l, a, "\t"); m[a[1]] = a[2] } }
      /^>/   { h = substr($0, 2); split(h, b, " "); if (b[1] in m) sub(">" b[1], ">" m[b[1]]) }
             { print }
    ' "$ref" > "${ref}.renamed"
    mv "${ref}.renamed" "$ref"
    log "改名映射表: $map"
  fi

  need samtools
  log "samtools faidx"
  samtools faidx "$ref"

  if [[ ! -s "$ref.bwt" ]] || (( FORCE )); then
    need bwa
    log "bwa index（单线程，约 40 分钟，请耐心）"
    bwa index "$ref"
  fi

  local dict="${ref%.fa}.dict"
  if [[ ! -s "$dict" ]] || (( FORCE )); then
    if command -v gatk >/dev/null 2>&1; then
      log "gatk CreateSequenceDictionary"
      gatk CreateSequenceDictionary -R "$ref" -O "$dict"
    else
      log "samtools dict（没装 gatk 时的等价替代）"
      samtools dict "$ref" -o "$dict"
    fi
  fi
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
  deploy_ref "$@"
fi
