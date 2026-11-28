#!/usr/bin/env bash
# =============================================================================
#  src/00_common.sh —— 本流程特有的公共部分
#
#  通用基础设施（日志 / 命令执行 / 工具定位 / 断点状态 / 解释器选择）在
#  pipelines/_shared/common.sh，所有流程共用一份，改那里等于所有流程一起改。
#
#  这里只留"换到别的流程就不成立"的东西：
#    参考基因组怎么探测（变异检测要 .bwt/.dict，ATAC 要 bowtie2 索引）
#    GATK bundle 的命名规则
#    BAM 有没有 read group（GATK 的硬要求）
#    参考基因组的 contig 命名风格（chr17 vs 17，用来查 VCF/BED 是否对得上）
#
#  约定：只定义函数，source 时不产生任何副作用（不建目录、不改 PATH）。
# =============================================================================

VERSION="2.1.0"
PROG="germline-wes-wgs"

_SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../_shared" && pwd)"
# shellcheck source=../../_shared/common.sh
source "$_SHARED_DIR/common.sh"
unset _SHARED_DIR

# ---------------------------------------------------------------- 参考基因组
# 在目录里找参考基因组。不依赖固定文件名，也不依赖软链接：
#   Broad 版叫 Homo_sapiens_assembly38.fasta，Ensembl 版叫
#   Homo_sapiens.GRCh38.dna.primary_assembly.fa，用户自己重命名成什么都行。
# 打分规则（谁分高选谁）：
#   +2 有 .fai        +2 有 .bwt        +1 有 .dict
#   +1 名字里带 assembly38 / GRCh38 / hg38
# 同分取字典序第一个，保证多次运行结果一致。
discover_ref() {
  local dir="${1:-$PROJECT_ROOT/reference/genome}"
  [[ -d "$dir" ]] || return 0

  local f score best="" best_score=-1
  for f in "$dir"/*.fa "$dir"/*.fasta "$dir"/*.fna; do
    [[ -f "$f" ]] || continue
    score=0
    [[ -s "$f.fai" ]] && score=$(( score + 2 ))
    [[ -s "$f.bwt" ]] && score=$(( score + 2 ))
    [[ -s "${f%.fa}.dict" || -s "${f%.fasta}.dict" || -s "${f%.fna}.dict" ]] && score=$(( score + 1 ))
    case "$(basename "$f")" in
      *assembly38*|*GRCh38*|*hg38*|*GRCH38*) score=$(( score + 1 )) ;;
    esac
    if (( score > best_score )); then best="$f"; best_score=$score; fi
  done
  printf '%s' "$best"
}

# ---------------------------------------------------------------- GATK bundle
# GATK bundle 必需的 4 个已知位点。.vcf 和 .vcf.gz 都算（GATK 两种都吃），
# 下载脚本、状态表、预检共用这一份清单，避免三处判断不一致。
GATK_BUNDLE_FILES="dbsnp_146.hg38 hapmap_3.3.hg38 Mills_and_1000G_gold_standard.indels.hg38 1000G_phase1.snps.high_confidence.hg38"

# 已知位点：可能是 .vcf.gz（官方下载的），也可能是解压后的 .vcf
# （GATK 自己建的 .idx）。两种都认，.gz 优先。
bundle_file() {
  local n="$1" f
  [[ -n "${BUNDLE:-}" ]] || return 1
  for f in "$BUNDLE/$n.vcf.gz" "$BUNDLE/$n.vcf"; do
    [[ -s "$f" ]] && { printf '%s' "$f"; return 0; }
  done
  return 1
}

# 未压缩 VCF 要 .idx（GATK 建），压缩的要 .tbi/.csi。缺了就补，否则
# BQSR 会报 "An index is required"。
ensure_vcf_index() {
  local f="$1"
  [[ -s "$f" ]] || return 0
  (( ${DRY_RUN:-0} )) && return 0
  case "$f" in
    *.vcf.gz) [[ -s "$f.tbi" || -s "$f.csi" ]] && return 0 ;;
    *.vcf)    [[ -s "$f.idx" ]] && return 0 ;;
    *)        return 0 ;;
  esac
  warn "缺少索引，补建: $(basename "$f")"
  gatk IndexFeatureFile -I "$f" >/dev/null 2>&1 \
    || warn "补建失败（不影响流程继续，但 GATK 可能会报错）: $f"
}

# ---------------------------------------------------------------- 本流程专用
# 判断 BAM 是否带 @RG。不能写成 'samtools view -H | grep -q'：
# pipefail 下 grep 提前退出会让 samtools 收到 SIGPIPE，判断结果反掉
bam_has_rg() {
  local hdr="${TMPD}/.header.$$.sam"
  samtools view -H "$1" > "$hdr" 2>/dev/null || return 1
  local rc=0
  grep -q '^@RG' "$hdr" || rc=1
  rm -f "$hdr"
  return $rc
}

# 参考基因组的 contig 命名风格（chr17 vs 17），用于一致性检查
ref_style() {
  local f="${REF}.fai"
  [[ -s "$f" ]] || { printf 'unknown'; return 0; }
  if head -1 "$f" | grep -q '^chr'; then printf 'ucsc'; else printf 'ensembl'; fi
}
