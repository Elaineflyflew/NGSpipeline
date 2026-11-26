#!/usr/bin/env bash
# =============================================================================
#  src/00_common.sh —— 公共基础设施
#
#  约定：本文件只定义函数，source 时不产生任何副作用（不建目录、不改 PATH）。
#  所有模块都依赖这里的 log / run / need / stage_done，所以它必须最先加载。
# =============================================================================

VERSION="2.0.0"
PROG="germline-wes-wgs"

# ---------------------------------------------------------------- 日志
log()  { printf '[%s] %s\n'      "$(date '+%m-%d %H:%M:%S')" "$*"; }
warn() { printf '[%s] WARN  %s\n' "$(date '+%m-%d %H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] ERROR %s\n' "$(date '+%m-%d %H:%M:%S')" "$*" >&2; exit 1; }

# 分节标题，让日志有层次
section() { printf '\n[%s] ===== %s =====\n' "$(date '+%m-%d %H:%M:%S')" "$*"; }

# ---------------------------------------------------------------- 路径工具
# 返回第一个存在的路径；都不存在则输出空串（不报错，交给调用方决定）
first_existing() {
  local p
  for p in "$@"; do
    [[ -e "$p" ]] && { printf '%s' "$p"; return 0; }
  done
  printf ''
}

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

# GATK bundle 必需的 4 个已知位点。.vcf 和 .vcf.gz 都算（GATK 两种都吃），
# 下载脚本、状态表、预检共用这一份清单，避免三处判断不一致。
GATK_BUNDLE_FILES="dbsnp_146.hg38 hapmap_3.3.hg38 Mills_and_1000G_gold_standard.indels.hg38 1000G_phase1.snps.high_confidence.hg38"

# GATK bundle 的已知位点：可能是 .vcf.gz（官方下载的），
# 也可能是解压后的 .vcf（GATK 自己建的 .idx）。两种都认，.gz 优先。
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

# ---------------------------------------------------------------- 工具定位
# 找可执行文件：先看 PATH，再到项目自带目录和常见的 conda/mamba 环境里找。
#
# 为什么要这一步：用户机器上往往已经装好了 gatk / vep / bcftools（可能在某个
# 没激活的 conda env 里），逼他们 conda 再装一遍既浪费空间又可能换版本。
# 找到就用，找不到才提示装。
find_tool() {
  local t="$1" p
  p=$(command -v "$t" 2>/dev/null) && { printf '%s' "$p"; return 0; }

  local home="${HOME:-/nonexistent}"
  for p in "$PROJECT_ROOT"/extratools/bin \
           "${CONDA_PREFIX:-/nonexistent}/bin" \
           /home/biotools/miniforge3/envs/*/bin \
           /opt/conda/envs/*/bin \
           "$home"/miniforge3/envs/*/bin \
           "$home"/mambaforge/envs/*/bin \
           "$home"/miniconda3/envs/*/bin \
           "$home"/anaconda3/envs/*/bin; do
    [[ -x "$p/$t" ]] && { printf '%s' "$p/$t"; return 0; }
  done
  return 1
}

# ---------------------------------------------------------------- 前置检查
# 三种模式：
#   dry-run    只记录"需要什么"，不检查
#   --check    缺了就记一笔继续走，最后统一汇报（预检的意义就是列出缺什么）
#   正常运行   缺了立刻退出
MISSING_TOOLS=0
need() {
  (( ${DRY_RUN:-0} )) && return 0

  local p
  if p=$(find_tool "$1"); then
    # 找到但不在 PATH 里（比如在别的 conda env 下）：临时加进去，别让它报错
    case ":$PATH:" in
      *":$(dirname "$p"):"*) ;;
      *) export PATH="$(dirname "$p"):$PATH" ;;
    esac
    return 0
  fi

  if (( ${CHECK_ONLY:-0} )); then
    warn "缺少命令: $1"
    MISSING_TOOLS=$(( MISSING_TOOLS + 1 ))
    return 0
  fi
  die "找不到命令 $1 —— 跑 'deploy/deploy.sh doctor' 看缺什么，再 'deploy/deploy.sh env' 补装"
}

# 必填参数：dry-run 下只提示，方便在没有数据的机器上预览命令
require_arg() {
  [[ -n "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] 缺少 %s（实际运行会报错）\n' "$2"; return 0; }
  die "$2 未指定：$3"
}

checkdir() {
  [[ -d "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] 目录不存在（实际运行会报错）: %s\n' "$1"; return 0; }
  die "$2 不存在: $1 —— 先跑 deploy/deploy.sh"
}

# ---------------------------------------------------------------- 命令执行
# 回显：只在含空格时加引号，避免 %q 把 > 和逗号转义成天书
show() {
  local a
  for a in "$@"; do
    if [[ "$a" == *" "* || -z "$a" ]]; then printf "'%s' " "$a"; else printf '%s ' "$a"; fi
  done
  printf '\n'
}

# 打印并执行；--dry-run 只打印
run() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] '; show "$@"; return 0; fi
  printf '+ '; show "$@"
  "$@"
}

# 管道用：$1 是完整命令字符串（仅用于本脚本内部拼好的命令，不接受外部输入）
runpipe() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] %s\n' "$1"; return 0; fi
  printf '+ %s\n' "$1"
  eval "$1"
}

# ---------------------------------------------------------------- 阶段状态
# STAGES 形如 "qc,trim,align"，用逗号包围做精确匹配，避免 qc 匹配到 qc2
stage_enabled() { [[ ",${STAGES:-}," == *",$1,"* ]]; }
stage_done()    { [[ -f "$STATUS/$1.done" ]] && (( ! ${FORCE:-0} )); }
# dry-run 绝不落标记，否则正式跑会被误跳过
mark_done()     { (( ${DRY_RUN:-0} )) && return 0; touch "$STATUS/$1.done"; }

# 统一调度入口：已完成则跳过
run_stage() {
  local name="$1" fn="$2"
  if stage_done "$name"; then
    log "阶段 [$name] 已完成，跳过（--force 可强制重跑）"
    return 0
  fi
  "$fn"
}

# ---------------------------------------------------------------- 其他
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
