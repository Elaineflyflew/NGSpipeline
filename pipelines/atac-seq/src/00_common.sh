#!/usr/bin/env bash
# =============================================================================
#  src/00_common.sh —— 本流程特有的公共部分
#
#  通用基础设施（日志 / 命令执行 / 工具定位 / 断点状态 / 解释器选择）在
#  pipelines/_shared/common.sh，所有流程共用一份，改那里等于所有流程一起改。
#
#  这里只留"换到别的流程就不成立"的东西：
#    参考基因组和 bowtie2 索引怎么探测（ATAC 用 bowtie2，WES 用 bwa）
#    bowtie2 索引完整性检查（缺 rev 文件时 bowtie2 是 segfault，不是报错）
#    指标五列格式（metric），atac-seq 的 report 依赖它
#    contig 过滤表达式（依赖 DROP_CONTIGS，只有本流程有这个配置项）
#
#  约定：只定义函数，source 时不产生任何副作用（不建目录、不改 PATH）。
# =============================================================================

VERSION="1.1.0"
PROG="atac"

_SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../_shared" && pwd)"
# shellcheck source=../../_shared/common.sh
source "$_SHARED_DIR/common.sh"
unset _SHARED_DIR

# ---------------------------------------------------------------- 参考数据探测
# 在 reference/genome 里找基因组 fasta。不依赖固定文件名：
# UCSC 版叫 hg38.fa，Ensembl 版叫 Homo_sapiens.GRCh38.dna.primary_assembly.fa，
# 用户自己重命名成什么都行。打分规则（谁分高选谁）：
#   +2 有 .fai     +1 名字里带 hg38 / GRCh38 / assembly38
# 同分取字典序第一个，保证多次运行结果一致。
discover_genome() {
  local dir="${1:-$PROJECT_ROOT/reference/genome}"
  [[ -d "$dir" ]] || return 0

  local f score best="" best_score=-1
  for f in "$dir"/*.fa "$dir"/*.fasta "$dir"/*.fna; do
    [[ -f "$f" ]] || continue
    score=0
    [[ -s "$f.fai" ]] && score=$(( score + 2 ))
    case "$(basename "$f")" in
      *hg38*|*GRCh38*|*GRCH38*|*assembly38*) score=$(( score + 1 )) ;;
    esac
    if (( score > best_score )); then best="$f"; best_score=$score; fi
  done
  printf '%s' "$best"
}

# bowtie2 索引前缀：reference/genome/bt2/<name>.1.bt2 → 输出 .../bt2/<name>
# 找不到带 .1.bt2 的前缀就返回空。
discover_bt2() {
  local dir="${1:-$PROJECT_ROOT/reference/genome/bt2}"
  [[ -d "$dir" ]] || return 0
  local f
  for f in "$dir"/*.1.bt2 "$dir"/*.1.bt2l; do
    [[ -f "$f" ]] || continue
    f="${f%.1.bt2}"; f="${f%.1.bt2l}"
    printf '%s' "$f"
    return 0
  done
  printf ''
}

# 在目录里找第一个匹配的文件（用于 GTF / 黑名单这类"有一个就行"的参考）
discover_one() {
  local dir="$1"; shift
  local pat
  for pat in "$@"; do
    local f
    for f in "$dir"/$pat; do
      [[ -f "$f" ]] && { printf '%s' "$f"; return 0; }
    done
  done
  printf ''
}

# ---------------------------------------------------------------- 索引完整性
# bowtie2 索引是 6 个文件一组（1/2/3/4 + rev.1/rev.2）。缺 rev 的那两个时
# bowtie2 不会给出"索引不全"，而是直接 Segmentation fault —— 所以必须提前查全。
need_bt2() {
  local p="$1" ext s
  local -a miss=()
  [[ -n "$p" ]] || return 0

  if [[ -f "$p.1.bt2" ]]; then
    ext="bt2"
  elif [[ -f "$p.1.bt2l" ]]; then
    ext="bt2l"
  else
    (( ${DRY_RUN:-0} )) && { printf '[dry-run] 没有 bowtie2 索引: %s.1.bt2\n' "$p"; return 0; }
    if (( ${CHECK_ONLY:-0} )); then
      warn "没有 bowtie2 索引: $p.1.bt2"
      MISSING_TOOLS=$(( MISSING_TOOLS + 1 ))
      return 0
    fi
    die "没有 bowtie2 索引: $p.1.bt2 —— 跑 deploy/deploy.sh bt2"
  fi

  for s in 2 3 4 rev.1 rev.2; do
    [[ -f "$p.$s.$ext" ]] || miss+=("$p.$s.$ext")
  done
  (( ${#miss[@]} == 0 )) && return 0

  if (( ${DRY_RUN:-0} )); then
    printf '[dry-run] bowtie2 索引缺 %s 个: %s\n' "${#miss[@]}" "${miss[*]}"
    return 0
  fi
  if (( ${CHECK_ONLY:-0} )); then
    warn "bowtie2 索引缺 ${#miss[@]} 个: ${miss[*]}"
    MISSING_TOOLS=$(( MISSING_TOOLS + 1 ))
    return 0
  fi
  die "bowtie2 索引不完整，缺 ${#miss[@]} 个: ${miss[*]} —— bowtie2 会直接 segfault，先把这几个文件补上"
}

# ---------------------------------------------------------------- 指标落盘
# 统一五列：step / key / value / unit / benchmark
# 下游 collect_metrics.py 只认这个格式，所以写指标一律走这里。
metric() {
  local step="$1" key="$2" value="$3" unit="${4:-}" bench="${5:-}"
  (( ${DRY_RUN:-0} )) && return 0
  mkdir -p "$METRIC"
  local f="$METRIC/${step}.tsv"
  # 同一个 step+key 重跑时覆盖而不是追加。纯追加的话，断点续跑几次之后
  # summary.tsv 里同一行会出现两遍（实测 04_bam_filter / 07_fragbin 就重复过）。
  if [[ -s "$f" ]]; then
    awk -F'\t' -v s="$step" -v k="$key" '!($1==s && $2==k)' "$f" > "$f.tmp" \
      && mv "$f.tmp" "$f"
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$step" "$key" "$value" "$unit" "$bench" >> "$f"
}

# ---------------------------------------------------------------- 本流程专用
# 把 DROP_CONTIGS（含 * 通配）翻成 samtools -e 表达式。
# 例：("chrM" "chrUn*" "*_random") → rname != "chrM" && rname !~ "chrUn.*" && rname !~ ".*_random"
contig_expr() {
  local c parts=()
  for c in ${DROP_CONTIGS[@]+"${DROP_CONTIGS[@]}"}; do
    [[ -z "$c" ]] && continue
    if [[ "$c" == *"*"* ]]; then
      parts+=("rname !~ \"${c//\*/.*}\"")
    else
      parts+=("rname != \"$c\"")
    fi
  done
  if (( ${#parts[@]} == 0 )); then printf 'true'; return 0; fi

  # 注意不能用 "${parts[*]}" 加 IFS=' && '：[*] 只认 IFS 的第一个字符，
  # 分隔符会退化成空格，拼出来的表达式语法就错了
  local out="" p
  for p in "${parts[@]}"; do
    out="${out:+$out && }$p"
  done
  printf '%s' "$out"
}

# TSS BED 是所有样本共用的，做一次就够；存在且非空就不重做
tss_bed_ready() { [[ -s "$TSS_BED" ]]; }
