#!/usr/bin/env bash
# =============================================================================
#  src/00_common.sh —— 公共基础设施
#
#  约定：本文件只定义函数，source 时不产生任何副作用（不建目录、不改 PATH）。
#  所有模块都依赖这里的 log / run / need / stage_done，所以它必须最先加载。
# =============================================================================

VERSION="1.0.0"
PROG="atac"

# ---------------------------------------------------------------- 日志
log()  { printf '[%s] %s\n'       "$(date '+%m-%d %H:%M:%S')" "$*"; }
ok()   { printf '[%s]   ✔ %s\n'   "$(date '+%m-%d %H:%M:%S')" "$*"; }
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

# ---------------------------------------------------------------- 工具定位
# 找可执行文件：先看 PATH，再到项目自带目录和常见的 conda/mamba 环境里找。
#
# 为什么要这一步：用户机器上往往已经装好了 samtools / macs2（可能在某个
# 没激活的 conda env 里），逼他们再装一遍既浪费空间又可能换版本。
find_tool() {
  local t="$1" p
  p=$(command -v "$t" 2>/dev/null) && { printf '%s' "$p"; return 0; }

  local home="${HOME:-/nonexistent}"
  for p in "$PROJECT_ROOT"/extratools/bin \
           "${CONDA_PREFIX:-/nonexistent}/bin" \
           /home/biotools/miniforge3/envs/*/bin \
           /bioapps/bin \
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
#   --check    缺了就记一笔继续走，最后统一汇报
#   正常运行   缺了立刻退出
MISSING_TOOLS=0
# dry-run 下缺的参考数据另记一笔：不算失败，但汇总时不能说"齐全"
SOFT_MISSING=0
need() {
  (( ${DRY_RUN:-0} )) && return 0

  local p
  if p=$(find_tool "$1"); then
    # 找到但不在 PATH 里（比如在别的 conda env 下）：临时加进去
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

# 必填文件：dry-run 下只提示，方便在没有数据的机器上预览命令
need_file() {
  [[ -s "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] 文件不存在（实际运行会报错）: %s\n' "$1"; return 0; }
  die "$2 不存在: $1 —— 先跑 deploy/deploy.sh"
}

# 必填变量（比如 BT2 这种由自动探测填的路径）
need_var() {
  [[ -n "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] %s\n' "$2"; return 0; }
  die "$2"
}

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

# ---------------------------------------------------------------- 工具兼容
# FastQC 0.12 起把 --noextract 删掉了（默认就不解压，要解压才写 --extract），
# 0.11 及以前反过来。探一次帮助文本，两种版本都能跑。
FASTQC_EXTRA=()
fastqc_probe() {
  FASTQC_EXTRA=()
  command -v fastqc >/dev/null 2>&1 || return 0
  if fastqc -h 2>&1 | grep -q -- '--noextract'; then
    FASTQC_EXTRA=(--noextract)
  fi
  return 0
}

# 选一个 python3 解释器。优先挑 import 得动 matplotlib 的那个——⑦ 的片段分布图
# 和 ⑫ 的注释饼图要用它；实在没有就退回普通 python3，图跳过、数值不受影响。
pick_python() {
  local c
  for c in ${PYTHON:-} "$(find_tool python3 || true)" /usr/bin/python3; do
    [[ -n "$c" ]] || continue
    [[ -x "$c" ]] || continue
    if "$c" -c 'import matplotlib' >/dev/null 2>&1; then
      printf '%s' "$c"
      return 0
    fi
  done
  find_tool python3 || printf 'python3'
}

# ---------------------------------------------------------------- 挑 R
# 为什么需要这个：课程服务器上同时装着两套 R —— 系统 R（/usr/bin/Rscript，
# 4.6，Bioconductor 3.24）和 conda R（/bioapps/bin/Rscript，4.4，Bioconductor
# 3.20）。两套的 libPaths 不同，同一个包经常只在其中一套里有。而 /bioapps/bin
# 排在 PATH 最前面，所以"直接 Rscript"实际用的是 conda R，可能恰好缺你要的包。
# 做法：把所有能解析到的 Rscript 都试一遍，选第一个"这些包全都装得上"的。
# 全都不满足时输出空串并返回 1，由调用方决定是降级还是报错。
_r_has_pkgs() {
  local rs="$1"; shift
  local list="" p
  for p in "$@"; do list+="\"$p\","; done
  "$rs" -e "for (p in c(${list%,})) if (!requireNamespace(p, quietly=TRUE)) quit(status=1)" \
      >/dev/null 2>&1
}

pick_rscript() {
  local cands=() c seen=" "
  while IFS= read -r c; do
    [[ -n "$c" && -x "$c" ]] || continue
    [[ "$seen" == *" $c "* ]] && continue
    seen+="$c "
    cands+=("$c")
  done < <(type -a -p Rscript 2>/dev/null || true; printf '%s\n' /usr/bin/Rscript /bioapps/bin/Rscript)
  for c in "${cands[@]}"; do
    if _r_has_pkgs "$c" "$@"; then printf '%s' "$c"; return 0; fi
  done
  return 1
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
# 失败时不靠 ERR trap 兜（那样只能报出 "run() 第 159 行: $@"，看不出是哪条命令），
# 在这里把整条命令原样回显出来再退出。
run() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] '; show "$@"; return 0; fi
  printf '+ '; show "$@"
  local rc=0
  "$@" || rc=$?
  (( rc == 0 )) || fail_cmd "$rc" "$@"
}

# 管道用：$1 是完整命令字符串（只用于本脚本内部拼好的命令，不接受外部输入）
runpipe() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] %s\n' "$1"; return 0; fi
  printf '+ %s\n' "$1"
  local rc=0
  eval "$1" || rc=$?
  (( rc == 0 )) || fail_cmd "$rc" "$1"
}

# 统一的外部命令失败出口
fail_cmd() {
  local rc="$1"; shift
  printf '\n[FAIL] 命令退出码 %s: ' "$rc" >&2
  show "$@" >&2
  die "上一条命令失败，日志: ${LOGS:-./results/00_logs}"
}

# ---------------------------------------------------------------- 非致命执行
# 给"锦上添花"的步骤用（IGV 截图、可选的交叉验证 QC）：这些步骤失败不该让
# 主链白跑。语义是"永远返回 0"，调用方自己按产物判断成没成 —— 因为很多图形
# 程序退出码是 0 但根本没干活，光看 rc 判断不了。
try_run() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] '; show "$@"; return 0; fi
  printf '+ '; show "$@"
  local rc=0
  "$@" || rc=$?
  (( rc == 0 )) || warn "命令退出码 $rc（不中断主链）: $(show "$@")"
  return 0
}

try_runpipe() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] %s\n' "$1"; return 0; fi
  printf '+ %s\n' "$1"
  local rc=0
  eval "$1" || rc=$?
  (( rc == 0 )) || warn "命令退出码 $rc（不中断主链）: $1"
  return 0
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

# 从 samtools stats 的 SN 行取值：stats_field <file> <key>
# SN 行格式：SN\traw total sequences:\t20000000\t# excluding...
stats_field() {
  awk -F'\t' -v k="$2" '$1=="SN" && index($2, k":")==1 {print $3; exit}' "$1"
}

# ---------------------------------------------------------------- 其他
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

# FastQC 的输出名由输入文件名推出：x_1.fastq.gz → x_1_fastqc.html
# 这里给出同一个 stem，用于在 --noextract 模式下回读 zip 里的 fastqc_data.txt
fastqc_stem() {
  local b
  b=$(basename "$1")
  b="${b%.gz}"
  b="${b%.fastq}"; b="${b%.fq}"
  printf '%s' "$b"
}

# 把任意标签转成安全的指标键名：Promoter (<=1kb) → promoter_1kb
slug() {
  printf '%s' "$1" \
    | tr '[:upper:]' '[:lower:]' \
    | sed -E 's/[^a-z0-9]+/_/g; s/^_+//; s/_+$//'
}

# 从 FastQC zip 里读一个字段：fastqc_field <zip> <字段名>
fastqc_field() {
  local z="$1" k="$2"
  [[ -s "$z" ]] || return 1
  command -v unzip >/dev/null 2>&1 || return 1
  unzip -p "$z" '*/fastqc_data.txt' 2>/dev/null \
    | awk -F'\t' -v k="$k" '$1==k {print $2; exit}'
}

# TSS BED 是所有样本共用的，做一次就够；存在且非空就不重做
tss_bed_ready() { [[ -s "$TSS_BED" ]]; }
