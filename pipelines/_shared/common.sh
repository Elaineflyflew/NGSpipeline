#!/usr/bin/env bash
# =============================================================================
#  pipelines/_shared/common.sh —— 所有流程共用的基础设施
#
#  这里只放"与数据类型无关"的东西：日志、命令执行、工具定位、断点状态。
#
#  收录判据只有一条：**把这个函数原样搬到另一个流程，一行都不用改。**
#  只要它开始关心"参考基因组长什么样"或者"这一步要产出什么"，就说明它属于
#  某个流程的 src/00_common.sh，不属于这里。
#
#  因此下面这些**不**在共用层：
#    discover_*       各流程探测的参考数据结构不同（bwa 索引 vs bowtie2 索引）
#    stage_*          业务逻辑，一个流程一套
#    bundle_file()    GATK bundle 的命名规则，只有变异检测流程关心
#    metric()         指标五列格式是 atac-seq 定的，别的流程未必照抄
#    contig_expr()    依赖 DROP_CONTIGS 这个 atac-seq 才有的配置项
#
#  约定：只定义函数和少量默认变量，source 时不产生任何副作用
#  （不建目录、不改 PATH、不打印任何东西）。
#  所有流程的 src/00_common.sh 都 source 它，且必须早于其它模块。
# =============================================================================

# ---------------------------------------------------------------- 日志
# ok() 前面留的空格是为了和 log() 的 "[11-28 22:40:18] " 对齐，
# 让 ✔ / WARN / ERROR 三种行在终端里竖直对齐
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

# ---------------------------------------------------------------- 工具定位
# 找可执行文件：先看 PATH，再到项目自带目录和常见的 conda/mamba 环境里找。
#
# 为什么要这一步：用户机器上往往已经装好了 gatk / samtools / macs2（可能在某个
# 没激活的 conda env 里），逼他们再装一遍既浪费空间又可能换版本。找到就用，
# 找不到才由 need() 决定是报错还是记一笔。
#
# 这里**故意不写死课程镜像的 /bioapps/bin**：那是"这个流程信任哪些目录"的策略，
# 属于各流程 config 里的 EXTRA_PATH（20_context.sh 会把它加进 PATH）。写进共用层
# 有两个坏处：一是和 config 重复，二是让测试无法用桩目录模拟"某个工具没装"——
# 删掉桩里的 gatk，真 gatk 还是会从 /bioapps/bin 被扫出来，7 项 deploy 测试全废。
#
# 想加自己的搜索目录，设 TOOL_SEARCH_DIRS（空格分隔）。
find_tool() {
  local t="$1" p
  p=$(command -v "$t" 2>/dev/null) && { printf '%s' "$p"; return 0; }

  # 测试用：只看 PATH，不去翻 conda 环境。冒烟测试要在桩目录里模拟"某个工具
  # 没装"，而服务器上 /home/biotools/.../envs/*/bin 里有真 gatk —— 不关掉这一步，
  # 删掉桩里的 gatk 也照样能扫到，测试就永远测不出东西。
  (( ${TOOL_SEARCH_DISABLE:-0} )) && return 1

  local home="${HOME:-/nonexistent}"
  local -a dirs=()
  # shellcheck disable=SC2206  # 这里就是要按空格拆成多个目录
  [[ -n "${TOOL_SEARCH_DIRS:-}" ]] && dirs+=($TOOL_SEARCH_DIRS)
  dirs+=("$PROJECT_ROOT"/extratools/bin
         "${CONDA_PREFIX:-/nonexistent}/bin"
         /home/biotools/miniforge3/envs/*/bin
         /opt/conda/envs/*/bin
         "$home"/miniforge3/envs/*/bin
         "$home"/mambaforge/envs/*/bin
         "$home"/miniconda3/envs/*/bin
         "$home"/anaconda3/envs/*/bin)

  for p in "${dirs[@]}"; do
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
# dry-run 下缺的**参考数据**另记一笔：dry-run 本来就不该因为没数据而失败，
# 但预检汇总时不能说"齐全"。用这个计数单独区分"缺工具"和"缺数据"。
SOFT_MISSING=0
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

# 必填文件：dry-run 下只提示，方便在没有数据的机器上预览命令
need_file() {
  [[ -s "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] 文件不存在（实际运行会报错）: %s\n' "$1"; return 0; }
  die "$2 不存在: $1 —— 先跑 deploy/deploy.sh"
}

# 必填变量（比如自动探测填进来的路径）
need_var() {
  [[ -n "$1" ]] && return 0
  (( ${DRY_RUN:-0} )) && { printf '[dry-run] %s\n' "$2"; return 0; }
  die "$2"
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
# 失败时不靠 ERR trap 兜（那样只能报出 "run() 第 159 行: $@"，看不出是哪条命令），
# 在这里把整条命令原样回显出来再退出。
run() {
  if (( ${DRY_RUN:-0} )); then printf '[dry-run] '; show "$@"; return 0; fi
  printf '+ '; show "$@"
  local rc=0
  "$@" || rc=$?
  (( rc == 0 )) || fail_cmd "$rc" "$@"
}

# 管道用：$1 是完整命令字符串（仅用于本脚本内部拼好的命令，不接受外部输入）
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
# 给"锦上添花"的步骤用（截图、可选的交叉验证 QC）：这些步骤失败不该让主链
# 白跑。语义是"永远返回 0"，调用方自己按产物判断成没成 —— 因为很多图形程序
# 退出码是 0 但根本没干活，光看 rc 判断不了。
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

# 选一个 python3 解释器。优先挑 import 得动 matplotlib 的那个——画图步骤要用它；
# 实在没有就退回普通 python3，图跳过、数值不受影响。
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
# 为什么需要这个：课程服务器上同时装着两套 R —— 系统 R（/usr/bin/Rscript）和
# conda R（/bioapps/bin/Rscript）。两套的 libPaths 不同，同一个 Bioconductor 包
# 经常只在其中一套里有。而 /bioapps/bin 排在 PATH 最前面，所以"直接 Rscript"
# 实际用的是 conda R，可能恰好缺你要的包。
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

# ---------------------------------------------------------------- 小工具
# 把任意标签转成安全的指标键名：Promoter (<=1kb) → promoter_1kb
slug() {
  printf '%s' "$1" \
    | tr '[:upper:]' '[:lower:]' \
    | sed -E 's/[^a-z0-9]+/_/g; s/^_+//; s/_+$//'
}

# 从 samtools stats 的 SN 行取值：stats_field <file> <key>
# SN 行格式：SN\traw total sequences:\t20000000\t# excluding...
stats_field() {
  awk -F'\t' -v k="$2" '$1=="SN" && index($2, k":")==1 {print $3; exit}' "$1"
}

# FastQC 的输出名由输入文件名推出：x_1.fastq.gz → x_1_fastqc.html
# 给出同一个 stem，用于在 --noextract 模式下回读 zip 里的 fastqc_data.txt
fastqc_stem() {
  local b
  b=$(basename "$1")
  b="${b%.gz}"
  b="${b%.fastq}"; b="${b%.fq}"
  printf '%s' "$b"
}

# 从 FastQC zip 里读一个字段：fastqc_field <zip> <字段名>
fastqc_field() {
  local z="$1" k="$2"
  [[ -s "$z" ]] || return 1
  command -v unzip >/dev/null 2>&1 || return 1
  unzip -p "$z" '*/fastqc_data.txt' 2>/dev/null \
    | awk -F'\t' -v k="$k" '$1==k {print $2; exit}'
}
