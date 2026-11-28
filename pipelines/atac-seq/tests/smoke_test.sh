#!/usr/bin/env bash
# =============================================================================
#  tests/smoke_test.sh —— 冒烟测试
#
#  不需要任何生信工具、不需要参考数据：全部用 --dry-run 验证脚手架。
#  真数据必须上服务器跑，这里只保证"参数解析、路径派生、阶段跳过、命令拼装"
#  这些不需要数据就能验证的部分是对的。
#
#  用法: bash tests/smoke_test.sh
# =============================================================================
set -uo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PIPE="bin/atac-seq"
PASS=0; FAIL=0
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

ok()   { PASS=$((PASS+1)); printf '  ok   %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '  FAIL %s\n' "$1"; }
check(){ if eval "$2" >/dev/null 2>&1; then ok "$1"; else bad "$1"; fi; }

R1="$TMP/Case1_R1.fastq.gz"; R2="$TMP/Case1_R2.fastq.gz"
: > "$R1"; : > "$R2"
OUT="$TMP/out"

# dry-run 的输出先落文件再 grep。
# 不能直接 `cmd | grep -q`：pipefail 下 grep 命中就退出，上游写 stdout 时
# 收到 SIGPIPE，整条管道返回非 0，判断会反掉。
dry_to() {
  local log="$1"; shift
  "$PIPE" "$@" --dry-run > "$log" 2>&1 || true
  tr '\n' ' ' < "$log"
}

# ---------------------------------------------------------------- 0. 语法
echo "[0] 语法与编码"
for f in bin/atac-seq src/*.sh deploy/*.sh tests/*.sh; do
  if bash -n "$f" 2>/dev/null; then ok "bash -n $f"; else bad "bash -n $f"; fi
done

badpunct=$(LC_ALL=C grep -nE '\$[A-Za-z_][A-Za-z0-9_]*[^ -~]' src/*.sh bin/atac-seq deploy/*.sh 2>/dev/null || true)
if [[ -z "$badpunct" ]]; then ok "没有 \${VAR} 紧贴全角标点"; else bad "变量名紧贴全角标点：$badpunct"; fi

# ---------------------------------------------------------------- 1. 帮助
echo "[1] 命令行接口"
check "--help 能退出 0"        "$PIPE --help"
check "--help 有用法段"        "$PIPE --help | grep -q '用法'"
check "--version"             "$PIPE --version | grep -q 'atac-seq'"
NSTAGE=$($PIPE --list-stages 2>/dev/null | grep -cE '^  (qc|trim|align|bam|shift|fraction|peak|metrics|bigwig|annotate|atacqc|report) ')
check "--list-stages 列出 12 个阶段" "[[ $NSTAGE -eq 12 ]]"
check "未知参数会报错"         "! $PIPE --not-a-flag"

# ---------------------------------------------------------------- 2. dry-run
echo "[2] dry-run 全流程"
LOG="$TMP/dry.log"
DRY=$(dry_to "$LOG" -1 "$R1" -2 "$R2" -o "$OUT" -t 4)

check "dry-run 能跑完"           "grep -q '全部阶段结束' '$LOG'"
check "打印了 bowtie2"           "grep -q 'bowtie2 --local' '$LOG'"
check "打印了 alignmentSieve"    "grep -q 'alignmentSieve' '$LOG'"
check "打印了 macs2 三次"        "[[ \$(grep -c 'macs2 callpeak' '$LOG') -eq 3 ]]"
check "打印了 bedtools 去黑名单"  "grep -q 'bedtools intersect -v' '$LOG'"
check "打印了 bamCoverage"       "grep -q 'bamCoverage' '$LOG'"
check "建了 10 个输出目录"        "[[ \$(ls -d $OUT/0*_* $OUT/10_report 2>/dev/null | wc -l) -ge 10 ]]"
check "dry-run 不落 .done"       "[[ -z \"\$(ls $OUT/.status/*.done 2>/dev/null)\" ]]"
check "NFR / mono / multi 三档"   "[[ \$(grep -cE '\.(NFR|mono|multi)\.bam' '$LOG') -ge 3 ]]"

# ---------------------------------------------------------------- 3. 阶段控制
echo "[3] 阶段控制"
L3="$TMP/stage.log"
dry_to "$L3" -1 "$R1" -2 "$R2" -o "$OUT" --stages peak >/dev/null
check "--stages peak 不跑 fastqc"  "! grep -q 'fastqc' '$L3'"
check "非法阶段名会报错"            "! $PIPE -1 $R1 -2 $R2 -o "$OUT" --stages notastage --dry-run"

BAMF="$TMP/Case1.srt2.bam"; : > "$BAMF"
L4="$TMP/bam.log"
dry_to "$L4" --bam "$BAMF" -o "$OUT" >/dev/null
check "--bam 默认不跑 bowtie2"     "! grep -q 'bowtie2 --local' '$L4'"
check "--bam 仍会 call peak"       "grep -q 'macs2 callpeak' '$L4'"

# ---------------------------------------------------------------- 4. 参数传递
echo "[4] 参数传递"
L5="$TMP/p.log"
dry_to "$L5" -1 "$R1" -2 "$R2" -o "$OUT" -t 12 --min-mapq 30 --keep-mt --nfr-max 150 >/dev/null
check "-t 影响 bowtie2 线程"      "grep -qE 'bowtie2 .*-p 12' '$L5'"
check "--min-mapq 影响 samtools"   "grep -qE 'samtools view .*-q 30' '$L5'"
check "--keep-mt 关掉过滤"         "grep -q 'drop=0' '$L5'"
check "--nfr-max 影响分档阈值"     "grep -q 'nfrmax=150' '$L5'"

L6="$TMP/nomotif.log"
dry_to "$L6" -1 "$R1" -2 "$R2" -o "$OUT" --no-motif >/dev/null
check "--no-motif 跳过 motif"      "grep -q 'RUN_MOTIF=0' '$L6'"

# ---------------------------------------------------------------- 5. 物种配置
echo "[5] 物种配置"
check "默认物种是 hg38"           "grep -q 'hg38' '$LOG'"
check "不存在的物种会报错"         "! $PIPE -1 $R1 -2 $R2 -o $OUT -g mm10 --dry-run"

# ---------------------------------------------------------------- 6. deploy
echo "[6] deploy"
check "deploy --help"             "bash deploy/deploy.sh --help | grep -q 'doctor'"
check "deploy status 能跑"        "bash deploy/deploy.sh status"
check "deploy doctor 能跑"        "bash deploy/deploy.sh doctor"

# ---------------------------------------------------------------- 7. 错误处理
echo "[7] 错误处理"
check "不给输入会报错"            "! $PIPE --dry-run"
check "R1 单独给会报错"           "! $PIPE -1 $R1 --dry-run"

echo
printf '通过 %d 项，失败 %d 项\n' "$PASS" "$FAIL"
(( FAIL == 0 )) || exit 1
