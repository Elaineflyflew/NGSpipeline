#!/usr/bin/env bash
# =============================================================================
#  src/70_preflight.sh —— 预检
#
#  只做"能不能跑"的判断，不碰数据。分四类：
#    1. 工具：先按启用阶段收集，再去 PATH 之外的地方找（别的 conda env 里）
#    2. 参考基因组索引是否齐全（fai / dict / bwt）
#    3. 命名一致性（chr17 vs 17）—— 这是最隐蔽的坑：不报错，但结果全错
#    4. 数据库是否就位（bundle / vep cache / funcotator）
#
#  设计原则：已经在机器上的东西直接用，不要求用户为它再装一遍。
# =============================================================================

# 本流程需要的外部命令，按启用阶段和 caller/annotator 组合收集，去重后逐行输出。
# 用子 shell 输出而不是全局数组，调用方可以 while read 循环，也方便测试。
required_tools() {
  local -a out=()
  local c a t

  _rt_add() {
    for t in ${out[@]+"${out[@]}"}; do [[ "$t" == "$1" ]] && return 0; done
    out+=("$1")
  }

  stage_enabled qc     && _rt_add fastqc
  stage_enabled trim   && _rt_add fastp
  stage_enabled align  && { _rt_add bwa; _rt_add samtools; }
  stage_enabled bam    && _rt_add samtools
  stage_enabled report && _rt_add bcftools

  for c in "${CALLERS[@]}"; do
    case "$c" in
      bcftools)    _rt_add bcftools ;;
      gatk)        _rt_add gatk ;;
      deepvariant) [[ -n "$DV_SIF" ]] || _rt_add run_deepvariant ;;
    esac
  done
  for a in ${ANNOTATORS[@]+"${ANNOTATORS[@]}"}; do
    case "$a" in
      vep)        _rt_add vep; _rt_add bgzip; _rt_add tabix ;;
      annovar)    _rt_add table_annovar.pl ;;
      funcotator) _rt_add gatk ;;
    esac
  done

  for t in ${out[@]+"${out[@]}"}; do printf '%s\n' "$t"; done
}

# 把"找到但不在 PATH 里"的工具目录加进来。
# 典型场景：gatk 装在 /home/biotools/miniforge3/envs/gatk/bin 但从没 activate 过。
resolve_tools() {
  local t p d
  local -a dirs=()

  while read -r t; do
    [[ -z "$t" ]] && continue
    command -v "$t" >/dev/null 2>&1 && continue
    p=$(find_tool "$t") || continue
    d=$(dirname "$p")
    case " ${dirs[*]-} " in *" $d "*) continue ;; esac
    dirs+=("$d")
    log "找到 $t: $p"
  done < <(required_tools)

  if (( ${#dirs[@]} )); then
    local i
    for (( i = ${#dirs[@]} - 1; i >= 0; i-- )); do
      export PATH="${dirs[i]}:$PATH"
    done
    log "已加入 PATH: ${dirs[*]}"
  fi
}

preflight() {
  section "预检"

  printf '样本=%s\n参考=%s\n输出=%s\n线程=%s\ncaller=%s\nannotator=%s\n阶段=%s\n\n' \
    "$SAMPLE" "$REF" "$OUTDIR" "$THREADS" "${CALLERS[*]}" "$ANNOT_DESC" "$STAGES"

  # ---- 1. 工具：先找（PATH 之外的也算），再判断缺不缺 ----
  resolve_tools

  local t p
  while read -r t; do
    [[ -z "$t" ]] && continue
    if (( DRY_RUN )); then
      printf '[dry-run] 需要命令: %s\n' "$t"
      continue
    fi
    if p=$(find_tool "$t"); then
      (( CHECK_ONLY )) && printf '  OK   %-18s %s\n' "$t" "$p"
    elif (( CHECK_ONLY )); then
      printf '  MISS %-18s\n' "$t"
      warn "缺少命令: $t"
      MISSING_TOOLS=$(( MISSING_TOOLS + 1 ))
    else
      die "找不到命令 $t —— 跑 'deploy/deploy.sh doctor' 看缺什么，再 'deploy/deploy.sh env' 补装"
    fi
  done < <(required_tools)

  # ---- 2. 参考索引 ----
  [[ -f "${REF}.fai" ]] || warn "缺少 ${REF}.fai（samtools faidx）"
  local dict="${REF%.fa}.dict"; [[ "$REF" == *.fasta ]] && dict="${REF%.fasta}.dict"
  [[ -f "$dict" ]] || warn "缺少 $(basename "$dict")（gatk CreateSequenceDictionary）"
  if stage_enabled align; then
    [[ -f "${REF}.bwt" ]] || warn "缺少 ${REF}.bwt（bwa index）—— deploy/deploy.sh ref 会建好"
  fi

  # ---- 3. 命名一致性 ----
  local style; style=$(ref_style)
  if [[ "$style" == "ensembl" ]] && stage_enabled call; then
    case "$CALLER" in
      gatk|all|bcftools)
        warn "参考用的是 Ensembl 命名（17 而非 chr17）：GATK bundle 和 dbsnp 都是 chr 前缀，BQSR 会找不到染色体"
        warn "解决办法：deploy/deploy.sh ref --source broad（推荐），或 bcftools annotate --rename-chrs 改名"
        ;;
    esac
  fi
  if [[ -n "$BED" && -s "$BED" && -s "${REF}.fai" ]]; then
    local bed_chr ref_chr
    bed_chr=$(head -1 "$BED" | cut -f1)
    ref_chr=$(head -1 "${REF}.fai" | cut -f1)
    if [[ "$bed_chr" == chr* && "$ref_chr" != chr* ]] || [[ "$bed_chr" != chr* && "$ref_chr" == chr* ]]; then
      warn "BED 第一条染色体是 '$bed_chr'，参考是 '$ref_chr' —— 命名不一致，-L/-R 会匹配不到区域"
    fi
  fi

  # ---- 4. 数据库 ----
  if (( ${#ANNOTATORS[@]} )); then
    local miss=0 a
    for a in "${ANNOTATORS[@]}"; do
      case "$a" in
        vep)        [[ -d "$VEP_CACHE/homo_sapiens" ]] || { warn "VEP cache 缺失: $VEP_CACHE/homo_sapiens"; miss=1; } ;;
        annovar)    [[ -d "$ANNOVAR_DB" ]] || { warn "ANNOVAR humandb 缺失: $ANNOVAR_DB"; miss=1; } ;;
        funcotator) [[ -d "$FUNCOTATOR_DATA" ]] || { warn "Funcotator 数据源缺失: ${FUNCOTATOR_DATA:-<未设置>}"; miss=1; } ;;
      esac
    done
    (( miss )) && warn "缺失的注释库会导致对应注释被跳过（不影响变异检测结果）"
  fi
  if [[ "$CALLER" == "gatk" || "$CALLER" == "all" ]]; then
    # 逐个点清，缺哪个说哪个；.vcf / .vcf.gz 都算
    local n=0 miss="" nm f
    for nm in $GATK_BUNDLE_FILES; do
      if f=$(bundle_file "$nm"); then
        n=$(( n + 1 ))
        # 未压缩的 VCF 要 .idx（GATK 建的那种），缺了 BQSR 会直接报错
        if [[ "$f" == *.vcf && ! -s "$f.idx" ]] && (( ! DRY_RUN )); then
          warn "$(basename "$f") 没有 .idx，运行时自动补建（gatk IndexFeatureFile）"
        fi
      else
        miss="$miss $nm"
      fi
    done
    if (( n == 4 )); then
      log "GATK bundle: 4/4 已知位点就绪"
    else
      warn "GATK bundle 缺:${miss} —— 缺 dbsnp/Mills/1000G 跳过 BQSR，缺 hapmap 跳过 Tranche"
    fi
  fi
  if [[ "$CALLER" == "deepvariant" || "$CALLER" == "all" ]]; then
    case "$DV_MODEL" in WES|WGS) ;; *) die "--dv-model 只能是 WES 或 WGS，收到: $DV_MODEL" ;; esac
  fi

  if (( CHECK_ONLY && MISSING_TOOLS > 0 )); then
    warn "预检完成：缺失 ${MISSING_TOOLS} 个命令（见上）"
  else
    log "预检通过"
  fi
}
