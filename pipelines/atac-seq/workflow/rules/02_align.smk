################################################################################
# 02_align.smk —— ③ bowtie2 ④ 过滤 ⑤ Tn5 校正 ⑥ 重排 ⑦ 片段分箱
################################################################################

rule bowtie2:
    """③ 比对。SAM 是作业图里的产物；keep_sam=false 时标为 temp，跑完 ④ 自动删。

    --local 允许局部比对（reads 常在低复杂度处截断），--very-fast-local 提速；
    -X 2000 容纳跨多个核小体的长片段。
    """
    input:
        r1 = res(D_TRIM, "{sample}.clean_R1.fastq.gz"),
        r2 = res(D_TRIM, "{sample}.clean_R2.fastq.gz"),
    output:
        sam = (temp(res(D_ALN, "{sample}.bowtie2.sam"))
               if not config.get("keep_sam", False)
               else res(D_ALN, "{sample}.bowtie2.sam")),
        txt = res(D_ALN, "{sample}.bowtie2.txt"),
        metrics = res(D_METRIC, "03_bowtie2.{sample}.tsv"),
    params:
        opts = config["bowtie2"]["extra"],
        insert_max = config["bowtie2"]["insert_max"],
        index = bowtie2_index,
    log:
        os.path.join(LOGS, "{sample}.bowtie2.log"),
    benchmark:
        bench("bowtie2", "{sample}"),
    threads: tcap(20)
    resources:
        mem_mb = 12000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        bowtie2 {params.opts} -X {params.insert_max} \
            -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            --rg-id {wildcards.sample} --rg SM:{wildcards.sample} \
            --rg PL:ILLUMINA --rg DS:SRR13849313 \
            -S {output.sam} 2> {output.txt}

        # 一次遍历同时拿总 reads 和 chrM reads
        counts=$(awk 'BEGIN{{t=0;m=0}} !/^@/{{t++; if($3=="chrM") m++}} END{{printf "%d %d", t, m}}' {output.sam})
        total=$(echo $counts | cut -d' ' -f1)
        chrm=$(echo $counts | cut -d' ' -f2)
        rate=$(grep -o '[0-9.]*% overall alignment rate' {output.txt} | grep -o '^[0-9.]*' || echo NA)

        {{
          printf '03_bowtie2\taligned_reads\t%s\treads\t\n' "$total"
          printf '03_bowtie2\tchrm_reads\t%s\treads\t\n' "$chrm"
          printf '03_bowtie2\tchrm_fraction\t%.4f\tratio\t<0.05\n' \
              $(awk -v a="$chrm" -v b="$total" 'BEGIN{{print (b>0)? a/b : 0}}')
          printf '03_bowtie2\toverall_alignment_rate\t%s\t%%\t>95\n' "$rate"
        }} > {output.metrics}
        """


rule bam_filter:
    """④ MAPQ≥20 + 去 chrM/组装碎片 + 排序 + 索引。

    用 samtools 的过滤表达式而不是 grep：grep 会误伤 read name，
    也拦不住 chrUn*/random 这些同样该丢的 contig。
    """
    input:
        sam = res(D_ALN, "{sample}.bowtie2.sam"),
    output:
        bam = res(D_BAM, "{sample}.srt.bam"),
        bai = res(D_BAM, "{sample}.srt.bam.bai"),
        flagstat = res(D_BAM, "{sample}.srt.flagstat.txt"),
        stats = res(D_BAM, "{sample}.srt.stats.txt"),
        metrics = res(D_METRIC, "04_bam_filter.{sample}.tsv"),
    params:
        expr = contig_filter(),
        mapq = MAPQ,
        mem = config.get("sort_mem", "2G"),
        tmpdir = res(D_BAM),
    log:
        os.path.join(LOGS, "{sample}.bam_filter.log"),
    benchmark:
        bench("bam_filter", "{sample}"),
    threads: tcap(16)
    resources:
        mem_mb = 16000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        samtools view -@ {threads} -b -q {params.mapq} -e '{params.expr}' {input.sam} \
            | samtools sort -@ {threads} -m {params.mem} -T {params.tmpdir}/.tmp_{wildcards.sample} \
                -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        samtools stats  -@ {threads} {output.bam} > {output.stats}

        kept=$(samtools view -c {output.bam})
        proper=$(awk '/properly paired/{{print $1; exit}}' {output.flagstat})
        dup=$(awk -F'\t' '$1=="SN" && index($2,"reads duplicated:")==1 {{print $3; exit}}' {output.stats})
        isize=$(awk -F'\t' '$1=="SN" && index($2,"insert size average:")==1 {{print $3; exit}}' {output.stats})

        {{
          printf '04_bam_filter\treads_after_filter\t%s\treads\t\n' "$kept"
          printf '04_bam_filter\tproperly_paired\t%s\treads\t\n' "$proper"
          printf '04_bam_filter\tduplicate_reads\t%s\treads\t\n' "$dup"
          printf '04_bam_filter\tinsert_size_mean\t%s\tbp\t~200\n' "$isize"
          printf '04_bam_filter\tduplicate_fraction\t%.4f\tratio\t<0.30\n' \
              $(awk -v d="$dup" -v t="$kept" 'BEGIN{{print (t>0)? d/t : 0}}')
        }} > {output.metrics}
        """


rule atac_shift:
    """⑤ Tn5 末端校正：+4 / -5。

    这一步只改坐标，不改 TLEN —— 所以 ⑦ 的分箱仍然反映真实片段长度。
    """
    input:
        bam = res(D_BAM, "{sample}.srt.bam"),
        bai = res(D_BAM, "{sample}.srt.bam.bai"),
    output:
        bam = res(D_BAM, "{sample}.shift.bam"),
    params:
        flag = "--ATACshift" if config.get("atac_shift", True) else "",
    log:
        os.path.join(LOGS, "{sample}.atac_shift.log"),
    benchmark:
        bench("atac_shift", "{sample}"),
    threads: tcap(16)
    resources:
        mem_mb = 8000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        alignmentSieve -p {threads} {params.flag} -b {input.bam} -o {output.bam} 2> {log}
        """


rule resort:
    """⑥ 重排 + 建索引。

    alignmentSieve 只保留前 11 个字段，且坐标挪动后顺序不再有序 ——
    不重排直接建索引会得到错误的区间树，后面所有区间查询都会错。
    """
    input:
        bam = res(D_BAM, "{sample}.shift.bam"),
    output:
        bam = res(D_BAM, "{sample}.srt2.bam"),
        bai = res(D_BAM, "{sample}.srt2.bam.bai"),
        check = res(D_BAM, "{sample}.field_check.txt"),
        metrics = res(D_METRIC, "06_resort.{sample}.tsv"),
    params:
        mem = config.get("sort_mem", "2G"),
        tmpdir = res(D_BAM),
    log:
        os.path.join(LOGS, "{sample}.resort.log"),
    benchmark:
        bench("resort", "{sample}"),
    threads: tcap(16)
    resources:
        mem_mb = 16000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        samtools sort -@ {threads} -m {params.mem} -T {params.tmpdir}/.tmp2_{wildcards.sample} \
            -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}

        # head 提前关管道会让 samtools 收到 SIGPIPE，pipefail 下要兜住
        nf=$(samtools view {output.bam} 2>/dev/null | head -n 1 | awk '{{print NF}}' || true)
        printf 'first_record_fields\t%s\n' "$nf" > {output.check}

        printf '06_resort\trecord_fields\t%s\tcolumns\t11\n' "$nf" > {output.metrics}
        printf '06_resort\treads_after_shift\t%s\treads\t\n' \
            $(samtools view -c {output.bam}) >> {output.metrics}
        """


rule frag_bin:
    """⑦ 按 |TLEN| 分 NFR / 单核小体 / 多核小体。

    必须在已偏移的 srt2.bam 上分：147 bp 阈值对应"缠一个核小体"，
    偏移把两端各挪了 4–5 bp，用未偏移的 BAM 分档会让 NFR 档系统性偏少。
    """
    input:
        bam = res(D_BAM, "{sample}.srt2.bam"),
        bai = res(D_BAM, "{sample}.srt2.bam.bai"),
    output:
        nfr  = res(D_BAM, "{sample}_nucleosome_free.bam"),
        nfr_bai = res(D_BAM, "{sample}_nucleosome_free.bam.bai"),
        mono = res(D_BAM, "{sample}_mono_nucleosome.bam"),
        mono_bai = res(D_BAM, "{sample}_mono_nucleosome.bam.bai"),
        poly = res(D_BAM, "{sample}_multi_nucleosome.bam"),
        poly_bai = res(D_BAM, "{sample}_multi_nucleosome.bam.bai"),
        hist = res(D_FRAG, "{sample}.fragsize_hist.tsv"),
        stats = res(D_FRAG, "{sample}.fragbin_stats.tsv"),
        png = res(D_FRAG, "{sample}.fragsize_dist.png"),
        metrics = res(D_METRIC, "07_fragbin.{sample}.tsv"),
    params:
        nfr_max = config["tlen"]["nfr_max"],
        mono_max = config["tlen"]["mono_max"],
        outdir = res(D_FRAG),
        bamdir = res(D_BAM),
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "{sample}.frag_bin.log"),
    benchmark:
        bench("frag_bin", "{sample}"),
    threads: tcap(8)
    resources:
        mem_mb = 8000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        python3 {params.scripts}/frag_bin.py \
            --bam {input.bam} --outdir {params.outdir} --bamdir {params.bamdir} \
            --sample {wildcards.sample} \
            --nfr-max {params.nfr_max} --mono-max {params.mono_max} \
            --threads {threads} 2> {log}

        # fragbin_stats.tsv → 标准 5 列指标
        awk -F'\t' 'NR>1 {{
            printf "07_fragbin\t%s_fraction\t%s\tratio\t0.3-0.6\n", $1, $3
            printf "07_fragbin\t%s_fragments\t%s\treads\t\n", $1, $2
        }}' {output.stats} > {output.metrics}
        """
