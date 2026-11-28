################################################################################
# 01_preprocess.smk —— ① FastQC（原始） ② fastp（接头 + 质量）
################################################################################

rule fastqc_raw:
    """① 原始质检：接头污染、质量分布、GC、duplication。"""
    input:
        r1 = res(D_INPUT, "{sample}_R1.fastq.gz"),
        r2 = res(D_INPUT, "{sample}_R2.fastq.gz"),
    output:
        html1 = res(D_QC, "{sample}_R1_fastqc.html"),
        html2 = res(D_QC, "{sample}_R2_fastqc.html"),
        zip1  = res(D_QC, "{sample}_R1_fastqc.zip"),
        zip2  = res(D_QC, "{sample}_R2_fastqc.zip"),
    params:
        outdir = res(D_QC),
    log:
        os.path.join(LOGS, "{sample}.fastqc_raw.log"),
    benchmark:
        bench("fastqc_raw", "{sample}"),
    threads: tcap(4)
    resources:
        mem_mb = 2000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        fastqc -t {threads} --noextract -o {params.outdir} {input.r1} {input.r2} 2> {log}
        """


rule fastp:
    """② 接头修剪 + 质量过滤。

    ATAC 的两个注意点：
      · 不开 --trim_poly_g（非 NovaSeq 建库没有 polyG 偏好）
      · -q 不压到 20 以上，否则 Tn5 插入位点附近的低质量碱基被削掉，切割位点会偏移
    """
    input:
        r1 = res(D_INPUT, "{sample}_R1.fastq.gz"),
        r2 = res(D_INPUT, "{sample}_R2.fastq.gz"),
    output:
        r1   = res(D_TRIM, "{sample}.clean_R1.fastq.gz"),
        r2   = res(D_TRIM, "{sample}.clean_R2.fastq.gz"),
        json = res(D_TRIM, "{sample}.fastp.json"),
        html = res(D_TRIM, "{sample}.fastp.html"),
        metrics = res(D_METRIC, "02_fastp.{sample}.tsv"),
    params:
        extra = config["fastp"]["extra"],
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "{sample}.fastp.log"),
    benchmark:
        bench("fastp", "{sample}"),
    threads: tcap(16)
    resources:
        mem_mb = 6000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              {params.extra} -w {threads} \
              -j {output.json} -h {output.html} 2> {log}

        python3 {params.scripts}/fastp_stats.py \
            --json {output.json} --out {output.metrics} --step 02_fastp
        """


rule fastqc_clean:
    """修剪后再看一遍：确认接头掉了、质量上来了，且没有过度修剪。"""
    input:
        r1 = res(D_TRIM, "{sample}.clean_R1.fastq.gz"),
        r2 = res(D_TRIM, "{sample}.clean_R2.fastq.gz"),
    output:
        html1 = res(D_QC, "clean", "{sample}_clean_R1_fastqc.html"),
        html2 = res(D_QC, "clean", "{sample}_clean_R2_fastqc.html"),
        zip1  = res(D_QC, "clean", "{sample}_clean_R1_fastqc.zip"),
        zip2  = res(D_QC, "clean", "{sample}_clean_R2_fastqc.zip"),
    params:
        outdir = res(D_QC, "clean"),
    log:
        os.path.join(LOGS, "{sample}.fastqc_clean.log"),
    benchmark:
        bench("fastqc_clean", "{sample}"),
    threads: tcap(4)
    resources:
        mem_mb = 2000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        fastqc -t {threads} --noextract -o {params.outdir} {input.r1} {input.r2} 2> {log}
        """
