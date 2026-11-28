################################################################################
# 04_qc.smk —— ⑩ FRiP ⑪ bigWig + TSS 富集 ⑭ ATACseqQC
################################################################################

rule frip:
    """⑩ FRiP：samtools 与 bedtools 两条独立实现各算一遍。

    两个数应该几乎一致；差得多说明 BED 与 BAM 的染色体命名/坐标体系/排序
    对不上，这时候不该继续解读 FRiP，而该先查文件。
    """
    input:
        bam = lambda w: frac_bam(w),
        peaks = res(D_PEAK, "{sample}.{frac}.filtered.narrowPeak"),
    output:
        frip = res(D_FRIP, "{sample}", "{frac}.frip.tsv"),
        metrics = res(D_METRIC, "10_frip.{sample}.{frac}.tsv"),
    params:
        bedtools = "--bedtools" if config.get("frip_with_bedtools", True) else "",
        outdir = res(D_FRIP, "{sample}"),
        scripts = SCRIPTS,
    wildcard_constraints:
        frac = "all|nfr|mono",
    log:
        os.path.join(LOGS, "{sample}.{frac}.frip.log"),
    benchmark:
        bench("frip", "{sample}.{frac}"),
    threads: tcap(4)
    resources:
        mem_mb = 4000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        python3 {params.scripts}/frip.py \
            --bam {input.bam} --peaks {input.peaks} --label {wildcards.frac} \
            --outdir {params.outdir} {params.bedtools} \
            --threads {threads} 2> {log}

        awk -F'\t' -v OFS='\t' -v f={wildcards.frac} \
            '{{print "10_frip", f"_"$2, $3, $4, $5}}' \
            {output.frip} > {output.metrics}
        """


rule plot_fingerprint:
    """⑩ plotFingerprint：曲线偏离对角线越多，signal-to-noise 越好。

    这是判断"peak 是不是真的"最直观的一张图——如果曲线几乎贴对角线，
    说明 reads 在基因组上是均匀铺开的，根本没有富集。
    """
    input:
        bam = res(D_BAM, "{sample}.srt2.bam"),
    output:
        png = res(D_FRIP, "{sample}.plotFingerprint.png"),
        qc = res(D_FRIP, "{sample}.plotFingerprint.qc.tsv"),
    log:
        os.path.join(LOGS, "{sample}.plotFingerprint.log"),
    benchmark:
        bench("plot_fingerprint", "{sample}"),
    threads: tcap(8)
    resources:
        mem_mb = 6000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        plotFingerprint -b {input.bam} --labels {wildcards.sample} \
            --plotFile {output.png} --outQualityMetrics {output.qc} \
            --numberOfSamples 100000 -p {threads} 2> {log}
        """


rule bam_coverage:
    """⑪ bamCoverage → bigWig（CPM 归一化）。

    --ignoreDuplicates + 黑名单：PCR 重复和高拷贝区都会让信号图出现"假高峰"。
    --extendReads 对双端数据按真实片段延伸，这才是 ATAC 的信号单位。
    """
    input:
        bam = res(D_BAM, "{sample}.srt2.bam"),
        bai = res(D_BAM, "{sample}.srt2.bam.bai"),
    output:
        bw = res(D_COV, "{sample}.bamCoverage.bw"),
    params:
        binsize = config["bwBinSize"],
        norm = config["normalizeUsing"],
        blacklist = "--blackListFileName %s" % blacklist_bed,
        ignore = ("--ignoreForNormalization %s" % ignoreForNormalization
                  if ignoreForNormalization else ""),
    log:
        os.path.join(LOGS, "{sample}.bamCoverage.log"),
    benchmark:
        bench("bam_coverage", "{sample}"),
    threads: tcap(16)
    resources:
        mem_mb = 8000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        bamCoverage -p {threads} --ignoreDuplicates {params.blacklist} {params.ignore} \
            -bs {params.binsize} --normalizeUsing {params.norm} --extendReads \
            -b {input.bam} -o {output.bw} 2> {log}
        """


rule tss_prep:
    """⑪ 从 GTF 取每个基因的 TSS（带链方向）。

    只做一次，所有样本共用；所以这条 rule 没有 {sample} 通配符。
    """
    input:
        gtf = genes_gtf,
    output:
        bed = res(D_TSS, "genes.tss.bed"),
    params:
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "tss_prep.log"),
    benchmark:
        bench("tss_prep"),
    threads: 1
    resources:
        mem_mb = 2000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        python3 {params.scripts}/tss_prep.py --gtf {input.gtf} --out {output.bed} 2> {log}
        """


rule tss_matrix:
    """⑪ TSS ±2 kb 的信号矩阵。

    --referencePoint TSS 会按 BED 第 6 列的链方向翻转，负链基因不会算反。
    """
    input:
        bw = res(D_COV, "{sample}.bamCoverage.bw"),
        bed = res(D_TSS, "genes.tss.bed"),
    output:
        matrix = res(D_TSS, "{sample}.tss_matrix.gz"),
    params:
        flank = config["tss"]["flank"],
        binsize = config["tss"]["bin_size"],
    log:
        os.path.join(LOGS, "{sample}.computeMatrix.log"),
    benchmark:
        bench("tss_matrix", "{sample}"),
    threads: tcap(8)
    resources:
        mem_mb = 10000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        computeMatrix reference-point -S {input.bw} -R {input.bed} \
            -b {params.flank} -a {params.flank} -bs {params.binsize} \
            --referencePoint TSS --missingDataAsZero \
            -p {threads} -o {output.matrix} 2> {log}
        """


rule tss_profile:
    """⑪ 出图。图给人看，数给阈值比，两样都要有。"""
    input:
        matrix = res(D_TSS, "{sample}.tss_matrix.gz"),
    output:
        pdf = res(D_TSS, "{sample}.tss_profile.pdf"),
        data = res(D_TSS, "{sample}.tss_profile_data.tsv"),
    log:
        os.path.join(LOGS, "{sample}.plotProfile.log"),
    benchmark:
        bench("tss_profile", "{sample}"),
    threads: 1
    resources:
        mem_mb = 4000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        plotProfile -m {input.matrix} -o {output.pdf} \
            --outFileNameData {output.data} \
            --plotTitle "TSS enrichment ({wildcards.sample})" \
            --yAxisLabel "mean coverage (CPM)" 2> {log}
        """


rule tss_score:
    """⑪ TSS 富集倍数 = TSS ±50 bp 均值 / 两端 flank 均值。

    和 ATACseqQC 的 TSSEscore 是同一个定义，⑭ 可以拿来做交叉验证。
    """
    input:
        matrix = res(D_TSS, "{sample}.tss_matrix.gz"),
    output:
        score = res(D_TSS, "{sample}.tss_score.tsv"),
        metrics = res(D_METRIC, "11_coverage.{sample}.tsv"),
    params:
        flank = config["tss"]["flank"],
        window = config["tss"]["score_window"],
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "{sample}.tss_score.log"),
    benchmark:
        bench("tss_score", "{sample}"),
    threads: 1
    resources:
        mem_mb = 4000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        python3 {params.scripts}/tss_score.py --matrix {input.matrix} \
            --flank {params.flank} --window {params.window} \
            --out {output.score} 2> {log}

        score=$(awk -F'\t' '$1=="tss_enrichment_score"{{print $2}}' {output.score})
        printf '11_coverage\ttss_enrichment_score\t%s\tfold\t>6\n' "$score" > {output.metrics}
        """


rule atacqc:
    """⑭ ATACseqQC：片段分布 + TSSE。

    用未偏移的 srt.bam：TSSE 和片段分布都基于真实片段长度，
    而 alignmentSieve 挪了坐标却不会同步改 TLEN，两者会自相矛盾。
    """
    input:
        bam = res(D_BAM, "{sample}.srt.bam"),
    output:
        metrics = res(D_ATACQC, "{sample}.atacqc_metrics.tsv"),
        pdf = res(D_ATACQC, "{sample}.fragSizeDist.pdf"),
        summary = res(D_METRIC, "14_atacqc.{sample}.tsv"),
    params:
        outdir = res(D_ATACQC),
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "{sample}.atacqc.log"),
    benchmark:
        bench("atacqc", "{sample}"),
    threads: tcap(4)
    resources:
        mem_mb = 12000,
    conda:
        "../envs/atacqc.yaml"
    shell:
        """
        Rscript {params.scripts}/atacqc.R --bam {input.bam} \
            --outdir {params.outdir} --sample {wildcards.sample} 2> {log}

        awk -F'\t' 'NR>1 {{print "14_atacqc", $1, $2, "", ""}}' \
            {output.metrics} > {output.summary}
        """
