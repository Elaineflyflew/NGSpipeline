################################################################################
# 06_report.smk —— ⑮ IGV 截图 + 汇总报告
################################################################################

rule igv_snapshot:
    """⑮ MYC 位点截图。

    用 batch 脚本而不是手动点：可复现、可进 CI。
    服务器无图形环境时用 xvfb-run；没有 IGV 就退回 pyGenomeTracks 出同一区域。
    """
    input:
        bw = res(D_COV, "{sample}.bamCoverage.bw"),
        peaks = res(D_PEAK, "{sample}.all.filtered.narrowPeak"),
    output:
        png = res(D_IGV, "{sample}.MYC.png"),
        batch = res(D_IGV, "{sample}.igv_batch.txt"),
    params:
        locus = config["igv"]["locus"],
        width = config["igv"]["width"],
        height = config["igv"]["height"],
        genome = homer_genome,
        tracks = res(D_IGV, "{sample}") + ".tracks.ini",
    log:
        os.path.join(LOGS, "{sample}.igv.log"),
    benchmark:
        bench("igv_snapshot", "{sample}"),
    threads: 1
    resources:
        mem_mb = 4000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        cat > {output.batch} <<EOF
        new
        genome {params.genome}
        load {input.bw}
        load {input.peaks}
        goto {params.locus}
        sort position
        snapshot {output.png}
        exit
        EOF

        if command -v igv >/dev/null 2>&1; then
            if command -v xvfb-run >/dev/null 2>&1; then
                xvfb-run -a igv -b {output.batch} -g {params.genome} > {log} 2>&1
            else
                igv -b {output.batch} -g {params.genome} > {log} 2>&1
            fi
        elif command -v pyGenomeTracks >/dev/null 2>&1; then
            echo "没有 IGV，改用 pyGenomeTracks" > {log}
            cat > {params.tracks} <<EOF
        [bigwig]
        file = {input.bw}
        title = {wildcards.sample} (CPM)
        height = 6
        color = #8B0000

        [peaks]
        file = {input.peaks}
        title = MACS2 peaks
        height = 2
        color = #1f4e79
        EOF
            pyGenomeTracks --tracks {params.tracks} \
                --region {params.locus} -o {output.png} >> {log} 2>&1
        else
            echo "既没有 igv 也没有 pyGenomeTracks；batch 已生成，可本地跑" >&2
            exit 1
        fi
        """


rule collect_metrics:
    """汇总所有 metrics/*.tsv，按阈值打 PASS/WARN/FAIL。

    依赖挂在几个"末端产物"上：它们存在就说明前面每一步的指标都已落盘，
    比用 glob 猜依赖可靠（glob 在 DAG 构建时还看不到文件）。
    """
    input:
        expand(res(D_BAM, "{sample}.srt2.bam"), sample=SAMPLES),
        expand(res(D_FRAG, "{sample}.fragsize_dist.png"), sample=SAMPLES),
        expand(res(D_PEAK, "{sample}.{frac}.filtered.narrowPeak"),
               sample=SAMPLES, frac=FRACTIONS),
        expand(res(D_FRIP, "{sample}", "{frac}.frip.tsv"),
               sample=SAMPLES, frac=FRACTIONS),
        expand(res(D_TSS, "{sample}.tss_score.tsv"), sample=SAMPLES),
        *annot_metric_deps(),
    output:
        summary = res(D_REPORT, "summary.tsv"),
    params:
        metricdir = res(D_METRIC),
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "collect_metrics.log"),
    benchmark:
        bench("collect_metrics"),
    threads: 1
    resources:
        mem_mb = 2000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        python3 {params.scripts}/collect_metrics.py --metrics {params.metricdir} \
            --out {output.summary} 2> {log}
        """


rule multiqc:
    """把所有步骤的 QC 报告聚成一个 HTML（fastqc / fastp / samtools stats）。"""
    input:
        expand(res(D_QC, "{sample}_R1_fastqc.zip"), sample=SAMPLES),
        expand(res(D_QC, "{sample}_R2_fastqc.zip"), sample=SAMPLES),
        expand(res(D_TRIM, "{sample}.fastp.json"), sample=SAMPLES),
        expand(res(D_BAM, "{sample}.srt.stats.txt"), sample=SAMPLES),
    output:
        html = res(D_REPORT, "multiqc_report.html"),
    params:
        outdir = res(D_REPORT),
        results = RES,
    log:
        os.path.join(LOGS, "multiqc.log"),
    benchmark:
        bench("multiqc"),
    threads: 1
    resources:
        mem_mb = 4000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        multiqc -f -o {params.outdir} -n multiqc_report.html {params.results} 2> {log}
        """


rule report:
    """把 summary.tsv 和所有图拼成一份自包含 HTML。

    这份报告回答的不是"跑完了吗"，而是"跑对了吗"：每步的产物在哪、
    指标多少、和阈值比是什么结论。
    """
    input:
        summary = res(D_REPORT, "summary.tsv"),
    output:
        html = res(D_REPORT, "report.html"),
    params:
        sample = SAMPLES[0],
        outdir = res(D_REPORT),
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "report.log"),
    benchmark:
        bench("report"),
    threads: 1
    resources:
        mem_mb = 4000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        python3 {params.scripts}/make_report.py --summary {input.summary} \
            --sample {params.sample} --outdir {params.outdir} 2> {log}
        """
