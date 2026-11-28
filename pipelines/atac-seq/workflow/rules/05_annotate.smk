################################################################################
# 05_annotate.smk —— ⑫ HOMER peak 注释 ⑬ HOMER known motif 富集
#
# 这两步依赖 HOMER 的基因组索引。没装时 00_common.smk 会把它们从 DAG 里摘掉，
# 而不是让整条链在第 12 步崩掉——已经跑完的比对结果不该因为注释工具缺失而作废。
################################################################################

rule homer_annotate:
    """⑫ annotatePeaks.pl：每个 peak 落在 promoter / exon / intron / intergenic 哪里。

    -annStats 会单独写一份分布统计，比解析主表可靠。
    """
    input:
        peaks = res(D_PEAK, "{sample}.all.filtered.narrowPeak"),
    output:
        anno = res(D_ANN, "{sample}", "annotated_peaks.txt"),
        stats = res(D_ANN, "{sample}", "annostats.txt"),
        metrics = res(D_METRIC, "12_annotate.{sample}.tsv"),
    params:
        genome = homer_genome,
    log:
        os.path.join(LOGS, "{sample}.annotatePeaks.log"),
    benchmark:
        bench("homer_annotate", "{sample}"),
    threads: tcap(4)
    resources:
        mem_mb = 8000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        annotatePeaks.pl {input.peaks} {params.genome} \
            -annStats {output.stats} > {output.anno} 2> {log}

        n=$(wc -l < {output.anno})
        printf '12_annotate\tannotated_peaks\t%s\tpeaks\t\n' "$n" > {output.metrics}
        """


rule annot_dist:
    """⑫ 把 annStats 的文本块转成表 + 条形图。

    ATAC 的 peak 应以 distal intergenic + intron 为主；promoter 超过 60%
    就要怀疑 TSS 附近的高信号把 peak 全吸过去了。
    """
    input:
        stats = res(D_ANN, "{sample}", "annostats.txt"),
    output:
        tsv = res(D_ANN, "{sample}", "annotation_distribution.tsv"),
        png = res(D_ANN, "{sample}", "annotation_distribution.png"),
    params:
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "{sample}.annot_dist.log"),
    benchmark:
        bench("annot_dist", "{sample}"),
    threads: 1
    resources:
        mem_mb = 2000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        python3 {params.scripts}/annot_dist.py --annstats {input.stats} \
            --out {output.tsv} --png {output.png} 2> {log}
        """


rule toppeaks:
    """⑬ 按 -log10(q) 取 top N。

    全量 7 万 peak 跑 known motif 要 2 小时，取信号最强的 1 万既省时间，
    又不改变富集结论（motif 富集是被最强的那批 peak 驱动的）。
    """
    input:
        peaks = res(D_PEAK, "{sample}.{frac}.filtered.narrowPeak"),
    output:
        top = res(D_PEAK, "{sample}.{frac}.toppeaks.narrowPeak"),
    params:
        n = config["motif"]["topn"],
    wildcard_constraints:
        frac = "all|nfr|mono",
    log:
        os.path.join(LOGS, "{sample}.{frac}.toppeaks.log"),
    benchmark:
        bench("toppeaks", "{sample}.{frac}"),
    threads: 1
    resources:
        mem_mb = 2000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        # 不用 sort | head：head 提前退出会让 sort 收到 SIGPIPE，pipefail 下直接中断
        sort -k9,9nr {input.peaks} -o {output.top}.sorted
        head -n {params.n} {output.top}.sorted > {output.top}
        rm -f {output.top}.sorted
        """


rule homer_motif:
    """⑬ findMotifsGenome.pl -mknown：只比对已知 motif 库。

    不加 -mknown 的话 HOMER 会顺带做 de novo discovery，那一步要几小时。
    """
    input:
        top = res(D_PEAK, "{sample}.{frac}.toppeaks.narrowPeak"),
    output:
        known = res(D_MOTIF, "{sample}", "{frac}", "knownResults.txt"),
        html = res(D_MOTIF, "{sample}", "{frac}", "knownResults.html"),
        metrics = res(D_METRIC, "13_motif.{sample}.{frac}.tsv"),
    params:
        genome = homer_genome,
        size = config["motif"]["size"],
        known = homer_known_motifs,
        outdir = res(D_MOTIF, "{sample}", "{frac}"),
    wildcard_constraints:
        frac = "all|nfr|mono",
    log:
        os.path.join(LOGS, "{sample}.{frac}.motif.log"),
    benchmark:
        bench("homer_motif", "{sample}.{frac}"),
    threads: tcap(8)
    resources:
        mem_mb = 12000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        findMotifsGenome.pl {input.top} {params.genome} {params.outdir} \
            -size {params.size} -mask -mknown {params.known} > {log} 2>&1

        n=$(awk 'NR>1' {output.known} | wc -l)
        best=$(awk -F'\t' 'NR==2{{print $1}}' {output.known})
        {{
          printf '13_motif\tmotifs_tested_%s\t%s\tmotifs\t\n' {wildcards.frac} "$n"
          printf '13_motif\ttop_motif_%s\t%s\t\t\n' {wildcards.frac} "$best"
        }} > {output.metrics}
        """


rule motif_top:
    """⑬ 整理成一张能直接看的表（按 fold 排序，fold 比 p 值更能说明富集强度）。"""
    input:
        known = res(D_MOTIF, "{sample}", "{frac}", "knownResults.txt"),
    output:
        top = res(D_MOTIF, "{sample}", "{frac}", "top_motifs.tsv"),
    wildcard_constraints:
        frac = "all|nfr|mono",
    params:
        scripts = SCRIPTS,
    log:
        os.path.join(LOGS, "{sample}.{frac}.motif_top.log"),
    benchmark:
        bench("motif_top", "{sample}.{frac}"),
    threads: 1
    resources:
        mem_mb = 2000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        python3 {params.scripts}/motif_top.py --known {input.known} \
            --out {output.top} --top 15 2> {log}
        """
