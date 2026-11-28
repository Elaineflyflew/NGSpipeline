################################################################################
# 03_peaks.smk —— ⑧ MACS2 peak calling（三条路线） ⑨ 去黑名单
################################################################################

def frac_inputs(wildcards):
    bam = frac_bam(wildcards)
    return {"bam": bam, "bai": bam + ".bai"}


rule macs2:
    """⑧ 三条路线各 call 一次 peak：全部片段 / NFR / 单核小体。

    多核小体档不参与：跨多个核小体的片段两端落在不同开放区，peak 边界会被拉糊。
    这一档的用途是看核小体周期性（⑦⑭），不是定位开放区。

    -f BAMPE 用 reads 的真实跨度定义片段；此时 MACS2 会把 --shift/--extsize
    强制归零（OptValidator.py），传了不报错但静默无效，所以不传。
    """
    input:
        unpack(frac_inputs),
    output:
        peaks = res(D_PEAK, "{sample}.{frac}_peaks.narrowPeak"),
        xls = res(D_PEAK, "{sample}.{frac}_peaks.xls"),
        summits = res(D_PEAK, "{sample}.{frac}_summits.bed"),
        metrics = res(D_METRIC, "08_macs2.{sample}.{frac}.tsv"),
    params:
        genome = config["macs2"]["genome_size"],
        qvalue = config["macs2"]["qvalue"],
        extra = config["macs2"]["extra"],
        name = "{sample}.{frac}",
        outdir = res(D_PEAK),
    wildcard_constraints:
        frac = "all|nfr|mono",
    log:
        os.path.join(LOGS, "{sample}.{frac}.macs2.log"),
    benchmark:
        bench("macs2", "{sample}.{frac}"),
    threads: tcap(4)
    resources:
        mem_mb = 8000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        macs2 callpeak -t {input.bam} -f BAMPE {params.extra} \
            -g {params.genome} -q {params.qvalue} \
            -n {params.name} --outdir {params.outdir} 2> {log}

        n=$(wc -l < {output.peaks})
        printf '08_macs2\tpeaks_%s\t%s\tpeaks\t20k-120k\n' {wildcards.frac} "$n" > {output.metrics}
        """


rule blacklist:
    """⑨ 去掉 ENCODE 黑名单区域。

    着丝粒、端粒、rDNA 这些高拷贝区会产生大量堆叠的假信号。
    不去掉，"peak 数很多"只是假象，FRiP 也会被这些区域抬高。
    """
    input:
        peaks = res(D_PEAK, "{sample}.{frac}_peaks.narrowPeak"),
        blacklist = blacklist_bed,
    output:
        peaks = res(D_PEAK, "{sample}.{frac}.filtered.narrowPeak"),
        stats = res(D_PEAK, "{sample}.{frac}.blacklist_stats.tsv"),
        metrics = res(D_METRIC, "09_blacklist.{sample}.{frac}.tsv"),
    wildcard_constraints:
        frac = "all|nfr|mono",
    log:
        os.path.join(LOGS, "{sample}.{frac}.blacklist.log"),
    benchmark:
        bench("blacklist", "{sample}.{frac}"),
    threads: 1
    resources:
        mem_mb = 2000,
    conda:
        "../envs/atac.yaml"
    shell:
        """
        bedtools intersect -v -a {input.peaks} -b {input.blacklist} > {output.peaks}

        before=$(wc -l < {input.peaks})
        after=$(wc -l < {output.peaks})
        removed=$(awk -v a="$before" -v b="$after" 'BEGIN{{printf "%.2f", (a>0)? (a-b)*100/a : 0}}')

        printf 'before\tkept\tremoved_pct\n%s\t%s\t%s\n' "$before" "$after" "$removed" > {output.stats}

        {{
          printf '09_blacklist\tpeaks_%s_before_blacklist\t%s\tpeaks\t\n' {wildcards.frac} "$before"
          printf '09_blacklist\tpeaks_%s_after_blacklist\t%s\tpeaks\t\n' {wildcards.frac} "$after"
          printf '09_blacklist\tblacklist_removed_%s\t%s\t%%\t<1\n' {wildcards.frac} "$removed"
        }} > {output.metrics}
        """
