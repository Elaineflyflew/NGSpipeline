################################################################################
# 00_common.smk —— 全局常量、路径、小工具函数
#
# 对应 snakePipes 的 shared/rules/ 里那些被到处 include 的公共片段。
################################################################################

import os
from shutil import which

SAMPLES = sorted(config["samples"].keys())
FRACTIONS = ["all", "nfr", "mono"]

THREADS = int(config.get("threads", 8))
MAPQ = int(config.get("mapq", 20))

RES = "results"
LOGS = "logs"
BENCH = os.path.join(RES, ".benchmark")

D_INPUT  = "00_input"
D_QC     = "01_fastqc"
D_TRIM   = "02_fastp"
D_ALN    = "03_map"
D_BAM    = "04_bam"
D_FRAG   = "05_fragbin"
D_PEAK   = "06_peaks"
D_FRIP   = "07_frip"
D_COV    = "08_bigwig"
D_TSS    = "09_tss"
D_ANN    = "10_annot"
D_MOTIF  = "11_motif"
D_ATACQC = "12_atacqc"
D_IGV    = "13_igv"
D_METRIC = "metrics"
D_REPORT = "99_report"

SCRIPTS = os.path.join(PROJECT_ROOT, "src", "py")


def res(*parts):
    return os.path.join(RES, *parts)


def bench(rule, sample="all"):
    return os.path.join(BENCH, "%s.%s.txt" % (rule, sample))


def logf(name):
    # log 必须和 output 带同一组通配符，否则 Snakemake 会拒绝（防止两个 job 写同一日志）
    return os.path.join(LOGS, "{sample}." + name + ".log")


def tcap(n):
    """单 rule 线程封顶。"""
    return max(1, min(int(n), THREADS))


def contig_filter():
    """把 drop_contigs 翻成 samtools -e 表达式，支持 * 通配。"""
    parts = []
    for c in config.get("drop_contigs", []):
        c = str(c).strip()
        if not c:
            continue
        if "*" in c:
            parts.append('rname !~ "%s"' % c.replace(".", r"\.").replace("*", ".*"))
        else:
            parts.append('rname != "%s"' % c)
    return " && ".join(parts) if parts else "true"


def frac_bam(wildcards):
    """⑦ 的三份 BAM 分别喂给 ⑧ 的三条 peak calling 路线。"""
    return {
        "all":  res(D_BAM, "%s.srt2.bam" % wildcards.sample),
        "nfr":  res(D_BAM, "%s_nucleosome_free.bam" % wildcards.sample),
        "mono": res(D_BAM, "%s_mono_nucleosome.bam" % wildcards.sample),
    }[wildcards.frac]


def homer_available():
    return which("annotatePeaks.pl") is not None and which("findMotifsGenome.pl") is not None


def annot_targets():
    """⑫⑬ 依赖 HOMER；没装就把它们从 DAG 里摘掉，而不是让整条链失败。"""
    if not homer_available():
        return []
    t = expand(res(D_ANN, "{sample}", "annotation_distribution.tsv"), sample=SAMPLES)
    t += expand(res(D_MOTIF, "{sample}", "{frac}", "knownResults.txt"),
                sample=SAMPLES, frac=FRACTIONS)
    return t


def annot_metric_deps():
    """给 collect_metrics 挂依赖用：这些产物在，说明对应指标文件也已经写好。"""
    if not homer_available():
        return []
    return [expand(res(D_ANN, "{sample}", "annotation_distribution.tsv"), sample=SAMPLES),
            expand(res(D_MOTIF, "{sample}", "{frac}", "top_motifs.tsv"),
                   sample=SAMPLES, frac=FRACTIONS)]


def atacqc_targets():
    if not config.get("run_atacqc", True):
        return []
    return expand(res(D_ATACQC, "{sample}.atacqc_metrics.tsv"), sample=SAMPLES)


def igv_targets():
    if not config.get("igv", {}).get("enabled", True):
        return []
    return expand(res(D_IGV, "{sample}.MYC.png"), sample=SAMPLES)


def multiqc_targets():
    if not config.get("multiqc", True):
        return []
    return [res(D_REPORT, "multiqc_report.html")]


################################################################################
# 原始 FASTQ 统一命名
#
# FastQC 的输出名是从输入文件名推出来的。原始名带 "sub20M.SRR13849313_1" 这种
# 前缀，下游没法按确定的名字取产物。先软链成 <sample>_R1.fastq.gz，
# 后面所有 rule 就都能用确定的名字了。
################################################################################
rule stage_input:
    output:
        r1 = res(D_INPUT, "{sample}_R1.fastq.gz"),
        r2 = res(D_INPUT, "{sample}_R2.fastq.gz"),
    params:
        src1 = lambda w: config["samples"][w.sample]["r1"],
        src2 = lambda w: config["samples"][w.sample]["r2"],
    log:
        logf("stage_input")
    run:
        for src, dst in ((params.src1, output.r1), (params.src2, output.r2)):
            if not os.path.exists(src):
                raise FileNotFoundError("输入不存在: %s" % src)
            if not os.path.exists(dst):
                os.symlink(os.path.abspath(src), dst)
