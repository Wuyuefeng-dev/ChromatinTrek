# rules/peaks.smk  – MACS3 peak calling with per-assay flags

_MACS3_FLAGS = {
    "atac":   "--nomodel --shift -100 --extsize 200",
    "chip":   "",
    "cuttag": "--nomodel --shift -75 --extsize 150",
}

rule macs3_callpeak:
    input:
        bam=f"{ALIGN_DIR}/{{sample}}/{{sample}}.filtered.bam",
    output:
        peak=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.narrowPeak",
        xls=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.xls",
    params:
        outdir=f"{PEAK_DIR}/{{sample}}",
        name="{sample}",
        gsize=GCFG.macs3_gsize,
        qval=CFG.get("qvalue", 0.05),
        flags=lambda wc: _MACS3_FLAGS.get(
            SAMPLES_DF.loc[wc.sample,"assay"].lower(), "--nomodel"),
    log: f"{LOG_DIR}/{{sample}}/macs3.log"
    resources:
        mem_mb=16000, cpus=4, runtime="02:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "macs3 callpeak -t {input.bam} -f BAMPE "
        "-g {params.gsize} -n {params.name} "
        "--outdir {params.outdir} -q {params.qval} "
        "{params.flags} &> {log}"


rule filter_peaks_blacklist:
    input:
        peak=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.narrowPeak",
    output:
        peak=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.filtered.narrowPeak",
    params:
        blacklist=GCFG.blacklist,
    log: f"{LOG_DIR}/{{sample}}/peak_filter.log"
    resources:
        mem_mb=4000, cpus=1, runtime="00:20:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "bedtools intersect -v -a {input.peak} -b {params.blacklist} "
        "> {output.peak} 2>{log}"
