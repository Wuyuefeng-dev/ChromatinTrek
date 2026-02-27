# rules/qc.smk  – FastQC, Trim Galore, MultiQC

rule fastqc_raw:
    input:
        r1=lambda wc: SAMPLES_DF.loc[wc.sample, "fastq_r1"],
    output:
        html=f"{QC_DIR}/{{sample}}/pre_trim/{{sample}}_R1_fastqc.html",
    params:
        outdir=f"{QC_DIR}/{{sample}}/pre_trim",
        threads=4,
    log:    f"{LOG_DIR}/{{sample}}/fastqc_raw.log"
    resources:
        mem_mb=8000, cpus=4, runtime="01:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "fastqc --threads {params.threads} --outdir {params.outdir} "
        "{input.r1} &> {log}"


rule trim_galore:
    input:
        r1=lambda wc: SAMPLES_DF.loc[wc.sample, "fastq_r1"],
        r2=lambda wc: SAMPLES_DF.loc[wc.sample, "fastq_r2"],
    output:
        r1=f"{TRIM_DIR}/{{sample}}/{{sample}}_R1_trimmed.fq.gz",
        r2=f"{TRIM_DIR}/{{sample}}/{{sample}}_R2_trimmed.fq.gz",
    params:
        outdir=f"{TRIM_DIR}/{{sample}}",
        extra=CFG.get("trim_extra",""),
    log:    f"{LOG_DIR}/{{sample}}/trim_galore.log"
    resources:
        mem_mb=8000, cpus=4, runtime="01:30:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "trim_galore --paired --quality 20 --length 25 --cores 4 "
        "--fastqc --output_dir {params.outdir} {params.extra} "
        "{input.r1} {input.r2} &> {log}"


rule multiqc:
    input:
        expand(f"{QC_DIR}/{{sample}}/pre_trim/{{sample}}_R1_fastqc.html", sample=SAMPLES),
        expand(f"{TRIM_DIR}/{{sample}}/{{sample}}_R1_trimmed.fq.gz", sample=SAMPLES),
    output:
        f"{QC_DIR}/multiqc_report.html",
    params:
        search=f"{QC_DIR} {TRIM_DIR}",
        outdir=QC_DIR,
    log:    f"{LOG_DIR}/multiqc.log"
    resources:
        mem_mb=4000, cpus=1, runtime="00:30:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "multiqc --outdir {params.outdir} --force {params.search} &> {log}"
