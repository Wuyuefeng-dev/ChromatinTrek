# rules/align.smk  – Bowtie2, samtools dedup/filter/flagstat

# Assay-specific Bowtie2 flags
_BT2_FLAGS = {
    "atac":   "--no-mixed --no-discordant -X 2000",
    "chip":   "--no-mixed --no-discordant",
    "cuttag": "--no-mixed --no-discordant -X 700",
}

rule bowtie2_align:
    input:
        r1=f"{TRIM_DIR}/{{sample}}/{{sample}}_R1_trimmed.fq.gz",
        r2=f"{TRIM_DIR}/{{sample}}/{{sample}}_R2_trimmed.fq.gz",
    output:
        bam=f"{ALIGN_DIR}/{{sample}}/{{sample}}.raw.bam",
        bai=f"{ALIGN_DIR}/{{sample}}/{{sample}}.raw.bam.bai",
    params:
        index=GCFG.bowtie2_index,
        flags=lambda wc: _BT2_FLAGS.get(
            SAMPLES_DF.loc[wc.sample,"assay"].lower(), "--no-mixed --no-discordant"),
        threads=16,
    log: f"{LOG_DIR}/{{sample}}/bowtie2.log"
    resources:
        mem_mb=32000, cpus=16, runtime="04:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "bowtie2 -p {params.threads} -x {params.index} {params.flags} "
        "-1 {input.r1} -2 {input.r2} 2>{log} "
        "| samtools sort -@ {params.threads} -o {output.bam} - && "
        "samtools index {output.bam}"


rule markdup:
    input:
        bam=f"{ALIGN_DIR}/{{sample}}/{{sample}}.raw.bam",
    output:
        bam=f"{ALIGN_DIR}/{{sample}}/{{sample}}.dedup.bam",
        bai=f"{ALIGN_DIR}/{{sample}}/{{sample}}.dedup.bam.bai",
    params:
        threads=8,
    log: f"{LOG_DIR}/{{sample}}/markdup.log"
    resources:
        mem_mb=16000, cpus=8, runtime="02:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "samtools sort -n -@ {params.threads} {input.bam} "
        "| samtools fixmate -m - - "
        "| samtools sort -@ {params.threads} - "
        "| samtools markdup -r -@ {params.threads} - {output.bam} 2>{log} && "
        "samtools index {output.bam}"


rule filter_bam:
    input:
        bam=f"{ALIGN_DIR}/{{sample}}/{{sample}}.dedup.bam",
    output:
        bam=f"{ALIGN_DIR}/{{sample}}/{{sample}}.filtered.bam",
        bai=f"{ALIGN_DIR}/{{sample}}/{{sample}}.filtered.bam.bai",
    params:
        blacklist=GCFG.blacklist,
        mapq=CFG.get("mapq", 30),
        threads=8,
    log: f"{LOG_DIR}/{{sample}}/filter.log"
    resources:
        mem_mb=16000, cpus=8, runtime="01:30:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "samtools view -@ {params.threads} -q {params.mapq} -f 2 -b {input.bam} "
        "| bedtools intersect -v -abam stdin -b {params.blacklist} "
        "| samtools sort -@ {params.threads} -o {output.bam} - && "
        "samtools index {output.bam} 2>{log}"


rule flagstat:
    input:
        bam=f"{ALIGN_DIR}/{{sample}}/{{sample}}.filtered.bam",
    output:
        txt=f"{ALIGN_DIR}/{{sample}}/{{sample}}.flagstat.txt",
    log: f"{LOG_DIR}/{{sample}}/flagstat.log"
    resources:
        mem_mb=4000, cpus=2, runtime="00:20:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "samtools flagstat {input.bam} > {output.txt} 2>{log}"
