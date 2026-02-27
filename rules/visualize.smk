# rules/visualize.smk – deepTools: bamCoverage, computeMatrix, plotHeatmap, plotProfile

rule bamcoverage:
    input:
        bam=f"{ALIGN_DIR}/{{sample}}/{{sample}}.filtered.bam",
        bai=f"{ALIGN_DIR}/{{sample}}/{{sample}}.filtered.bam.bai",
    output:
        bw=f"{BW_DIR}/{{sample}}.bw",
    params:
        eff_gsize=GCFG.effective_genome_size,
        blacklist=GCFG.blacklist,
        threads=16,
    log: f"{LOG_DIR}/{{sample}}/bamcoverage.log"
    resources:
        mem_mb=16000, cpus=16, runtime="02:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "bamCoverage -b {input.bam} -o {output.bw} "
        "--binSize 10 --normalizeUsing RPGC "
        "--effectiveGenomeSize {params.eff_gsize} "
        "--extendReads --blackListFileName {params.blacklist} "
        "-p {params.threads} &> {log}"


rule compute_matrix:
    input:
        bw=f"{BW_DIR}/{{sample}}.bw",
        peaks=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.filtered.narrowPeak",
    output:
        matrix=f"{MATRIX_DIR}/{{sample}}_matrix.gz",
    params:
        upstream=3000,
        downstream=3000,
        threads=16,
    log: f"{LOG_DIR}/{{sample}}/computematrix.log"
    resources:
        mem_mb=32000, cpus=16, runtime="03:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "computeMatrix reference-point --referencePoint center "
        "-b {params.upstream} -a {params.downstream} "
        "-S {input.bw} -R {input.peaks} "
        "--binSize 50 -p {params.threads} "
        "-o {output.matrix} &> {log}"


rule plot_heatmap:
    input:
        matrix=f"{MATRIX_DIR}/{{sample}}_matrix.gz",
    output:
        png=f"{MATRIX_DIR}/{{sample}}_heatmap.png",
    params:
        title="{sample} Signal at Peaks",
    log: f"{LOG_DIR}/{{sample}}/plotheatmap.log"
    resources:
        mem_mb=8000, cpus=4, runtime="00:30:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "plotHeatmap -m {input.matrix} -out {output.png} "
        "--colorMap Blues --sortRegions descend "
        '--plotTitle "{params.title}" --dpi 300 &> {log}'


rule plot_profile:
    input:
        matrix=f"{MATRIX_DIR}/{{sample}}_matrix.gz",
    output:
        png=f"{MATRIX_DIR}/{{sample}}_profile.png",
    log: f"{LOG_DIR}/{{sample}}/plotprofile.log"
    resources:
        mem_mb=4000, cpus=2, runtime="00:20:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "plotProfile -m {input.matrix} -out {output.png} --dpi 300 &> {log}"
