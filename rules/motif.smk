# rules/motif.smk – HOMER findMotifsGenome + annotatePeaks

rule homer_find_motifs:
    input:
        peak=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.filtered.narrowPeak",
    output:
        known=f"{MOTIF_DIR}/{{sample}}/knownResults.txt",
    params:
        outdir=f"{MOTIF_DIR}/{{sample}}",
        genome=GCFG.homer_genome,
        size=200,
        cpus=16,
    log: f"{LOG_DIR}/{{sample}}/homer_motifs.log"
    resources:
        mem_mb=32000, cpus=16, runtime="04:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "findMotifsGenome.pl {input.peak} {params.genome} {params.outdir} "
        "-size {params.size} -p {params.cpus} -mask &> {log}"


rule homer_annotate_peaks:
    input:
        peak=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.filtered.narrowPeak",
        motif_dir=f"{MOTIF_DIR}/{{sample}}",
    output:
        txt=f"{MOTIF_DIR}/{{sample}}/{{sample}}_homer_annot.txt",
    params:
        genome=GCFG.homer_genome,
    log: f"{LOG_DIR}/{{sample}}/homer_annot.log"
    resources:
        mem_mb=16000, cpus=4, runtime="01:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "annotatePeaks.pl {input.peak} {params.genome} "
        "-m {input.motif_dir} > {output.txt} 2>{log}"
