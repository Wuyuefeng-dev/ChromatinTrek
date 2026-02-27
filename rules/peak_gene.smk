# rules/peak_gene.smk – HOMER annotatePeaks for gene association

rule annotate_peaks_gene:
    input:
        peak=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.filtered.narrowPeak",
    output:
        tsv=f"{GENE_DIR}/{{sample}}_annotated.tsv",
    params:
        genome=GCFG.homer_genome,
    log: f"{LOG_DIR}/{{sample}}/annotate_peaks.log"
    resources:
        mem_mb=16000, cpus=4, runtime="01:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    shell:
        "annotatePeaks.pl {input.peak} {params.genome} > {output.tsv} 2>{log}"


rule peak_gene_stats:
    """Summarise genomic annotation fractions using Python."""
    input:
        tsv=f"{GENE_DIR}/{{sample}}_annotated.tsv",
    output:
        png=f"{GENE_DIR}/{{sample}}_annotation_plot.png",
    params:
        sample="{sample}",
    log: f"{LOG_DIR}/{{sample}}/peak_gene_stats.log"
    resources:
        mem_mb=8000, cpus=2, runtime="00:30:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    run:
        import pandas as pd
        from chromatintrek.peak_gene import parse_annotation_column, plot_genomic_annotation, plot_tss_distance
        df = pd.read_csv(input.tsv, sep="\t")
        df.columns = [c.strip() for c in df.columns]
        if "Annotation" in df.columns:
            df["annotation_class"] = parse_annotation_column(df["Annotation"])
        else:
            df["annotation_class"] = "Intergenic"
            
        # 1. ChIPseeker equivalent: plotAnnoBar
        plot_genomic_annotation(
            {params.sample: df},
            annotation_col="annotation_class",
            out_png=output.png,
        )
        # 2. ChIPseeker equivalent: plotDistToTSS
        dist_col = next((c for c in df.columns if "distance" in c.lower() and "tss" in c.lower()), None)
        if dist_col:
            plot_tss_distance(
                df, distance_col=dist_col,
                out_png=output.tss_png,
                title=f"{params.sample} Distance to TSS",
            )
        else:
            open(output.tss_png, "w").close()
