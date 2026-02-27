# rules/go.smk – Gene Ontology enrichment via gseapy

rule go_enrichment:
    """
    Extract gene list from HOMER annotation and run GO/KEGG enrichment.
    Produces a TSV and bar plot.
    """
    input:
        tsv=f"{GENE_DIR}/{{sample}}_annotated.tsv",
    output:
        tsv=f"{GO_DIR}/{{sample}}/go_enrichment.tsv",
        barplot=f"{GO_DIR}/{{sample}}/go_barplot.png",
        dotplot=f"{GO_DIR}/{{sample}}/go_dotplot.png",
    params:
        outdir=f"{GO_DIR}/{{sample}}",
        organism=GCFG.go_organism,
        cutoff=CFG.get("qvalue", 0.05),
    log: f"{LOG_DIR}/{{sample}}/go_enrichment.log"
    resources:
        mem_mb=8000, cpus=4, runtime="01:00:00",
        slurm_partition=lambda wc: CFG.get("slurm_partition","short"),
        slurm_account=lambda wc: CFG.get("slurm_account",""),
    run:
        import pandas as pd
        from chromatintrek.go_enrichment import run_enrichment, plot_go_barplot, plot_go_dotplot

        df = pd.read_csv(input.tsv, sep="\t")
        df.columns = [c.strip() for c in df.columns]

        gene_col = next(
            (c for c in df.columns if "gene name" in c.lower() or c.lower() == "gene_name"),
            None,
        )
        genes = []
        if gene_col:
            genes = df[gene_col].dropna().astype(str).unique().tolist()
            genes = [g for g in genes if g not in ("", "NA", "nan")]

        if len(genes) < 5:
            # No genes to test – write empty outputs
            pd.DataFrame().to_csv(output.tsv, sep="\t", index=False)
            open(output.barplot, "w").close()
            open(output.dotplot, "w").close()
        else:
            result = run_enrichment(
                genes, organism=params.organism,
                out_dir=params.outdir, cutoff_padj=params.cutoff,
            )
            if not result.empty:
                plot_go_barplot(result, output.barplot,
                                title=f"{wildcards.sample} GO Enrichment")
                plot_go_dotplot(result, output.dotplot,
                                title=f"{wildcards.sample} GO Enrichment")
