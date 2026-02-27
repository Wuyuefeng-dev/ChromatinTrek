#!/usr/bin/env python3
"""
demo/generate_demo_figures.py
=============================
Generates synthetic publication-ready plots demonstrating the output
of ChromatinTrek (without requiring Snakemake execution or HOMER tools).

These plots will be saved to `demo/figures/` and embedded into the README.

Outputs generated:
1. Genomic annotation stacked bar chart
2. Distance-to-TSS histogram
3. TF Enrichment dot plot
4. GO Enrichment bar plot
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Add repo root to path to import chromatintrek modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from chromatintrek.peak_gene import plot_genomic_annotation, plot_tss_distance
from chromatintrek.motif_tf import plot_tf_dotplot
from chromatintrek.go_enrichment import plot_go_barplot

OUT_DIR = Path(__file__).parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _genomic_features() -> None:
    # Synthetic annotation fractions
    def gen_df(prom, intron, intergenic):
        n = 1000
        cats = (
            ["Promoter"] * int(n * prom) +
            ["Intron"] * int(n * intron) +
            ["Intergenic"] * int(n * intergenic) +
            ["Exon"] * int(n * (1.0 - prom - intron - intergenic))
        )
        return pd.DataFrame({"annotation_class": cats})

    dfs = {
        "ATAC-seq": gen_df(0.40, 0.30, 0.20),
        "ChIP-seq (TF)": gen_df(0.60, 0.15, 0.15),
        "CUT&Tag": gen_df(0.50, 0.25, 0.15),
    }

    plot_genomic_annotation(
        dfs, out_png=OUT_DIR / "demo_genomic_annotation.png",
        figsize=(8, 5)
    )


def _tss_distance() -> None:
    # Synthetic distance distribution (heavy cluster at 0)
    rng = np.random.default_rng(42)
    dist_tss = np.concatenate([
        rng.normal(0, 500, 4000),      # core promoter
        rng.normal(0, 5000, 2000),     # proximal
        rng.uniform(-50000, 50000, 1000) # distal
    ])
    df = pd.DataFrame({"distance_to_tss": dist_tss})
    plot_tss_distance(
        df,
        out_png=OUT_DIR / "demo_tss_distance.png",
        title="ATAC-seq Distance to Nearest TSS"
    )


def _tf_enrichment() -> None:
    # Synthetic HOMER knownResults parsing output
    data = {
        "tf_name": ["CTCF", "Fosl2", "JunB", "Fra1", "BATF", "AP-1", "TEAD", "RUNX1"],
        "pvalue": [1e-150, 1e-85, 1e-70, 1e-65, 1e-45, 1e-40, 1e-25, 1e-15],
        "enrichment": [8.5, 6.2, 5.8, 5.1, 4.0, 3.8, 2.5, 2.1],
        "pct_targets_val": [45.2, 32.1, 28.5, 25.0, 18.2, 16.5, 12.0, 9.5],
    }
    df = pd.DataFrame(data)
    df["log_pvalue"] = -np.log10(df["pvalue"])

    plot_tf_dotplot(
        df,
        out_png=OUT_DIR / "demo_tf_enrichment.png",
        title="Top Enriched Motifs (ATAC-seq peaks)"
    )


def _go_enrichment() -> None:
    # Synthetic gseapy output
    data = {
        "Category": [
            "GO_Biological_Process_2023", "GO_Biological_Process_2023", "GO_Biological_Process_2023",
            "GO_Molecular_Function_2023", "GO_Molecular_Function_2023",
            "KEGG_2021_Human", "KEGG_2021_Human"
        ],
        "term": [
            "regulation of transcription (GO:0006355)",
            "chromatin organization (GO:0006325)",
            "cell cycle (GO:0007049)",
            "DNA-binding transcription factor activity (GO:0003700)",
            "chromatin binding (GO:0003682)",
            "Pathways in cancer",
            "Cell cycle"
        ],
        "padj": [1e-25, 1e-18, 1e-12, 1e-20, 1e-15, 1e-10, 1e-8],
    }
    df = pd.DataFrame(data)

    plot_go_barplot(
        df, out_png=OUT_DIR / "demo_go_enrichment.png",
        title="GO/KEGG Enrichment of Target Genes"
    )


def main():
    print("Generating demo figures...")
    _genomic_features()
    _tss_distance()
    _tf_enrichment()
    _go_enrichment()
    print(f"Figures saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
