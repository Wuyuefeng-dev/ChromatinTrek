"""
chromatintrek.go_enrichment
============================
Gene Ontology and pathway enrichment analysis via gseapy.

Functions
---------
run_enrichment()     Run GO/KEGG enrichment for a gene list
plot_go_barplot()    Horizontal bar chart of top enriched terms
plot_go_dotplot()    Dot plot coloured by p-value, sized by gene count
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Gene set → organism code mapping
# ---------------------------------------------------------------------------

_ORGANISM_GENESETS: Dict[str, Dict[str, List[str]]] = {
    "Human": {
        "gene_sets": [
            "GO_Biological_Process_2023",
            "GO_Molecular_Function_2023",
            "GO_Cellular_Component_2023",
            "KEGG_2021_Human",
        ],
        "enrichr_organism": "Human",
    },
    "Mouse": {
        "gene_sets": [
            "GO_Biological_Process_2023",
            "GO_Molecular_Function_2023",
            "GO_Cellular_Component_2023",
            "KEGG_2021_Mouse",
        ],
        "enrichr_organism": "Mouse",
    },
}


# ---------------------------------------------------------------------------
# Enrichment
# ---------------------------------------------------------------------------

def run_enrichment(
    gene_list: List[str],
    organism: str = "Human",
    gene_sets: Optional[List[str]] = None,
    out_dir: Optional[Union[str, Path]] = None,
    cutoff_padj: float = 0.05,
    min_genes: int = 5,
) -> pd.DataFrame:
    """
    Run ORA (over-representation analysis) via gseapy.enrichr.

    Parameters
    ----------
    gene_list   : List of HGNC / MGI gene symbols.
    organism    : "Human" | "Mouse"  (auto-resolved from species config)
    gene_sets   : Gene set libraries to query (defaults to GO + KEGG).
    out_dir     : If provided, save result TSVs here.
    cutoff_padj : Adjusted p-value cutoff.
    min_genes   : Minimum overlap to retain a term.

    Returns
    -------
    DataFrame with columns: Gene_set, Term, P-value, Adjusted P-value,
                             Overlap, Genes, Category
    """
    try:
        import gseapy as gp
    except ImportError:
        raise ImportError("gseapy is required: pip install gseapy")

    org_info = _ORGANISM_GENESETS.get(organism, _ORGANISM_GENESETS["Human"])
    if gene_sets is None:
        gene_sets = org_info["gene_sets"]
    enrichr_org = org_info["enrichr_organism"]

    results = []
    for gs in gene_sets:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=gs,
                organism=enrichr_org,
                outdir=None,
                verbose=False,
            )
            df = enr.results.copy()
            df["Category"] = gs
            results.append(df)
        except Exception as e:
            log.warning("gseapy enrichr failed for %s: %s", gs, e)

    if not results:
        return pd.DataFrame()

    combined = pd.concat(results, ignore_index=True)
    # Standardise column names across gseapy versions
    col_map = {
        "Adjusted P-value": "padj",
        "P-value": "pvalue",
        "Term": "term",
        "Overlap": "overlap",
        "Genes": "genes",
    }
    combined = combined.rename(columns=col_map)

    if "padj" in combined.columns:
        combined = combined[combined["padj"] <= cutoff_padj]

    # Filter by minimum gene overlap
    def _overlap_count(s: str) -> int:
        try:
            num, _ = s.split("/")
            return int(num)
        except Exception:
            return 0

    if "overlap" in combined.columns:
        combined["overlap_count"] = combined["overlap"].apply(_overlap_count)
        combined = combined[combined["overlap_count"] >= min_genes]

    combined = combined.sort_values("padj").reset_index(drop=True)

    if out_dir:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        combined.to_csv(Path(out_dir) / "go_enrichment.tsv", sep="\t", index=False)
        log.info("GO results saved to %s", out_dir)

    return combined


# ---------------------------------------------------------------------------
# Bar plot
# ---------------------------------------------------------------------------

_CATEGORY_COLORS = {
    "GO_Biological_Process_2023": "#4DBBD5",
    "GO_Molecular_Function_2023": "#E64B35",
    "GO_Cellular_Component_2023": "#00A087",
    "KEGG_2021_Human":            "#3C5488",
    "KEGG_2021_Mouse":            "#3C5488",
}


def plot_go_barplot(
    df: pd.DataFrame,
    out_png: Union[str, Path],
    top_n: int = 15,
    title: str = "GO Enrichment",
    figsize: Tuple[float, float] = (8, 6),
    dpi: int = 150,
) -> Path:
    """
    Horizontal bar chart: x = -log10(padj), coloured by gene set category.
    """
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    plot_df = df.head(top_n).copy()
    plot_df["-log10(padj)"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
    plot_df = plot_df.sort_values("-log10(padj)", ascending=True)

    colors = [_CATEGORY_COLORS.get(c, "#8491B4") for c in plot_df["Category"]]
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    bars = ax.barh(plot_df["term"], plot_df["-log10(padj)"],
                   color=colors, edgecolor="white", linewidth=0.3, height=0.7)
    ax.set_xlabel("-log₁₀(adjusted p-value)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.axvline(x=-np.log10(0.05), color="grey", linestyle="--",
               linewidth=0.8, alpha=0.7, label="padj = 0.05")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, loc="lower right")

    # Legend for categories
    seen = {}
    handles = []
    for cat, col in _CATEGORY_COLORS.items():
        short = cat.split("_")[1].replace("Process", "BP").replace("Function", "MF").replace("Component", "CC")
        if cat in plot_df["Category"].values and cat not in seen:
            handles.append(mpatches.Patch(color=col, label=short))
            seen[cat] = True
    if handles:
        ax.legend(handles=handles, fontsize=8, loc="lower right", title="Gene Set")

    plt.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("GO barplot: %s", out_png)
    return Path(out_png)


# ---------------------------------------------------------------------------
# Dot plot
# ---------------------------------------------------------------------------

def plot_go_dotplot(
    df: pd.DataFrame,
    out_png: Union[str, Path],
    top_n: int = 20,
    title: str = "GO Enrichment",
    figsize: Tuple[float, float] = (7, 8),
    dpi: int = 150,
) -> Path:
    """
    Dot plot: x = gene ratio, y = term, size = gene count, colour = -log10(padj).
    """
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    plot_df = df.head(top_n).copy()
    plot_df["-log10(padj)"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
    if "overlap" in plot_df.columns:
        plot_df["overlap_count"] = plot_df["overlap"].apply(
            lambda s: int(str(s).split("/")[0]) if "/" in str(s) else 5)
    else:
        plot_df["overlap_count"] = 10

    plot_df = plot_df.sort_values("-log10(padj)", ascending=True)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    sc = ax.scatter(
        plot_df["-log10(padj)"],
        range(len(plot_df)),
        s=plot_df["overlap_count"] * 15,
        c=plot_df["-log10(padj)"],
        cmap="YlOrRd",
        alpha=0.85,
        edgecolors="grey",
        linewidths=0.3,
    )
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["term"], fontsize=8)
    ax.set_xlabel("-log₁₀(adjusted p-value)", fontsize=10)
    ax.set_title(title, fontsize=13, fontweight="bold")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label("-log₁₀(padj)", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("GO dotplot: %s", out_png)
    return Path(out_png)
