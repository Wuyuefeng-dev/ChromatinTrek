"""
chromatintrek.peak_gene
=======================
Peak-to-gene annotation, TSS-distance statistics, and genomic feature plots.

Uses pybedtools / PyRanges for fast interval intersection and a GTF file
for gene coordinate information.

Key outputs
-----------
* Annotated peak table with nearest gene, Ensembl ID, distance-to-TSS
* Genomic feature distribution (promoter, exon, intron, intergenic, …)
* TSS-distance violin / bar plot
* Peak count breakdown bar chart
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .utils import makedirs, require_tools, run_cmd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HOMER-based annotation (requires HOMER)
# ---------------------------------------------------------------------------

def annotate_peaks_homer(
    peak_bed: Union[str, Path],
    out_txt: Union[str, Path],
    genome: str = "hg38",
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """
    Annotate peaks with HOMER annotatePeaks.pl and return a DataFrame.

    Columns include: PeakID, chr, start, end, Annotation,
                     Distance to TSS, Nearest PromoterID, Gene Name, Gene ID
    """
    require_tools("annotatePeaks.pl")
    makedirs(Path(out_txt).parent)
    cmd = f"annotatePeaks.pl {peak_bed} {genome} > {out_txt}"
    log.info("annotatePeaks: %s -> %s", peak_bed, out_txt)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    df = pd.read_csv(out_txt, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    return df


# ---------------------------------------------------------------------------
# Bedtools-based annotation (no HOMER required)
# ---------------------------------------------------------------------------

def annotate_peaks_bedtools(
    peak_bed: Union[str, Path],
    tss_bed: Union[str, Path],
    out_txt: Union[str, Path],
    distance_cutoffs: Tuple[int, ...] = (500, 1000, 5000, 10000, 50000),
    dry_run: bool = False,
) -> pd.DataFrame:
    """
    Annotate peaks with the nearest gene TSS using bedtools closest.

    Parameters
    ----------
    peak_bed         : Input narrowPeak / BED file.
    tss_bed          : 6-column BED of TSS positions with gene name in col4.
    distance_cutoffs : Distance bins (bp) for summary statistics.

    Returns
    -------
    DataFrame with columns: peak_chr, peak_start, peak_end, gene_name,
                             distance_to_tss, annotation_class
    """
    require_tools("bedtools")
    makedirs(Path(out_txt).parent)
    tmp = str(out_txt) + ".tmp"
    cmd = f"bedtools closest -a {peak_bed} -b {tss_bed} -D ref > {tmp}"
    run_cmd(cmd, dry_run=dry_run)
    if dry_run:
        return pd.DataFrame()

    df = pd.read_csv(tmp, sep="\t", header=None)
    df = df.iloc[:, [0, 1, 2, 9, -1]]
    df.columns = ["chr", "start", "end", "gene_name", "distance_to_tss"]

    # Annotate by distance
    def _classify(d: float) -> str:
        ad = abs(d)
        if ad <= 1000:
            return "Promoter (≤1kb)"
        elif ad <= 5000:
            return "Promoter (1-5kb)"
        elif ad <= 10000:
            return "Distal (5-10kb)"
        elif ad <= 50000:
            return "Distal (10-50kb)"
        else:
            return "Intergenic (>50kb)"

    df["annotation_class"] = df["distance_to_tss"].apply(_classify)
    df.to_csv(out_txt, sep="\t", index=False)
    return df


# ---------------------------------------------------------------------------
# Genomic feature annotation from HOMER output
# ---------------------------------------------------------------------------

_ANNO_CATEGORIES = {
    "promoter": "Promoter",
    "tss": "Promoter",
    "5'utr": "5' UTR",
    "3'utr": "3' UTR",
    "exon": "Exon",
    "intron": "Intron",
    "intergenic": "Intergenic",
    "non-coding": "Non-coding",
}


def parse_annotation_column(annotation_series: pd.Series) -> pd.Series:
    """Map HOMER Annotation column strings to simplified feature categories."""
    def _map(val: str) -> str:
        if not isinstance(val, str):
            return "Other"
        v = val.lower()
        for key, label in _ANNO_CATEGORIES.items():
            if key in v:
                return label
        return "Other"
    return annotation_series.map(_map)


# ---------------------------------------------------------------------------
# Plotting: genomic feature distribution
# ---------------------------------------------------------------------------

FEATURE_COLORS = {
    "Promoter": "#E64B35",
    "5' UTR": "#F39B7F",
    "3' UTR": "#FCCF77",
    "Exon": "#4DBBD5",
    "Intron": "#91D1C2",
    "Intergenic": "#8491B4",
    "Non-coding": "#B09C85",
    "Other": "#ADB4BC",
}


def plot_genomic_annotation(
    annot_dfs: Dict[str, pd.DataFrame],
    annotation_col: str = "annotation_class",
    out_png: Union[str, Path] = "genomic_annotation.png",
    figsize: Tuple[float, float] = (14, 5),
    dpi: int = 150,
) -> Path:
    """
    Stacked bar chart of genomic feature distribution for multiple assays.

    Parameters
    ----------
    annot_dfs : dict {assay_label: DataFrame with *annotation_col*}
    """
    makedirs(Path(out_png).parent)
    categories = list(FEATURE_COLORS.keys())
    assays = list(annot_dfs.keys())

    # Compute fractions
    fracs = {}
    for label, df in annot_dfs.items():
        vc = df[annotation_col].value_counts(normalize=True)
        fracs[label] = {cat: vc.get(cat, 0.0) for cat in categories}

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    bottoms = np.zeros(len(assays))
    for cat in categories:
        vals = [fracs[a][cat] * 100 for a in assays]
        ax.bar(assays, vals, bottom=bottoms,
               color=FEATURE_COLORS[cat], label=cat, edgecolor="white", linewidth=0.4)
        bottoms += np.array(vals)

    ax.set_ylabel("Fraction of Peaks (%)", fontsize=11)
    ax.set_ylim(0, 105)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.8,
              bbox_to_anchor=(1.18, 1))
    ax.set_title("Genomic Feature Distribution", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Genomic annotation plot: %s", out_png)
    return Path(out_png)


# ---------------------------------------------------------------------------
# Plotting: TSS distance distribution
# ---------------------------------------------------------------------------

def plot_tss_distance(
    annot_df: pd.DataFrame,
    distance_col: str = "distance_to_tss",
    out_png: Union[str, Path] = "tss_distance.png",
    title: str = "Distance to Nearest TSS",
    color: str = "#4C72B0",
    dpi: int = 150,
) -> Path:
    """Histogram of peak distances to the nearest TSS."""
    makedirs(Path(out_png).parent)
    distances = annot_df[distance_col].abs().clip(upper=100_000)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=dpi)
    ax.hist(distances, bins=100, color=color, edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.set_xlabel("Distance to TSS (bp)", fontsize=11)
    ax.set_ylabel("Number of Peaks", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_png)
