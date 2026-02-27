"""
chromatintrek.motif_tf
======================
Motif finding and TF enrichment using HOMER.

Steps
-----
1. findMotifsGenome.pl   – de-novo and known motif enrichment
2. annotatePeaks.pl      – annotate peaks with motif presence
3. parse_homer_known()   – parse knownResults.txt into a DataFrame
4. plot_tf_dotplot()     – bubble / dot plot of top enriched TFs

The parse and plot functions work on the HOMER output files and produce
publication-quality figures using matplotlib only (no tool required).
"""

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

from .utils import makedirs, require_tools, run_cmd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HOMER findMotifsGenome
# ---------------------------------------------------------------------------

def find_motifs(
    peak_bed: Union[str, Path],
    out_dir: Union[str, Path],
    genome: str = "hg38",
    size: int = 200,
    mask: bool = True,
    cpus: int = 16,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Run HOMER findMotifsGenome.pl on a peak BED file.

    Known motifs output: <out_dir>/knownResults.txt
    De-novo motifs:      <out_dir>/homerMotifs.all.motifs
    """
    require_tools("findMotifsGenome.pl")
    makedirs(out_dir)
    mask_flag = "-mask" if mask else ""
    cmd = (
        f"findMotifsGenome.pl {peak_bed} {genome} {out_dir} "
        f"-size {size} -p {cpus} {mask_flag}"
    )
    log.info("HOMER findMotifsGenome: %s -> %s", peak_bed, out_dir)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_dir)


# ---------------------------------------------------------------------------
# HOMER annotatePeaks
# ---------------------------------------------------------------------------

def annotate_peaks_homer(
    peak_bed: Union[str, Path],
    genome: str = "hg38",
    out_txt: Optional[Union[str, Path]] = None,
    motif_dir: Optional[Union[str, Path]] = None,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Run HOMER annotatePeaks.pl to annotate peaks with genome features and motifs.
    """
    require_tools("annotatePeaks.pl")
    if out_txt is None:
        out_txt = Path(peak_bed).with_suffix(".homer_annot.txt")
    makedirs(Path(out_txt).parent)
    motif_flag = f"-m {motif_dir}" if motif_dir else ""
    cmd = (
        f"annotatePeaks.pl {peak_bed} {genome} {motif_flag} > {out_txt}"
    )
    log.info("HOMER annotatePeaks: %s", peak_bed)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_txt)


# ---------------------------------------------------------------------------
# Parse HOMER knownResults.txt
# ---------------------------------------------------------------------------

def parse_homer_known(
    known_results_txt: Union[str, Path],
    top_n: int = 20,
    pvalue_cutoff: float = 0.05,
) -> pd.DataFrame:
    """
    Parse HOMER knownResults.txt into a tidy DataFrame.

    Columns returned:
        tf_name, pvalue, log_pvalue, enrichment, pct_targets, pct_background
    """
    df = pd.read_csv(known_results_txt, sep="\t", comment="#")
    df.columns = [c.strip() for c in df.columns]

    # HOMER column names vary slightly by version – harmonise
    col_map = {
        "Motif Name": "tf_name",
        "P-value": "pvalue",
        "Log P-value": "log_pvalue",
        "% of Target Sequences with Motif": "pct_targets",
        "% of Background Sequences with Motif": "pct_background",
    }
    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    # Strip trailing motif version numbers  e.g. "CTCF(Zf)/GM12878-CTCF-ChIP-Seq/Homer" -> "CTCF"
    if "tf_name" in df.columns:
        df["tf_name"] = df["tf_name"].str.split("(").str[0].str.strip()

    if "pvalue" in df.columns:
        df = df[df["pvalue"] <= pvalue_cutoff]
        df["log_pvalue"] = -np.log10(df["pvalue"].clip(lower=1e-300))

    # Enrichment ratio
    if "pct_targets" in df.columns and "pct_background" in df.columns:
        pct_t = df["pct_targets"].str.rstrip("%").astype(float)
        pct_b = df["pct_background"].str.rstrip("%").astype(float).clip(lower=0.001)
        df["enrichment"] = pct_t / pct_b
        df["pct_targets_val"] = pct_t

    return df.head(top_n).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plot TF enrichment dot plot
# ---------------------------------------------------------------------------

def plot_tf_dotplot(
    df: pd.DataFrame,
    out_png: Union[str, Path],
    title: str = "TF Enrichment",
    top_n: int = 20,
    figsize: Tuple[int, int] = (7, 8),
    dpi: int = 150,
) -> Path:
    """
    Bubble / dot plot of TF enrichment results.

    x-axis : enrichment ratio (% targets / % background)
    y-axis : TF name (sorted by enrichment)
    size   : % of target peaks with motif
    color  : -log10(p-value)
    """
    makedirs(Path(out_png).parent)
    df = df.head(top_n).copy()
    df = df.sort_values("enrichment", ascending=True)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    scatter = ax.scatter(
        df["enrichment"],
        df["tf_name"],
        s=df["pct_targets_val"] * 10,
        c=df["log_pvalue"],
        cmap="plasma",
        alpha=0.85,
        edgecolors="white",
        linewidths=0.3,
    )
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("-log₁₀(p-value)", fontsize=9)

    ax.set_xlabel("Enrichment (% targets / % background)", fontsize=10)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.axvline(1, color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("TF dotplot saved: %s", out_png)
    return Path(out_png)
