"""
chromatintrek.peak_calling
==========================
Peak calling via MACS3 with per-assay parameter presets.

Assay presets
-------------
atac    : --nomodel --shift -100 --extsize 200      (open chromatin)
chip    : model-based or --nomodel; narrow or broad  (histone marks / TFs)
cuttag  : --nomodel --shift -75  --extsize 150      (CUT&Tag sub-nucleosomal)
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .utils import makedirs, require_tools, run_cmd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# MACS3 parameter presets
# ---------------------------------------------------------------------------

_PRESETS: Dict[str, Dict[str, Any]] = {
    "atac": {
        "nomodel": True,
        "shift": -100,
        "extsize": 200,
        "broad": False,
    },
    "chip": {
        "nomodel": False,
        "shift": 0,
        "extsize": 200,
        "broad": False,
    },
    "cuttag": {
        "nomodel": True,
        "shift": -75,
        "extsize": 150,
        "broad": False,
    },
}


def get_macs3_params(assay: str, cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Merge default preset with per-assay config overrides."""
    preset = _PRESETS.get(assay, _PRESETS["atac"]).copy()
    overrides = cfg.get("macs3", {}).get(assay, {})
    preset.update({k: v for k, v in overrides.items() if v is not None})
    return preset


# ---------------------------------------------------------------------------
# MACS3 callpeak
# ---------------------------------------------------------------------------

def call_peaks(
    treatment_bam: Union[str, Path],
    out_dir: Union[str, Path],
    sample_name: str,
    assay: str = "atac",
    control_bam: Optional[Union[str, Path]] = None,
    genome_size: str = "hs",
    qvalue: float = 0.05,
    params: Optional[Dict[str, Any]] = None,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Call peaks using MACS3.

    Parameters
    ----------
    treatment_bam : Path to treatment BAM.
    out_dir       : Output directory for MACS3 results.
    sample_name   : Prefix for output files.
    assay         : 'atac' | 'chip' | 'cuttag'
    control_bam   : Optional input/IgG control BAM (ChIP-seq).
    genome_size   : MACS3 -g flag (hs, mm, ce, dm, or integer).
    qvalue        : q-value cutoff.
    params        : Dict of extra MACS3 options (overrides preset).

    Returns
    -------
    Path to the narrowPeak (or broadPeak) file.
    """
    require_tools("macs3")
    makedirs(out_dir)

    p = params or {}
    nomodel_flag = "--nomodel" if p.get("nomodel", True) else ""
    shift_flag = f"--shift {p.get('shift', -100)}" if p.get("nomodel", True) else ""
    extsize_flag = f"--extsize {p.get('extsize', 200)}" if p.get("nomodel", True) else ""
    broad_flag = "--broad" if p.get("broad", False) else ""
    broad_cutoff = f"--broad-cutoff {p.get('broad_cutoff', 0.05)}" if p.get("broad", False) else ""
    ctrl_flag = f"-c {control_bam}" if control_bam else ""
    extra = p.get("extra", "")

    peak_format = "BAMPE"   # Paired-end; use BAM for SE

    cmd = (
        f"macs3 callpeak "
        f"-t {treatment_bam} {ctrl_flag} "
        f"-f {peak_format} -g {genome_size} -n {sample_name} "
        f"--outdir {out_dir} -q {qvalue} "
        f"{nomodel_flag} {shift_flag} {extsize_flag} "
        f"{broad_flag} {broad_cutoff} {extra}"
    )
    log.info("MACS3 (%s): %s", assay, treatment_bam)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)

    ext = "broadPeak" if p.get("broad", False) else "narrowPeak"
    return Path(out_dir) / f"{sample_name}_peaks.{ext}"


# ---------------------------------------------------------------------------
# Peak filtering helpers
# ---------------------------------------------------------------------------

def filter_peaks(
    peak_bed: Union[str, Path],
    blacklist: Union[str, Path],
    out_bed: Union[str, Path],
    dry_run: bool = False,
) -> Path:
    """Remove peaks overlapping blacklist regions using bedtools."""
    require_tools("bedtools")
    makedirs(Path(out_bed).parent)
    cmd = f"bedtools intersect -v -a {peak_bed} -b {blacklist} > {out_bed}"
    log.info("Filter peaks against blacklist: %s", blacklist)
    run_cmd(cmd, dry_run=dry_run)
    return Path(out_bed)


def read_peaks(peak_bed: Union[str, Path]) -> pd.DataFrame:
    """
    Read a narrowPeak / broadPeak / BED file into a DataFrame.

    Columns: chr, start, end, name, score, strand, signal, pvalue, qvalue, [summit]
    """
    cols = ["chr", "start", "end", "name", "score", "strand",
            "signal", "pvalue", "qvalue", "summit"]
    df = pd.read_csv(peak_bed, sep="\t", header=None, comment="#")
    df.columns = cols[: len(df.columns)]
    return df


def peak_stats(peak_bed: Union[str, Path]) -> Dict[str, Any]:
    """Return summary statistics for a peak file."""
    df = read_peaks(peak_bed)
    widths = df["end"] - df["start"]
    return {
        "n_peaks": len(df),
        "median_width_bp": int(widths.median()),
        "mean_width_bp": float(widths.mean()),
        "total_bp": int(widths.sum()),
    }
