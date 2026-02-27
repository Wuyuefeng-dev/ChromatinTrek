"""
chromatintrek.peak_visualization
=================================
Signal visualization using deepTools: bigwig generation, signal matrix,
heatmap, average profile, and fingerprint (library complexity / enrichment).
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .utils import makedirs, require_tools, run_cmd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# BAM → BigWig
# ---------------------------------------------------------------------------

def bam_to_bigwig(
    bam: Union[str, Path],
    out_bw: Union[str, Path],
    effective_genome_size: int = 2913022398,  # hg38
    normalize: str = "RPGC",
    bin_size: int = 10,
    threads: int = 16,
    blacklist: Optional[Union[str, Path]] = None,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Convert a BAM file to a normalised bigWig using deepTools bamCoverage.

    Parameters
    ----------
    normalize : RPKM | CPM | BPM | RPGC (default RPGC for ATAC/CUT&Tag)
    """
    require_tools("bamCoverage")
    makedirs(Path(out_bw).parent)
    bl_flag = f"--blackListFileName {blacklist}" if blacklist else ""
    cmd = (
        f"bamCoverage -b {bam} -o {out_bw} "
        f"--binSize {bin_size} --normalizeUsing {normalize} "
        f"--effectiveGenomeSize {effective_genome_size} "
        f"--extendReads -p {threads} {bl_flag}"
    )
    log.info("bamCoverage: %s -> %s", bam, out_bw)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_bw)


# ---------------------------------------------------------------------------
# computeMatrix (reference-point or scale-regions)
# ---------------------------------------------------------------------------

def compute_matrix(
    bigwig_files: List[Union[str, Path]],
    regions_bed: Union[str, Path],
    out_matrix: Union[str, Path],
    mode: str = "reference-point",
    upstream: int = 3000,
    downstream: int = 3000,
    bin_size: int = 50,
    reference_point: str = "center",
    threads: int = 16,
    sample_labels: Optional[List[str]] = None,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Run deepTools computeMatrix to generate a signal matrix around genomic regions.

    Parameters
    ----------
    mode            : 'reference-point' (default) or 'scale-regions'
    reference_point : 'center' | 'TES' | 'TSS' (for reference-point mode)
    """
    require_tools("computeMatrix")
    makedirs(Path(out_matrix).parent)
    bw_str = " ".join(str(b) for b in bigwig_files)
    label_flag = (
        f"--samplesLabel {' '.join(sample_labels)}" if sample_labels else ""
    )

    if mode == "reference-point":
        mode_flags = (
            f"reference-point --referencePoint {reference_point} "
            f"-b {upstream} -a {downstream}"
        )
    else:
        mode_flags = f"scale-regions -b {upstream} -a {downstream}"

    cmd = (
        f"computeMatrix {mode_flags} "
        f"-S {bw_str} -R {regions_bed} "
        f"--binSize {bin_size} -p {threads} "
        f"-o {out_matrix} {label_flag}"
    )
    log.info("computeMatrix [%s]: %s", mode, regions_bed)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_matrix)


# ---------------------------------------------------------------------------
# plotHeatmap
# ---------------------------------------------------------------------------

def plot_heatmap(
    matrix: Union[str, Path],
    out_png: Union[str, Path],
    color_map: str = "Blues",
    sort_regions: str = "descend",
    title: str = "",
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """Render a sorted signal heatmap from a deepTools matrix gz file."""
    require_tools("plotHeatmap")
    makedirs(Path(out_png).parent)
    title_flag = f'--plotTitle "{title}"' if title else ""
    cmd = (
        f"plotHeatmap -m {matrix} -out {out_png} "
        f"--colorMap {color_map} --sortRegions {sort_regions} "
        f"--dpi 300 {title_flag}"
    )
    log.info("plotHeatmap -> %s", out_png)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_png)


# ---------------------------------------------------------------------------
# plotProfile (average signal)
# ---------------------------------------------------------------------------

def plot_profile(
    matrix: Union[str, Path],
    out_png: Union[str, Path],
    title: str = "",
    per_group: bool = False,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """Render an average signal profile (line plot) from a deepTools matrix."""
    require_tools("plotProfile")
    makedirs(Path(out_png).parent)
    per_group_flag = "--perGroup" if per_group else ""
    title_flag = f'--plotTitle "{title}"' if title else ""
    cmd = (
        f"plotProfile -m {matrix} -out {out_png} "
        f"--dpi 300 {per_group_flag} {title_flag}"
    )
    log.info("plotProfile -> %s", out_png)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_png)


# ---------------------------------------------------------------------------
# plotFingerprint (library enrichment / complexity)
# ---------------------------------------------------------------------------

def plot_fingerprint(
    bam_files: List[Union[str, Path]],
    out_png: Union[str, Path],
    labels: Optional[List[str]] = None,
    threads: int = 8,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Run deepTools plotFingerprint to visualise cumulative enrichment.

    A good ChIP/ATAC sample shows a strong curve bending away from the
    diagonal (= the IgG / input control).
    """
    require_tools("plotFingerprint")
    makedirs(Path(out_png).parent)
    bam_str = " ".join(str(b) for b in bam_files)
    label_flag = f"--labels {' '.join(labels)}" if labels else ""
    cmd = (
        f"plotFingerprint -b {bam_str} -plot {out_png} "
        f"-p {threads} {label_flag} --dpi 300"
    )
    log.info("plotFingerprint -> %s", out_png)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_png)
