"""
chromatintrek.qc
================
Quality control: FastQC, Trim Galore, MultiQC wrappers.

All functions operate via subprocess calls to conda-installed tools.
Provide dry_run=True to print commands without executing.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .utils import makedirs, require_tools, run_cmd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# FastQC
# ---------------------------------------------------------------------------

def run_fastqc(
    fastq_files: List[Union[str, Path]],
    out_dir: Union[str, Path],
    threads: int = 4,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> None:
    """
    Run FastQC on one or more FASTQ files.

    Parameters
    ----------
    fastq_files : list of FASTQ paths (R1, and optionally R2)
    out_dir     : directory for FastQC output HTML/zip
    threads     : number of parallel threads
    logfile     : path to capture stdout + stderr
    dry_run     : print command without executing
    """
    require_tools("fastqc")
    makedirs(out_dir)
    fq_str = " ".join(str(f) for f in fastq_files)
    cmd = f"fastqc --threads {threads} --outdir {out_dir} {fq_str}"
    log.info("FastQC: %s", fq_str)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)


# ---------------------------------------------------------------------------
# Trim Galore
# ---------------------------------------------------------------------------

def run_trim_galore(
    fastq_r1: Union[str, Path],
    out_dir: Union[str, Path],
    fastq_r2: Optional[Union[str, Path]] = None,
    quality: int = 20,
    min_length: int = 25,
    cores: int = 4,
    extra: str = "",
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Dict[str, Path]:
    """
    Run Trim Galore to quality- and adapter-trim reads.

    Returns paths to trimmed FASTQ(s).
    """
    require_tools("trim_galore")
    makedirs(out_dir)
    paired = fastq_r2 is not None
    paired_flag = "--paired" if paired else ""
    r2_str = str(fastq_r2) if paired else ""
    cmd = (
        f"trim_galore {paired_flag} --quality {quality} --length {min_length} "
        f"--cores {cores} --fastqc --output_dir {out_dir} {extra} "
        f"{fastq_r1} {r2_str}"
    )
    log.info("Trim Galore: %s %s", fastq_r1, r2_str)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)

    stem = Path(fastq_r1).stem.replace(".fastq", "").replace(".fq", "")
    result: Dict[str, Path] = {"r1": Path(out_dir) / f"{stem}_trimmed.fq.gz"}
    if paired:
        stem2 = Path(fastq_r2).stem.replace(".fastq", "").replace(".fq", "")
        result["r2"] = Path(out_dir) / f"{stem2}_trimmed.fq.gz"
    return result


# ---------------------------------------------------------------------------
# MultiQC
# ---------------------------------------------------------------------------

def run_multiqc(
    search_dir: Union[str, Path],
    out_dir: Union[str, Path],
    title: str = "ChromatinTrek QC Report",
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> None:
    """
    Aggregate FastQC/Trim Galore logs with MultiQC into a single HTML report.

    Parameters
    ----------
    search_dir : directory to scan for tool outputs
    out_dir    : directory for MultiQC report
    title      : report title
    """
    require_tools("multiqc")
    makedirs(out_dir)
    cmd = f'multiqc --title "{title}" --outdir {out_dir} --force {search_dir}'
    log.info("MultiQC: scanning %s", search_dir)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)


# ---------------------------------------------------------------------------
# Convenience: full QC step for one sample
# ---------------------------------------------------------------------------

def qc_sample(
    sample: str,
    cfg: Dict[str, Any],
    fastq_r1: Union[str, Path],
    fastq_r2: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Dict[str, Path]:
    """
    Run FastQC -> Trim Galore -> FastQC (post-trim) for a single sample.

    Returns dict with trimmed FASTQ paths.
    """
    qc_dir = Path(cfg["qc_dir"]) / sample
    trim_dir = Path(cfg["trim_dir"]) / sample
    log_dir = Path(cfg["log_dir"]) / sample

    # Pre-trim FastQC
    fqs = [fastq_r1] + ([fastq_r2] if fastq_r2 else [])
    run_fastqc(
        fqs, qc_dir / "pre_trim",
        threads=cfg["fastqc"]["threads"],
        logfile=log_dir / "fastqc_pre.log",
        dry_run=dry_run,
    )

    # Trimming
    trim_params = cfg.get("trim_galore", {})
    trimmed = run_trim_galore(
        fastq_r1, trim_dir,
        fastq_r2=fastq_r2,
        quality=trim_params.get("quality", 20),
        min_length=trim_params.get("length", 25),
        cores=trim_params.get("cores", 4),
        extra=trim_params.get("extra", ""),
        logfile=log_dir / "trim_galore.log",
        dry_run=dry_run,
    )

    # Post-trim FastQC
    trimmed_fqs = list(trimmed.values())
    run_fastqc(
        trimmed_fqs, qc_dir / "post_trim",
        threads=cfg["fastqc"]["threads"],
        logfile=log_dir / "fastqc_post.log",
        dry_run=dry_run,
    )

    return trimmed
