"""
chromatintrek.alignment
=======================
Read alignment (Bowtie2) and BAM processing (samtools / Picard).

Supports assay-specific alignment flags:
  atac    --no-mixed --no-discordant -X 2000
  chip    --no-mixed --no-discordant
  cuttag  --no-mixed --no-discordant -X 700

Post-alignment steps: sort, index, mark duplicates, MAPQ filter,
remove reads overlapping blacklist regions, generate flagstat and
fragment-size distribution.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .utils import makedirs, require_tools, run_cmd

log = logging.getLogger(__name__)

# Assay-specific Bowtie2 flags
_BOWTIE2_FLAGS = {
    "atac":    "--no-mixed --no-discordant -X 2000",
    "chip":    "--no-mixed --no-discordant",
    "cuttag":  "--no-mixed --no-discordant -X 700",
}


# ---------------------------------------------------------------------------
# Bowtie2 alignment
# ---------------------------------------------------------------------------

def run_bowtie2(
    index: str,
    fastq_r1: Union[str, Path],
    out_bam: Union[str, Path],
    fastq_r2: Optional[Union[str, Path]] = None,
    assay: str = "atac",
    threads: int = 16,
    extra_flags: str = "",
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Align reads to the genome with Bowtie2 and pipe directly to a
    coordinate-sorted, indexed BAM via samtools.

    Returns path to the output BAM.
    """
    require_tools("bowtie2", "samtools")
    makedirs(Path(out_bam).parent)

    flags = _BOWTIE2_FLAGS.get(assay, "") + " " + extra_flags
    if fastq_r2:
        input_str = f"-1 {fastq_r1} -2 {fastq_r2}"
    else:
        input_str = f"-U {fastq_r1}"

    cmd = (
        f"bowtie2 -p {threads} -x {index} {flags} {input_str} "
        f"| samtools sort -@ {threads} -o {out_bam} - "
        f"&& samtools index {out_bam}"
    )
    log.info("Bowtie2 (%s): %s -> %s", assay, fastq_r1, out_bam)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_bam)


# ---------------------------------------------------------------------------
# Duplicate marking (samtools markdup)
# ---------------------------------------------------------------------------

def mark_duplicates(
    in_bam: Union[str, Path],
    out_bam: Union[str, Path],
    threads: int = 8,
    remove_dups: bool = True,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Mark (and optionally remove) PCR duplicates using samtools markdup.

    For ATAC-seq and CUT&Tag it is standard practice to remove duplicates.
    """
    require_tools("samtools")
    makedirs(Path(out_bam).parent)
    remove_flag = "-r" if remove_dups else ""
    tmp_sorted = str(out_bam) + ".namesorted.bam"

    cmd = (
        f"samtools sort -n -@ {threads} -o {tmp_sorted} {in_bam} && "
        f"samtools fixmate -m {tmp_sorted} - | "
        f"samtools sort -@ {threads} - | "
        f"samtools markdup {remove_flag} -@ {threads} - {out_bam} && "
        f"samtools index {out_bam} && "
        f"rm -f {tmp_sorted}"
    )
    log.info("Mark duplicates: %s -> %s", in_bam, out_bam)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_bam)


# ---------------------------------------------------------------------------
# MAPQ filtering & blacklist removal
# ---------------------------------------------------------------------------

def filter_bam(
    in_bam: Union[str, Path],
    out_bam: Union[str, Path],
    blacklist: Optional[Union[str, Path]] = None,
    mapq: int = 30,
    threads: int = 8,
    paired_end: bool = True,
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> Path:
    """
    Apply MAPQ filter and optionally remove reads overlapping a blacklist BED.

    For paired-end data, also filters unmapped and improperly paired reads.
    """
    require_tools("samtools")
    if blacklist:
        require_tools("bedtools")
    makedirs(Path(out_bam).parent)

    pe_flags = "-f 2" if paired_end else ""   # keep properly paired
    filter_cmd = f"samtools view -@ {threads} -q {mapq} {pe_flags} -b {in_bam}"

    if blacklist:
        cmd = (
            f"{filter_cmd} "
            f"| bedtools intersect -v -abam stdin -b {blacklist} "
            f"| samtools sort -@ {threads} -o {out_bam} - "
            f"&& samtools index {out_bam}"
        )
    else:
        cmd = (
            f"{filter_cmd} "
            f"| samtools sort -@ {threads} -o {out_bam} - "
            f"&& samtools index {out_bam}"
        )
    log.info("Filter BAM (mapq>=%d): %s -> %s", mapq, in_bam, out_bam)
    run_cmd(cmd, logfile=logfile, dry_run=dry_run)
    return Path(out_bam)


# ---------------------------------------------------------------------------
# Flagstat QC
# ---------------------------------------------------------------------------

def flagstat(
    bam: Union[str, Path],
    out_txt: Union[str, Path],
    threads: int = 4,
    dry_run: bool = False,
) -> None:
    """Write samtools flagstat summary to *out_txt*."""
    require_tools("samtools")
    makedirs(Path(out_txt).parent)
    cmd = f"samtools flagstat -@ {threads} {bam} > {out_txt}"
    log.info("flagstat: %s", bam)
    run_cmd(cmd, dry_run=dry_run)


# ---------------------------------------------------------------------------
# Full alignment pipeline for one sample
# ---------------------------------------------------------------------------

def align_sample(
    sample: str,
    cfg: Dict[str, Any],
    fastq_r1: Union[str, Path],
    fastq_r2: Optional[Union[str, Path]] = None,
    assay: str = "atac",
    dry_run: bool = False,
) -> Path:
    """
    Run the complete alignment pipeline for one sample:
      Bowtie2 → mark duplicates → MAPQ filter → blacklist removal → flagstat

    Returns path to the final filtered BAM.
    """
    align_dir = Path(cfg["align_dir"]) / sample
    log_dir = Path(cfg["log_dir"]) / sample
    makedirs(align_dir, log_dir)

    bt2_cfg = cfg.get("bowtie2", {})
    sam_cfg = cfg.get("samtools", {})
    threads_bt2 = bt2_cfg.get("threads", 16)
    threads_sam = sam_cfg.get("threads", 8)

    raw_bam = align_dir / f"{sample}.raw.bam"
    dedup_bam = align_dir / f"{sample}.dedup.bam"
    final_bam = align_dir / f"{sample}.filtered.bam"

    # Step 1: Align
    run_bowtie2(
        cfg["bowtie2_index"], fastq_r1, raw_bam,
        fastq_r2=fastq_r2, assay=assay, threads=threads_bt2,
        logfile=log_dir / "bowtie2.log", dry_run=dry_run,
    )

    # Step 2: Deduplicate
    mark_duplicates(
        raw_bam, dedup_bam, threads=threads_sam,
        remove_dups=True,
        logfile=log_dir / "markdup.log", dry_run=dry_run,
    )

    # Step 3: MAPQ filter + blacklist removal
    filter_bam(
        dedup_bam, final_bam,
        blacklist=cfg.get("blacklist"),
        mapq=sam_cfg.get("mapq", 30),
        threads=threads_sam,
        paired_end=fastq_r2 is not None,
        logfile=log_dir / "filter.log", dry_run=dry_run,
    )

    # Step 4: Flagstat
    flagstat(final_bam, align_dir / f"{sample}.flagstat.txt",
             threads=threads_sam, dry_run=dry_run)

    return final_bam
