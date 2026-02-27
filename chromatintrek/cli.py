"""
chromatintrek.cli
=================
Command-line interface for ChromatinTrek.

Usage
-----
chromatintrek run    --config config.yaml [--slurm] [--dry-run]
chromatintrek setup-genome  --species human|mouse [--genome-dir PATH]
chromatintrek check-env
"""

import argparse
import logging
import sys
from pathlib import Path

log = logging.getLogger("chromatintrek")


def _run(args: argparse.Namespace) -> None:
    """Launch the Snakemake pipeline."""
    import subprocess
    snakefile = Path(__file__).parent.parent / "Snakefile"
    cores = args.cores if hasattr(args, "cores") else 1
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--configfile", args.config,
        "--cores", str(cores),
        "--use-conda",
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    if args.dry_run:
        cmd.append("--dry-run")
    if args.slurm:
        profile_dir = Path(__file__).parent.parent / "profiles" / "slurm"
        cmd += ["--executor", "slurm", "--profile", str(profile_dir)]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def _setup_genome(args: argparse.Namespace) -> None:
    """Download and index reference genome files."""
    from .genome import resolve
    gcfg = resolve(args.species, genome_dir=args.genome_dir or "")
    print(f"Genome: {gcfg.genome}")
    print(f"Genome dir: {gcfg.genome_dir}")
    print("\nDownloading genome files (this may take a while)…")
    _download_genome(gcfg)


def _download_genome(gcfg) -> None:
    """Download genome FASTA, GTF, blacklist; build Bowtie2 index."""
    import subprocess, os
    genome = gcfg.genome
    gdir = Path(gcfg.genome_dir)
    gdir.mkdir(parents=True, exist_ok=True)

    urls = {
        "hg38": {
            "fasta": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            "blacklist": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
            "chrom_sizes": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
        },
        "mm10": {
            "fasta": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
            "blacklist": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
            "chrom_sizes": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes",
        },
    }
    if genome not in urls:
        print(f"No auto-download URL for genome {genome}. Please download manually.")
        return

    for name, url in urls[genome].items():
        dest = gdir / Path(url).name
        if dest.exists():
            print(f"  [skip] {dest.name} already exists")
            continue
        print(f"  Downloading {name} …")
        subprocess.run(["wget", "-q", "-O", str(dest), url], check=True)

    # Build Bowtie2 index
    fasta = gdir / f"{genome}.fa.gz"
    idx_prefix = gdir / "bowtie2" / genome
    idx_prefix.parent.mkdir(parents=True, exist_ok=True)
    if not (idx_prefix.parent / f"{genome}.1.bt2").exists():
        print("  Building Bowtie2 index (this takes ~1 hour) …")
        subprocess.run(
            ["bowtie2-build", "--threads", "8", str(fasta), str(idx_prefix)],
            check=True,
        )
    print(f"\n✓ Genome setup complete: {gdir}")


def _check_env(_args: argparse.Namespace) -> None:
    """Check that all required tools are on PATH."""
    from .utils import check_tools
    required = [
        "fastqc", "trim_galore", "bowtie2", "samtools",
        "macs3", "findMotifsGenome.pl", "annotatePeaks.pl",
        "bamCoverage", "computeMatrix", "plotHeatmap", "plotProfile",
        "multiqc", "bedtools",
    ]
    status = check_tools(*required)
    ok = all(status.values())
    print("ChromatinTrek environment check\n" + "─" * 40)
    for tool, found in status.items():
        mark = "✓" if found else "✗"
        print(f"  {mark}  {tool}")
    print("─" * 40)
    if ok:
        print("All tools found. ✓")
    else:
        missing = [t for t, f in status.items() if not f]
        print(f"\n⚠  Missing tools: {missing}")
        print("Run:  conda activate chromatintrek")
        sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="chromatintrek",
        description="ChromatinTrek – End-to-End Chromatin NGS Pipeline",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # run
    p_run = sub.add_parser("run", help="Execute the pipeline")
    p_run.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    p_run.add_argument("--cores", type=int, default=1, help="Local CPU cores (ignored with --slurm)")
    p_run.add_argument("--slurm", action="store_true", help="Submit jobs via SLURM")
    p_run.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    p_run.set_defaults(func=_run)

    # setup-genome
    p_sg = sub.add_parser("setup-genome", help="Download and index reference genome")
    p_sg.add_argument("--species", required=True, choices=["human", "mouse"])
    p_sg.add_argument("--genome-dir", default="", help="Override default genome directory")
    p_sg.set_defaults(func=_setup_genome)

    # check-env
    p_ce = sub.add_parser("check-env", help="Verify all required tools are installed")
    p_ce.set_defaults(func=_check_env)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args.func(args)


if __name__ == "__main__":
    main()
