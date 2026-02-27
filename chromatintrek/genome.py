"""
chromatintrek.genome
====================
Auto-resolve species → genome configuration.

Priority order (highest wins):
  1. Paths explicitly set in config.yaml  (genome_fasta, bowtie2_index, …)
  2. Paths under genome_dir (auto-managed by `chromatintrek setup-genome`)
  3. $CONDA_PREFIX/share/chromatintrek/genomes/<genome>/

When a bowtie2_index is already provided, setup-genome skips building it.
When a genome_fasta is already provided, setup-genome skips downloading it.
"""

import os
from pathlib import Path
from typing import Any, Dict, Optional


# ---------------------------------------------------------------------------
# Species → genome parameter table
# ---------------------------------------------------------------------------

_SPECIES_MAP: Dict[str, Dict[str, Any]] = {
    "human": {
        "genome":                "hg38",
        "homer_genome":          "hg38",
        "go_organism":           "Human",
        "macs3_gsize":           "hs",
        "effective_genome_size": 2_913_022_398,
    },
    "mouse": {
        "genome":                "mm10",
        "homer_genome":          "mm10",
        "go_organism":           "Mouse",
        "macs3_gsize":           "mm",
        "effective_genome_size": 2_652_783_500,
    },
}

# Relative paths under genome_dir
_RELATIVE_FILES = {
    "bowtie2_index":    "bowtie2/{genome}",        # prefix (no extension)
    "chromosome_sizes": "{genome}.chrom.sizes",
    "blacklist":        "{genome}-blacklist.v2.bed",
    "gtf":              "annotation/{genome}.gtf.gz",
    "tss_bed":          "annotation/{genome}_tss.bed",
    "genome_fasta":     "{genome}.fa.gz",
}

# Download URLs
_DOWNLOAD_URLS: Dict[str, Dict[str, str]] = {
    "hg38": {
        "genome_fasta":     "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "blacklist":        "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
        "chromosome_sizes": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    },
    "mm10": {
        "genome_fasta":     "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
        "blacklist":        "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
        "chromosome_sizes": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes",
    },
}


# ---------------------------------------------------------------------------
# GenomeConfig container
# ---------------------------------------------------------------------------

class GenomeConfig:
    """Resolved genome parameters used throughout the pipeline."""

    def __init__(self, **kwargs: Any):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def as_dict(self) -> Dict[str, Any]:
        return self.__dict__.copy()

    def __repr__(self) -> str:
        body = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"GenomeConfig({body})"


# ---------------------------------------------------------------------------
# Resolution logic
# ---------------------------------------------------------------------------

def resolve(
    species: str,
    genome_dir: str = "",
    user_overrides: Optional[Dict[str, Any]] = None,
) -> GenomeConfig:
    """
    Return a GenomeConfig for *species* ("human" | "mouse").

    User-supplied paths in *user_overrides* take the highest priority and
    suppress auto-download/build for those specific files.

    Parameters
    ----------
    species        : "human" or "mouse"
    genome_dir     : Root directory for auto-managed genome files.
                     Defaults to $CONDA_PREFIX/share/chromatintrek/genomes/<genome>/
    user_overrides : Dict of config keys → user-supplied paths.
                     Non-empty values override auto-resolved paths.
                     E.g. {"bowtie2_index": "/scratch/hg38/bowtie2/hg38",
                           "genome_fasta":  "/scratch/hg38/hg38.fa.gz"}
    """
    species = species.lower().strip()
    if species not in _SPECIES_MAP:
        raise ValueError(f"Unknown species '{species}'. Supported: {list(_SPECIES_MAP)}")

    params = _SPECIES_MAP[species].copy()
    genome = params["genome"]

    # Resolve genome_dir
    if not genome_dir:
        conda_prefix = os.environ.get("CONDA_PREFIX", "")
        if conda_prefix:
            genome_dir = str(Path(conda_prefix) / "share" / "chromatintrek" / "genomes" / genome)
        else:
            genome_dir = str(Path.home() / ".chromatintrek" / "genomes" / genome)
    params["genome_dir"] = genome_dir

    # Auto-resolve each file path from genome_dir
    for key, tmpl in _RELATIVE_FILES.items():
        params[key] = str(Path(genome_dir) / tmpl.format(genome=genome))

    # Apply user overrides (non-empty values win)
    if user_overrides:
        for k, v in user_overrides.items():
            if v and str(v).strip():
                params[k] = str(v).strip()

    return GenomeConfig(**params)


def from_config(cfg: Dict[str, Any]) -> GenomeConfig:
    """
    Build a GenomeConfig from a loaded pipeline config.yaml dict.

    User-supplied paths (genome_fasta, bowtie2_index, blacklist, etc.)
    are passed as overrides so they take priority over auto-resolved paths.
    """
    override_keys = [
        "genome_fasta", "bowtie2_index",
        "chromosome_sizes", "blacklist", "gtf", "tss_bed",
    ]
    overrides = {k: cfg.get(k, "") for k in override_keys}
    return resolve(
        species=cfg.get("species", "human"),
        genome_dir=cfg.get("genome_dir", ""),
        user_overrides=overrides,
    )


# ---------------------------------------------------------------------------
# Setup-genome helper (called by CLI)
# ---------------------------------------------------------------------------

def setup_genome(gcfg: GenomeConfig, threads: int = 8) -> None:
    """
    Download genome files and build Bowtie2 index if not already present.

    Files already specified by the user (non-default paths) are skipped.
    """
    import subprocess
    genome = gcfg.genome
    gdir = Path(gcfg.genome_dir)
    gdir.mkdir(parents=True, exist_ok=True)

    urls = _DOWNLOAD_URLS.get(genome, {})

    # ── Download FASTA ────────────────────────────────────────────────────
    fasta_path = Path(gcfg.genome_fasta)
    if fasta_path.exists():
        print(f"  [skip] genome FASTA already exists: {fasta_path}")
    elif "genome_fasta" in urls:
        print(f"  Downloading genome FASTA ({genome}) …")
        subprocess.run(["wget", "-q", "-O", str(fasta_path), urls["genome_fasta"]], check=True)
    else:
        print(f"  No download URL for {genome} FASTA. Provide path in config.yaml.")

    # ── Download blacklist ────────────────────────────────────────────────
    bl_path = Path(gcfg.blacklist)
    if bl_path.exists():
        print(f"  [skip] blacklist already exists: {bl_path}")
    elif "blacklist" in urls:
        print(f"  Downloading blacklist …")
        raw = str(bl_path) + ".gz"
        subprocess.run(["wget", "-q", "-O", raw, urls["blacklist"]], check=True)
        subprocess.run(["gunzip", "-f", raw], check=True)

    # ── Download chrom sizes ───────────────────────────────────────────────
    cs_path = Path(gcfg.chromosome_sizes)
    if cs_path.exists():
        print(f"  [skip] chrom.sizes already exists: {cs_path}")
    elif "chromosome_sizes" in urls:
        print(f"  Downloading chrom.sizes …")
        subprocess.run(["wget", "-q", "-O", str(cs_path), urls["chromosome_sizes"]], check=True)

    # ── Build Bowtie2 index ────────────────────────────────────────────────
    idx_prefix = Path(gcfg.bowtie2_index)
    idx_check = Path(str(idx_prefix) + ".1.bt2")
    if idx_check.exists():
        print(f"  [skip] Bowtie2 index already exists: {idx_prefix}")
    elif fasta_path.exists():
        print(f"  Building Bowtie2 index (this takes ~60 min) …")
        idx_prefix.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["bowtie2-build", "--threads", str(threads), str(fasta_path), str(idx_prefix)],
            check=True,
        )
    else:
        print(
            f"  ⚠  Cannot build Bowtie2 index: FASTA not found at {fasta_path}.\n"
            "     Provide 'genome_fasta' or 'bowtie2_index' in config.yaml."
        )

    print(f"\n✓ Genome setup complete: {gdir}")
