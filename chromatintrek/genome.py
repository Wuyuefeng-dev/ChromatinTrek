"""
chromatintrek.genome
====================
Auto-resolve species → genome configuration.

Calling resolve(species) returns a Namespace with all genome-level
parameters derived from just "human" or "mouse":

  genome          hg38 / mm10
  bowtie2_index   $CONDA_PREFIX/share/chromatintrek/genomes/<genome>/bowtie2/<genome>
  chromosome_sizes, blacklist, gtf, tss_bed  (same prefix)
  effective_genome_size
  homer_genome    hg38 / mm10
  go_organism     "Human" / "Mouse"
  macs3_gsize     "hs" / "mm"

Users can override any field via the advanced config keys.
"""

import os
from pathlib import Path
from typing import Any, Dict


# ---------------------------------------------------------------------------
# Species → genome parameter table
# ---------------------------------------------------------------------------

_SPECIES_MAP: Dict[str, Dict[str, Any]] = {
    "human": {
        "genome":                  "hg38",
        "homer_genome":            "hg38",
        "go_organism":             "Human",
        "macs3_gsize":             "hs",
        "effective_genome_size":   2_913_022_398,
    },
    "mouse": {
        "genome":                  "mm10",
        "homer_genome":            "mm10",
        "go_organism":             "Mouse",
        "macs3_gsize":             "mm",
        "effective_genome_size":   2_652_783_500,
    },
}

_GENOME_FILES = {
    "bowtie2_index":    "bowtie2/{genome}",
    "chromosome_sizes": "{genome}.chrom.sizes",
    "blacklist":        "{genome}-blacklist.v2.bed",
    "gtf":              "annotation/{genome}.gtf.gz",
    "tss_bed":          "annotation/{genome}_tss.bed",
}


class GenomeConfig:
    """Container for resolved genome parameters."""

    def __init__(self, **kwargs: Any):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def as_dict(self) -> Dict[str, Any]:
        return self.__dict__.copy()

    def __repr__(self) -> str:
        body = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"GenomeConfig({body})"


def resolve(species: str, genome_dir: str = "", user_overrides: Dict[str, Any] = None) -> GenomeConfig:
    """
    Return a GenomeConfig resolved from *species* ("human" | "mouse").

    Parameters
    ----------
    species        : "human" or "mouse"
    genome_dir     : root directory for genome files; defaults to
                     $CONDA_PREFIX/share/chromatintrek/genomes/<genome>/
    user_overrides : dict of keys that override auto-resolved values
                     (from advanced section of config.yaml)
    """
    species = species.lower().strip()
    if species not in _SPECIES_MAP:
        raise ValueError(
            f"Unknown species '{species}'. Supported: {list(_SPECIES_MAP)}"
        )
    params = _SPECIES_MAP[species].copy()
    genome = params["genome"]

    # Resolve genome directory
    if not genome_dir:
        conda_prefix = os.environ.get("CONDA_PREFIX", "")
        if conda_prefix:
            genome_dir = str(Path(conda_prefix) / "share" / "chromatintrek" / "genomes" / genome)
        else:
            genome_dir = str(Path.home() / ".chromatintrek" / "genomes" / genome)

    params["genome_dir"] = genome_dir

    # Resolve file paths
    for key, tmpl in _GENOME_FILES.items():
        rel = tmpl.format(genome=genome)
        params[key] = str(Path(genome_dir) / rel)

    # Override with any user-provided advanced values
    if user_overrides:
        for k, v in user_overrides.items():
            if v:   # only override if non-empty
                params[k] = v

    return GenomeConfig(**params)


def from_config(cfg: Dict[str, Any]) -> GenomeConfig:
    """Convenience: build GenomeConfig from a loaded pipeline config dict."""
    species = cfg.get("species", "human")
    genome_dir = cfg.get("genome_dir", "")
    overrides = {
        k: cfg.get(k)
        for k in ("bowtie2_index", "chromosome_sizes", "blacklist", "gtf", "tss_bed")
    }
    return resolve(species, genome_dir=genome_dir, user_overrides=overrides)
