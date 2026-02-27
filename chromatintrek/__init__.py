"""
chromatintrek
=============
End-to-end pipeline for ATAC-seq, ChIP-seq, and CUT&Tag chromatin analysis.

Modules
-------
utils              Shared utilities: config, logging, subprocess, sample sheet
qc                 Quality control: FastQC, Trim Galore, MultiQC
alignment          Read alignment (Bowtie2) and BAM processing (samtools/Picard)
peak_calling       Peak calling via MACS3 with per-assay parameter presets
peak_visualization Signal visualisation via deepTools (bigwig, heatmap, profile)
motif_tf           Motif finding and TF enrichment via HOMER
peak_gene          Peak-to-gene annotation and distance statistics
go_enrichment      Gene Ontology and pathway enrichment via gseapy
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("chromatintrek")
except PackageNotFoundError:
    __version__ = "0.1.0"

__author__ = "ChromatinTrek"
__license__ = "MIT"

from . import (
    alignment,
    go_enrichment,
    motif_tf,
    peak_calling,
    peak_gene,
    peak_visualization,
    qc,
    utils,
)

__all__ = [
    "utils",
    "qc",
    "alignment",
    "peak_calling",
    "peak_visualization",
    "motif_tf",
    "peak_gene",
    "go_enrichment",
]
