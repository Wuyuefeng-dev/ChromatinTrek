# =============================================================================
# ChromatinTrek – Main Snakefile
# =============================================================================
# Run locally:  snakemake --cores 8 --use-conda
# Run on SLURM: snakemake --executor slurm --profile profiles/slurm
# =============================================================================

from pathlib import Path
import sys

sys.path.insert(0, str(Path(workflow.basedir)))
from chromatintrek import utils
from chromatintrek.genome import from_config

# ── Load config & resolve genome ────────────────────────────────────────────
CFG   = config                                         # populated via --configfile
GCFG  = from_config(CFG)
SAMPLES_DF = utils.parse_sample_sheet(CFG["sample_sheet"])
SAMPLES    = SAMPLES_DF["sample"].tolist()
COMPARISONS = CFG.get("comparisons", [])

# ── Directories ─────────────────────────────────────────────────────────────
TRIM_DIR    = "data/trimmed"
ALIGN_DIR   = "data/aligned"
PEAK_DIR    = "results/peaks"
BW_DIR      = "results/bigwig"
MATRIX_DIR  = "results/matrix"
MOTIF_DIR   = "results/motifs"
GENE_DIR    = "results/peak_gene"
GO_DIR      = "results/go"
QC_DIR      = "results/qc"
LOG_DIR     = "logs"


# ── Rules ────────────────────────────────────────────────────────────────────
include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/peaks.smk"
include: "rules/visualize.smk"
include: "rules/motif.smk"
include: "rules/peak_gene.smk"
include: "rules/go.smk"


# ── Target: all ─────────────────────────────────────────────────────────────
rule all:
    input:
        # MultiQC
        expand(f"{QC_DIR}/multiqc_report.html"),
        # Per-sample peaks
        expand(f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.narrowPeak", sample=SAMPLES),
        # BigWigs
        expand(f"{BW_DIR}/{{sample}}.bw", sample=SAMPLES),
        # Peak heatmaps (peak centered and TSS centered)
        expand(f"{MATRIX_DIR}/{{sample}}_heatmap.png", sample=SAMPLES),
        expand(f"{MATRIX_DIR}/{{sample}}_profile.png", sample=SAMPLES),
        expand(f"{MATRIX_DIR}/{{sample}}_tss_heatmap.png", sample=SAMPLES),
        expand(f"{MATRIX_DIR}/{{sample}}_tss_profile.png", sample=SAMPLES),
        # Motif enrichment
        expand(f"{MOTIF_DIR}/{{sample}}/knownResults.txt", sample=SAMPLES),
        # Peak-gene annotation
        expand(f"{GENE_DIR}/{{sample}}_annotated.tsv", sample=SAMPLES),
        # Gene Ontology
        expand(f"{GO_DIR}/{{sample}}/go_enrichment.tsv", sample=SAMPLES),
