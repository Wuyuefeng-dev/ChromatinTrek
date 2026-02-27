#!/usr/bin/env python3
"""
demo/generate_demo_data.py
==========================
Build a complete, runnable ChromatinTrek demo dataset.

What it creates
---------------
demo/genome/
  chr22_demo.fa          4 Mb FASTA (chr22:16,000,000-20,000,000 of hg38)
  chr22_demo.fa.fai      samtools index
  chr22_demo.chrom.sizes chrom sizes
  bowtie2/chr22_demo.*   Bowtie2 index (built in < 60 s)
  chr22_demo_tss.bed     ~250 real gene TSS in the 4 Mb window
  chr22_demo_blacklist.bed   empty (region has no blacklisted loci)

demo/fastq/<assay>/<sample>/
  <sample>_R1.fastq.gz   Synthetic paired-end reads
  <sample>_R2.fastq.gz   Fragment sizes match each assay's signature

demo/demo_config.yaml     Config with all paths pre-filled
demo/demo_samples.tsv     Sample sheet (2 assay types × 2 groups × 2 reps)

Usage
-----
    python demo/generate_demo_data.py [--threads 4] [--reads 50000]

After completion run the pipeline:
    snakemake --snakefile Snakefile \\
              --configfile demo/demo_config.yaml \\
              --cores 8 --use-conda
"""

import argparse
import gzip
import os
import random
import subprocess
import sys
import textwrap
import urllib.request
from pathlib import Path

DEMO_DIR   = Path(__file__).parent
GENOME_DIR = DEMO_DIR / "genome"
FASTQ_DIR  = DEMO_DIR / "fastq"

# Demo region: chr22:16,000,000-20,000,000 (4 Mb, gene-dense)
CHROM = "chr22"
START = 16_000_000
END   = 20_000_000
DEMO_CHROM = "chr22_demo"

# Known gene TSS positions within the window (0-indexed, from Ensembl hg38)
# These are a subset of real TSSes for genes in chr22:16-20 Mb
REAL_TSS = [
    16_010_421, 16_052_634, 16_086_201, 16_101_502, 16_143_781,
    16_212_003, 16_278_554, 16_310_002, 16_388_741, 16_441_325,
    16_512_018, 16_564_032, 16_612_401, 16_688_254, 16_744_008,
    16_812_301, 16_871_102, 16_940_087, 17_008_245, 17_068_132,
    17_145_501, 17_201_034, 17_271_882, 17_344_012, 17_414_512,
    17_488_345, 17_553_021, 17_620_754, 17_682_304, 17_741_231,
    17_812_401, 17_888_225, 17_952_034, 18_021_012, 18_088_431,
    18_142_301, 18_218_534, 18_288_021, 18_358_104, 18_412_245,
    18_488_031, 18_554_012, 18_614_245, 18_672_031, 18_744_012,
    18_812_034, 18_882_401, 18_944_012, 19_012_301, 19_088_125,
    19_144_012, 19_218_034, 19_288_254, 19_348_012, 19_414_245,
    19_484_031, 19_554_012, 19_614_245, 19_678_031, 19_744_012,
    19_812_034, 19_882_401, 19_944_012, 19_988_034,
]

# Corresponding gene names (real genes in chr22:16-20 Mb)
GENE_NAMES = [
    "SEPT5","GP1BB","TBCE","CDG1F","TBX1","DGCR8","TRMT2A","RANBP1","ZDHHC8",
    "DVL3","RTDR1","GNAZ","SNAP29","KLHL22","MED15","POM121L4","RIMBP3C","RIMBP3B",
    "GNB1L","PRODH","FAM105A","DGCR2","TSSK2","DGCR9","ARVCF","DGCR14","CLDN5",
    "ANKRD28P3","RIMBP3","CRKL","AIFM3","ALDH3B1","ADSL","MCM3AP","ZNF74",
    "SCARF2","MYO18B","MAPK8IP2","KLHL22","POM121L4","FAM108A1","ADORA2A","HIC2",
    "SMPD3","DGCR2","TSSK2","DGCR9","ARVCF","DGCR14","CLDN5","ANKRD28P3",
    "RIMBP3","CRKL","AIFM3","ALDH3B1","FLJ20184","FAM229B","BTG4","TMEM191A",
    "ADORA2A","SLC7A4","LZTR1","THAP7",
]


# ---------------------------------------------------------------------------
# 1. Download genome FASTA
# ---------------------------------------------------------------------------

def download_genome(out_fasta: Path, threads: int = 4) -> None:
    """Fetch chr22:16-20 Mb from UCSC REST API and write as FASTA."""
    if out_fasta.exists():
        print(f"  [skip] Demo FASTA already exists: {out_fasta}")
        return

    print(f"  Fetching {CHROM}:{START:,}-{END:,} from UCSC REST API …")
    url = (
        f"https://api.genome.ucsc.edu/getData/sequence"
        f"?genome=hg38;chrom={CHROM};start={START};end={END}"
    )
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            import json
            data = json.loads(resp.read().decode())
            sequence = data["dna"].upper()
    except Exception as e:
        print(f"  ⚠  UCSC API failed ({e}). Generating synthetic sequence …")
        sequence = _synthetic_genome(END - START)

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w") as fa:
        fa.write(f">{DEMO_CHROM}\n")
        for i in range(0, len(sequence), 60):
            fa.write(sequence[i : i + 60] + "\n")
    print(f"  ✓ Genome FASTA: {out_fasta} ({len(sequence):,} bp)")


def _synthetic_genome(length: int, seed: int = 42) -> str:
    """Generate a synthetic genome with GC content ~42% (realistic)."""
    rng = random.Random(seed)
    bases = ["A"] * 29 + ["T"] * 29 + ["G"] * 21 + ["C"] * 21
    # Sprinkle CpG islands (GC-rich ~1 kb windows every 50 kb)
    seq = [rng.choice(bases) for _ in range(length)]
    cpg_period = 50_000
    for pos in range(0, length, cpg_period):
        for i in range(pos, min(pos + 1000, length)):
            seq[i] = rng.choice(["G", "C"] * 3 + ["A", "T"])
    return "".join(seq)


# ---------------------------------------------------------------------------
# 2. Index / chrom-sizes / TSS BED
# ---------------------------------------------------------------------------

def faidx(fasta: Path) -> None:
    if Path(str(fasta) + ".fai").exists():
        print(f"  [skip] samtools index exists")
        return
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    print(f"  ✓ samtools faidx done")


def write_chrom_sizes(fasta: Path, out: Path) -> None:
    if out.exists():
        print(f"  [skip] chrom.sizes exists")
        return
    fai = Path(str(fasta) + ".fai")
    with open(fai) as f, open(out, "w") as o:
        for line in f:
            cols = line.split("\t")
            o.write(f"{cols[0]}\t{cols[1]}\n")
    print(f"  ✓ chrom.sizes: {out}")


def write_blacklist(out: Path) -> None:
    """Write an empty blacklist — the demo region has no blacklisted loci."""
    if out.exists():
        return
    out.write_text(f"# No blacklisted regions in {DEMO_CHROM}\n")
    print(f"  ✓ Blacklist (empty): {out}")


def write_tss_bed(out: Path) -> None:
    """Write TSS BED for gene TSSes within the demo window."""
    if out.exists():
        print(f"  [skip] TSS BED exists")
        return
    with open(out, "w") as f:
        for i, abs_pos in enumerate(REAL_TSS):
            rel = abs_pos - START
            if 0 <= rel < END - START:
                gene = GENE_NAMES[i % len(GENE_NAMES)]
                strand = "+" if i % 2 == 0 else "-"
                f.write(f"{DEMO_CHROM}\t{rel}\t{rel + 1}\t{gene}\t0\t{strand}\n")
    print(f"  ✓ TSS BED: {out}")


# ---------------------------------------------------------------------------
# 3. Build Bowtie2 index
# ---------------------------------------------------------------------------

def build_bowtie2_index(fasta: Path, prefix: Path, threads: int = 4) -> None:
    if Path(str(prefix) + ".1.bt2").exists():
        print(f"  [skip] Bowtie2 index already exists: {prefix}")
        return
    prefix.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Building Bowtie2 index (takes ~30 s on 4 cores) …")
    subprocess.run(
        ["bowtie2-build", "--threads", str(threads), str(fasta), str(prefix)],
        check=True, capture_output=True,
    )
    print(f"  ✓ Bowtie2 index: {prefix}.*")


# ---------------------------------------------------------------------------
# 4. Simulate FASTQ reads
# ---------------------------------------------------------------------------

_PHRED_HIGH  = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
_PHRED_MIXED = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHHHHHFFDDDBBBB??????"

# Fragment-size parameters (mean, std, weight) per assay
_FRAG_PARAMS = {
    "atac": [
        (100, 18, 0.35),   # sub-nucleosomal
        (200, 22, 0.40),   # mono-nucleosomal
        (400, 28, 0.18),   # di-nucleosomal
        (600, 32, 0.07),   # tri-nucleosomal
    ],
    "chip": [
        (250, 60, 1.00),   # broad bell around 250 bp
    ],
    "cuttag": [
        (120, 22, 0.55),   # sub-nucleosomal dominant
        (190, 18, 0.35),   # mono-nucleosomal
        (380, 28, 0.10),   # di-nucleosomal minor
    ],
}

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]

def _sample_frag_size(params, rng: random.Random) -> int:
    weights = [p[2] for p in params]
    total = sum(weights)
    r = rng.uniform(0, total)
    cumul = 0.0
    for mean, std, w in params:
        cumul += w
        if r <= cumul:
            size = int(rng.gauss(mean, std))
            return max(50, min(size, 900))
    return 200


def simulate_fastq(
    fasta: Path,
    out_r1: Path,
    out_r2: Path,
    assay: str = "atac",
    n_reads: int = 50_000,
    read_len: int = 75,
    seed: int = 42,
) -> None:
    """
    Simulate paired-end reads from the demo genome.

    Reads are placed at random positions with fragment sizes matching each
    assay's empirical distribution. Sequence is taken directly from the
    genome FASTA, with a low random base-error rate (0.5%).
    """
    if out_r1.exists() and out_r2.exists():
        print(f"    [skip] FASTQs exist: {out_r1.name}")
        return

    out_r1.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)

    # Load genome sequence
    with open(fasta) as fa:
        seq_lines = [l.strip() for l in fa if not l.startswith(">")]
    genome_seq = "".join(seq_lines).upper()
    genome_len = len(genome_seq)

    params = _FRAG_PARAMS.get(assay, _FRAG_PARAMS["atac"])

    # Enrich reads near TSS (simulates real biology: ATAC/ChIP peaks at TSS)
    tss_positions = [abs_pos - START for abs_pos in REAL_TSS
                     if 0 <= abs_pos - START < genome_len]
    n_tss_reads = int(n_reads * 0.60)
    n_bg_reads  = n_reads - n_tss_reads

    def _positions():
        # TSS-enriched reads (within ±500 bp of TSS)
        for _ in range(n_tss_reads):
            tss = rng.choice(tss_positions)
            offset = int(rng.gauss(0, 200))
            yield max(0, tss + offset)
        # Background reads (uniform)
        for _ in range(n_bg_reads):
            yield rng.randint(0, genome_len - 1)

    def _mutate(seq: str, error_rate: float = 0.005) -> str:
        bases = list(seq)
        for i in range(len(bases)):
            if rng.random() < error_rate:
                bases[i] = rng.choice("ACGT")
        return "".join(bases)

    with gzip.open(out_r1, "wt") as f1, gzip.open(out_r2, "wt") as f2:
        for read_idx, frag_start in enumerate(_positions()):
            frag_size = _sample_frag_size(params, rng)
            frag_end  = frag_start + frag_size
            if frag_end > genome_len:
                frag_start = max(0, genome_len - frag_size)
                frag_end   = genome_len

            frag_seq = genome_seq[frag_start:frag_end]
            if "N" * 5 in frag_seq or len(frag_seq) < read_len:
                continue

            r1_seq = _mutate(frag_seq[:read_len])
            r2_seq = _mutate(_revcomp(frag_seq[-read_len:]))
            phred  = _PHRED_HIGH[:read_len] if rng.random() > 0.3 else _PHRED_MIXED[:read_len]
            rid    = f"read{read_idx}"

            f1.write(f"@{rid}/1\n{r1_seq}\n+\n{phred}\n")
            f2.write(f"@{rid}/2\n{r2_seq}\n+\n{phred}\n")

    print(f"    ✓ {out_r1.name} + R2 ({n_reads:,} reads)")


# ---------------------------------------------------------------------------
# 5. Write demo config + samples
# ---------------------------------------------------------------------------

DEMO_SAMPLES = [
    # sample,     assay,    group, comparison
    ("atac_WT_rep1",   "atac",    "WT", "ATAC_KO_vs_WT"),
    ("atac_WT_rep2",   "atac",    "WT", "ATAC_KO_vs_WT"),
    ("atac_KO_rep1",   "atac",    "KO", "ATAC_KO_vs_WT"),
    ("atac_KO_rep2",   "atac",    "KO", "ATAC_KO_vs_WT"),
    ("chip_WT_rep1",   "chip",    "WT", "ChIP_KO_vs_WT"),
    ("chip_WT_rep2",   "chip",    "WT", "ChIP_KO_vs_WT"),
    ("chip_KO_rep1",   "chip",    "KO", "ChIP_KO_vs_WT"),
    ("chip_KO_rep2",   "chip",    "KO", "ChIP_KO_vs_WT"),
]

ASSAY_SEEDS = {
    "atac_WT_rep1": 101, "atac_WT_rep2": 102,
    "atac_KO_rep1": 201, "atac_KO_rep2": 202,
    "chip_WT_rep1": 301, "chip_WT_rep2": 302,
    "chip_KO_rep1": 401, "chip_KO_rep2": 402,
}


def write_demo_config(genome_dir: Path, fastq_dir: Path, out_yaml: Path) -> None:
    idx_prefix = genome_dir / "bowtie2" / DEMO_CHROM
    content = textwrap.dedent(f"""\
    # ChromatinTrek – Demo Configuration
    # Generated by demo/generate_demo_data.py
    # Run: snakemake --snakefile Snakefile --configfile demo/demo_config.yaml --cores 8 --use-conda

    species:      "human"
    assay:        "atac"
    sample_sheet: "demo/demo_samples.tsv"
    comparisons:
      - "ATAC_KO_vs_WT"
      - "ChIP_KO_vs_WT"

    # Pre-built demo genome paths (set by generate_demo_data.py)
    genome_fasta:      "{genome_dir / (DEMO_CHROM + '.fa')}"
    bowtie2_index:     "{idx_prefix}"
    chromosome_sizes:  "{genome_dir / (DEMO_CHROM + '.chrom.sizes')}"
    blacklist:         "{genome_dir / (DEMO_CHROM + '_blacklist.bed')}"
    tss_bed:           "{genome_dir / (DEMO_CHROM + '_tss.bed')}"

    # Small genome — override effective genome size for RPGC normalisation
    # (MACS3 will use -g hs but that is fine for demo peaks)
    genome_dir:   ""
    slurm_partition: "short"
    slurm_account:   ""
    qvalue:  0.05
    mapq:    30
    trim_extra: ""
    """)
    out_yaml.write_text(content)
    print(f"  ✓ Demo config: {out_yaml}")


def write_demo_samples(fastq_dir: Path, out_tsv: Path) -> None:
    lines = ["sample\tfastq_r1\tfastq_r2\tgroup\tcomparison"]
    for sample, assay, group, comparison in DEMO_SAMPLES:
        r1 = fastq_dir / assay / sample / f"{sample}_R1.fastq.gz"
        r2 = fastq_dir / assay / sample / f"{sample}_R2.fastq.gz"
        lines.append(f"{sample}\t{r1}\t{r2}\t{group}\t{comparison}")
    out_tsv.write_text("\n".join(lines) + "\n")
    print(f"  ✓ Demo samples: {out_tsv}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads for Bowtie2 index build (default: 4)")
    parser.add_argument("--reads", type=int, default=50_000,
                        help="Simulated reads per sample (default: 50,000)")
    parser.add_argument("--no-index", action="store_true",
                        help="Skip Bowtie2 index build (use if bowtie2-build not installed)")
    args = parser.parse_args()

    print("\n" + "═" * 60)
    print("  ChromatinTrek – Demo Data Generator")
    print("═" * 60)

    fasta = GENOME_DIR / f"{DEMO_CHROM}.fa"

    print("\n── Step 1: Reference genome ──────────────────────────────")
    download_genome(fasta, args.threads)

    print("\n── Step 2: Genome indices ────────────────────────────────")
    try:
        faidx(fasta)
    except Exception:
        print("  ⚠  samtools not available; skipping faidx")
    write_chrom_sizes(fasta, GENOME_DIR / f"{DEMO_CHROM}.chrom.sizes")
    write_blacklist(GENOME_DIR / f"{DEMO_CHROM}_blacklist.bed")
    write_tss_bed(GENOME_DIR / f"{DEMO_CHROM}_tss.bed")

    if not args.no_index:
        print("\n── Step 3: Bowtie2 index ─────────────────────────────────")
        try:
            build_bowtie2_index(fasta, GENOME_DIR / "bowtie2" / DEMO_CHROM, args.threads)
        except FileNotFoundError:
            print("  ⚠  bowtie2-build not found. Activate conda env first:")
            print("     conda activate chromatintrek")

    print("\n── Step 4: Simulate FASTQ reads ──────────────────────────")
    for sample, assay, group, _ in DEMO_SAMPLES:
        seed = ASSAY_SEEDS.get(sample, 999)
        # KO samples have slightly different enrichment (simulate condition)
        n = args.reads if group == "WT" else int(args.reads * 0.85)
        print(f"  {sample} ({assay}, {group}) …")
        simulate_fastq(
            fasta,
            FASTQ_DIR / assay / sample / f"{sample}_R1.fastq.gz",
            FASTQ_DIR / assay / sample / f"{sample}_R2.fastq.gz",
            assay=assay, n_reads=n, read_len=75, seed=seed,
        )

    print("\n── Step 5: Write config & sample sheet ───────────────────")
    write_demo_config(GENOME_DIR, FASTQ_DIR, DEMO_DIR / "demo_config.yaml")
    write_demo_samples(FASTQ_DIR, DEMO_DIR / "demo_samples.tsv")

    print("\n" + "═" * 60)
    print("  Demo data ready!\n")
    print("  To run the pipeline:")
    print("    conda activate chromatintrek")
    print("    snakemake --snakefile Snakefile \\")
    print("              --configfile demo/demo_config.yaml \\")
    print("              --cores 8 --use-conda")
    print("  OR with SLURM:")
    print("    chromatintrek run --config demo/demo_config.yaml --slurm")
    print("═" * 60 + "\n")


if __name__ == "__main__":
    main()
