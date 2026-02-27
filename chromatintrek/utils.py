"""
chromatintrek.utils
===================
Shared utilities: config loading, logging, subprocess wrapper, sample sheet parsing.
"""

import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def load_config(path: Union[str, Path]) -> Dict[str, Any]:
    """Load pipeline configuration from a YAML file."""
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    return cfg


def validate_config(cfg: Dict[str, Any]) -> None:
    """Raise ValueError for any missing required config keys."""
    required = ["sample_sheet", "genome", "bowtie2_index"]
    missing = [k for k in required if not cfg.get(k)]
    if missing:
        raise ValueError(f"Missing required config fields: {missing}")


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(
    log_dir: Union[str, Path],
    name: str = "chromatintrek",
    level: int = logging.INFO,
) -> logging.Logger:
    """Configure and return a logger writing to file and stdout."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{name}.log"

    logger = logging.getLogger(name)
    logger.setLevel(level)
    if not logger.handlers:
        fmt = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        fh = logging.FileHandler(log_file)
        fh.setFormatter(fmt)
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(fmt)
        logger.addHandler(fh)
        logger.addHandler(sh)
    return logger


# ---------------------------------------------------------------------------
# Subprocess helper
# ---------------------------------------------------------------------------

def run_cmd(
    cmd: Union[str, List[str]],
    logfile: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
    check: bool = True,
) -> Optional[subprocess.CompletedProcess]:
    """
    Execute a shell command, optionally logging stdout/stderr to *logfile*.

    Parameters
    ----------
    cmd      : Command string or list of tokens.
    logfile  : Path to capture stdout + stderr.
    dry_run  : Print command without executing when True.
    check    : Raise CalledProcessError on non-zero exit.
    """
    cmd_str = cmd if isinstance(cmd, str) else " ".join(str(c) for c in cmd)
    if dry_run:
        print(f"[DRY-RUN] {cmd_str}")
        return None
    if logfile:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        with open(logfile, "w") as log:
            result = subprocess.run(
                cmd_str, shell=True, stdout=log, stderr=subprocess.STDOUT, check=check
            )
    else:
        result = subprocess.run(cmd_str, shell=True, check=check)
    return result


# ---------------------------------------------------------------------------
# File system
# ---------------------------------------------------------------------------

def makedirs(*paths: Union[str, Path]) -> None:
    """Create directories (and parents) if they do not exist."""
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Tool availability checks
# ---------------------------------------------------------------------------

def check_tools(*tool_names: str) -> Dict[str, bool]:
    """Return dict {tool_name: is_on_PATH} for each requested tool."""
    return {name: shutil.which(name) is not None for name in tool_names}


def require_tools(*tool_names: str) -> None:
    """Raise RuntimeError listing any tools not found on PATH."""
    missing = [t for t, ok in check_tools(*tool_names).items() if not ok]
    if missing:
        raise RuntimeError(
            f"Required tools not found on PATH: {missing}\n"
            "Activate the conda environment:  conda activate chromatintrek"
        )


# ---------------------------------------------------------------------------
# Sample sheet
# ---------------------------------------------------------------------------

def parse_sample_sheet(path: Union[str, Path]) -> pd.DataFrame:
    """
    Parse a tab-separated sample sheet.

    Required columns : sample, fastq_r1
    Optional columns : fastq_r2, assay (atac|chip|cuttag), control, group

    Returns a DataFrame indexed by *sample*.
    """
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    required = {"sample", "fastq_r1"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Sample sheet missing required columns: {missing}")

    for col, default in [("fastq_r2", ""), ("assay", "atac"), ("control", ""), ("group", "")]:
        if col not in df.columns:
            df[col] = default
    df["fastq_r2"] = df["fastq_r2"].fillna("")
    return df.set_index("sample", drop=False)


def get_samples(config: Dict[str, Any]) -> List[str]:
    """Return list of sample names from config's sample sheet."""
    return parse_sample_sheet(config["sample_sheet"])["sample"].tolist()


def is_paired_end(sample: str, config: Dict[str, Any]) -> bool:
    """Return True if sample has a second FASTQ (paired-end)."""
    df = parse_sample_sheet(config["sample_sheet"])
    return bool(df.loc[sample, "fastq_r2"])


def get_assay(sample: str, config: Dict[str, Any]) -> str:
    """Return assay type (atac|chip|cuttag) for a sample."""
    df = parse_sample_sheet(config["sample_sheet"])
    return df.loc[sample, "assay"].lower()


def get_control(sample: str, config: Dict[str, Any]) -> Optional[str]:
    """Return control sample name, or None if not set."""
    df = parse_sample_sheet(config["sample_sheet"])
    ctrl = df.loc[sample, "control"]
    return ctrl if ctrl else None


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def out_path(config: Dict[str, Any], key: str, *parts: str) -> Path:
    """Construct an output path: config[key] / parts."""
    return Path(config[key]).joinpath(*parts)
