"""Validation and Cell Ranger file contracts for optional ambient correction."""

import shutil
import subprocess
from pathlib import Path


MATRIX_FILES = ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
CELLBENDER_SETTINGS = {
    "default": {"epochs": 150, "learning_rate": 1e-4},
    "lower-lr": {"epochs": 300, "learning_rate": 5e-5},
}
# Default --tmp-dir values set by src/utils.py's default_tmp_dir(): Biowulf
# builds a node-local scratch path from $SLURM_JOBID, while Skyline's
# scratch is already exported per-node via $TMPDIR.
CELLBENDER_SLURM_TMP_DIRS = ("/lscratch/$SLURM_JOBID/", "$TMPDIR")


def gpu_available():
    """Checks whether a usable NVIDIA GPU is visible on this node by
    shelling out to `nvidia-smi`. Used to gate CellBender (which requires
    a GPU) under --mode local, where there is no job scheduler to request
    one from.
    @return <bool>: True if `nvidia-smi` is on PATH and reports at least
        one GPU, False otherwise (including if nvidia-smi is missing,
        times out, or errors, e.g. no driver loaded).
    """
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return False
    try:
        result = subprocess.run(
            [nvidia_smi, "--query-gpu=name", "--format=csv,noheader"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0 and bool(result.stdout.decode().strip())


def validate_ambient_options(options):
    method = options.get("ambient_rna_method", "none")
    if method not in ("none", "soupx", "cellbender"):
        raise ValueError("ambient_rna_method must be 'none', 'soupx', or 'cellbender'.")
    settings = options.get("cellbender_settings", "default")
    if settings not in CELLBENDER_SETTINGS:
        raise ValueError("cellbender_settings must be 'default' or 'lower-lr'.")
    if "cellbender_settings" in options and method != "cellbender":
        raise ValueError("--cellbender-settings requires --ambient-rna-method cellbender.")
    if method in ("soupx", "cellbender"):
        name = "SoupX" if method == "soupx" else "CellBender"
        if options["pipeline"] != "gex":
            raise ValueError(f"{name} requires --pipeline gex (gene expression only).")
        if options.get("aggregate", "") != "":
            raise ValueError(f"{name} cannot use Cell Ranger aggregation; omit --aggregate.")
    if method == "cellbender":
        if options["mode"] not in ("slurm", "local"):
            raise ValueError("CellBender requires --mode slurm or --mode local.")
        if options["mode"] == "slurm" and options["tmp_dir"] not in CELLBENDER_SLURM_TMP_DIRS:
            raise ValueError("CellBender requires the default scratch configuration; omit --tmp-dir.")
        if options["mode"] == "local" and not gpu_available():
            raise ValueError(
                "CellBender requires a GPU, but none was detected via nvidia-smi. "
                "--mode local runs on this node directly (no job scheduler to request "
                "a GPU from), so one must already be visible here. Run on a node with "
                "a GPU attached, or use --mode slurm to have CellBender's GPU requested "
                "through the cluster's job scheduler."
            )
    return method


def soupx_inputs(outs):
    outs = Path(outs)
    # These are the two GEX graph-clustering layouts supported by SoupX 1.6.2.
    modern = outs / "analysis" / "clustering" / "gene_expression_graphclust"
    legacy = outs / "analysis" / "clustering" / "graphclust"
    # Partial outputs must match the modern paths declared by the count rule.
    # Preserve SoupX's modern-directory precedence when both layouts exist.
    cluster_dir = modern
    if not modern.is_dir() and (legacy / "clusters.csv").is_file():
        cluster_dir = legacy
    return [
        str(outs / matrix / name)
        for matrix in ("raw_feature_bc_matrix", "filtered_feature_bc_matrix")
        for name in MATRIX_FILES
    ] + [str(cluster_dir / "clusters.csv")]


def validate_soupx_import(sample_dir):
    sample_dir = Path(sample_dir)
    missing = [p for p in soupx_inputs(sample_dir / "outs")
               if not Path(p).is_file()]
    if missing:
        raise ValueError("SoupX sample '{}': missing required input(s): {}".format(
            sample_dir.name, ", ".join(missing)
        ))


def cellbender_input(outs):
    return str(Path(outs) / "raw_feature_bc_matrix.h5")


def validate_cellbender_import(sample_dir):
    sample_dir = Path(sample_dir)
    raw = cellbender_input(sample_dir / "outs")
    if not Path(raw).is_file():
        raise ValueError(f"CellBender sample '{sample_dir.name}': missing required input: {raw}")
