"""Validation and Cell Ranger file contracts for optional ambient correction."""

from pathlib import Path


MATRIX_FILES = ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")


def validate_ambient_options(options):
    method = options.get("ambient_rna_method", "none")
    if method not in ("none", "soupx", "cellbender"):
        raise ValueError("ambient_rna_method must be 'none', 'soupx', or 'cellbender'.")
    if method in ("soupx", "cellbender"):
        name = "SoupX" if method == "soupx" else "CellBender"
        if options["pipeline"] != "gex":
            raise ValueError(f"{name} requires --pipeline gex (gene expression only).")
        if options.get("aggregate", "") != "":
            raise ValueError(f"{name} cannot use Cell Ranger aggregation; omit --aggregate.")
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
