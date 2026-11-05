#!/usr/bin/env python3
"""Export a CellBender 0.3.0 filtered v3 H5 without changing its counts/cell calls."""

import argparse
import gzip
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import scipy
from scipy.io import mmwrite
from scipy.sparse import csc_matrix


def convert(h5_path, output_dir, summary_path, sample, raw_input,
            settings="default", epochs=150, learning_rate=1e-4):
    # Read the genes-by-cells sparse matrix and its original annotations.
    with h5py.File(h5_path, "r") as handle:
        group = handle["matrix"]
        matrix = csc_matrix(
            (group["data"][:], group["indices"][:], group["indptr"][:]),
            shape=tuple(group["shape"][:]),
        )
        barcodes = group["barcodes"].asstr()[:]
        features = [group["features"][key].asstr()[:] for key in ("id", "name", "feature_type")]

    # Let SciPy validate the sparse structure; check counts and annotations here.
    matrix.check_format(full_check=True)
    if (min(matrix.shape) == 0 or matrix.dtype.kind not in "iuf"
            or not np.all(np.isfinite(matrix.data)) or np.any(matrix.data < 0)):
        raise ValueError("Expected a nonempty matrix of finite, nonnegative counts")
    for values, length in zip([barcodes, *features], [matrix.shape[1]] + [matrix.shape[0]] * 3):
        if values.shape != (length,):
            raise ValueError("Annotation lengths must match the matrix dimensions")
        if any(not value or any(c in value for c in "\t\r\n") for value in values):
            raise ValueError("Empty or invalid TSV annotation")

    # Write the three Cell Ranger files directly, without changing counts or names.
    output_dir, summary_path = Path(output_dir), Path(summary_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_dir / "matrix.mtx.gz", "wb") as stream:
        mmwrite(stream, matrix, field="integer" if matrix.dtype.kind in "iu" else "real",
                precision=17, symmetry="general")
    with gzip.open(output_dir / "barcodes.tsv.gz", "wt", encoding="utf-8") as stream:
        for barcode in barcodes:
            stream.write(barcode + "\n")
    with gzip.open(output_dir / "features.tsv.gz", "wt", encoding="utf-8") as stream:
        for annotation in zip(*features):
            stream.write("\t".join(annotation) + "\n")

    # Write the conversion record only after all three files finish successfully.
    record = {
        "state": "matrix_export_complete", "sample": sample,
        "raw_input": str(raw_input), "filtered_h5": str(h5_path),
        "matrix_dir": str(output_dir), "genes": matrix.shape[0],
        "cells": matrix.shape[1], "stored_entries": matrix.nnz,
        "module": "cellbender/0.3.0", "cuda": True, "fpr": 0.01,
        "cellbender_settings": settings, "epochs": epochs, "learning_rate": learning_rate,
        "expected_cells": "automatic", "total_droplets_included": "automatic",
        "python_executable": sys.executable, "python_version": sys.version,
        "h5py": h5py.__version__, "numpy": np.__version__, "scipy": scipy.__version__,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w", encoding="utf-8") as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="CellBender filtered v3 H5")
    parser.add_argument("--output", required=True, help="Cell Ranger MEX directory")
    parser.add_argument("--summary_path", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--raw-input", required=True)
    parser.add_argument("--settings", choices=["default", "lower-lr"], default="default")
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    args = parser.parse_args()
    convert(args.input, args.output, args.summary_path, args.sample, args.raw_input,
            args.settings, args.epochs, args.learning_rate)


if __name__ == "__main__":
    main()
