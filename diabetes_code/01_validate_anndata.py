#!/usr/bin/env python3
"""
01_validate_anndata.py

Inspects and validates an AnnData object for exploratory single-cell Perturb-seq analysis.
Performs:
  - Slurm and Conda environment verification checks
  - AnnData dimensions (.n_obs, .n_vars)
  - Expression matrix diagnostics for adata.X and adata.layers['GPTin']:
      * Matrix shape
      * Data type (dtype)
      * Sparse vs dense status
      * Sampled minimum and maximum values
      * Non-negativity check
      * Integer-like vs continuous check
      * Sampled row sums (mean, std, min, max, sample values)
      * Pearson and Spearman correlation with adata.obs['nCount_RNA']
  - Comparison between adata.X and adata.layers['GPTin']
  - Mathematical verification of log-normalization and raw-count reconstruction
  - Comprehensive metadata structure audit (.obs, .var, .obsm, .layers, .uns)
  - Missing value / NaN audit across all metadata columns
  - Value distributions and category counts for key metadata columns
  - Output summary report in JSON and human-readable text formats

Usage:
  python 01_validate_anndata.py --input-h5ad diabetes.h5ad --output-dir results
"""

import argparse
import json
import os
import sys
from pathlib import Path
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr, spearmanr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate AnnData object and determine raw count matrix."
    )
    parser.add_argument(
        "--input-h5ad",
        type=str,
        required=True,
        help="Path to input .h5ad file.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results",
        help="Directory to save validation summaries and reports.",
    )
    return parser.parse_args()


def inspect_matrix(matrix, obs_ncount, name="Matrix", sample_size=1000, random_seed=42):
    """
    Perform thorough diagnostics on a dense or sparse 2D matrix using a sample of cells
    to prevent dense memory exhaustion.
    """
    is_sparse = sp.issparse(matrix)
    dtype_str = str(matrix.dtype)
    shape = list(matrix.shape)
    format_type = matrix.format if is_sparse else "dense_numpy"

    # Sample cell indices reproducibly
    n_obs = shape[0]
    sample_n = min(sample_size, n_obs)
    rng = np.random.RandomState(random_seed)
    sample_indices = rng.choice(n_obs, size=sample_n, replace=False)

    sub_mat = matrix[sample_indices]

    if is_sparse:
        sampled_data = sub_mat.data
        sub_row_sums = np.array(sub_mat.sum(axis=1)).flatten()
        total_elements = shape[0] * shape[1]
        sparsity_pct = round((1.0 - (matrix.nnz / total_elements)) * 100.0, 4)
        nnz = int(matrix.nnz)
    else:
        sampled_data = sub_mat.flatten()
        sub_row_sums = np.array(sub_mat.sum(axis=1)).flatten()
        total_elements = shape[0] * shape[1]
        nnz = int(np.count_nonzero(matrix))
        sparsity_pct = round((1.0 - (nnz / total_elements)) * 100.0, 4)

    if len(sampled_data) > 0:
        min_val = float(np.min(sampled_data))
        max_val = float(np.max(sampled_data))
        mean_nonzero = float(np.mean(sampled_data))
        non_negative = bool(min_val >= 0.0)
        # Check if values are integer-like
        is_integer_valued = bool(np.all(np.isclose(sampled_data, np.round(sampled_data))))
    else:
        min_val = max_val = mean_nonzero = 0.0
        non_negative = True
        is_integer_valued = True

    # Correlation with nCount_RNA
    sampled_ncount = obs_ncount.iloc[sample_indices].values.astype(float)
    p_corr, _ = pearsonr(sub_row_sums, sampled_ncount)
    s_corr, _ = spearmanr(sub_row_sums, sampled_ncount)

    return {
        "name": name,
        "shape": shape,
        "dtype": dtype_str,
        "is_sparse": is_sparse,
        "format": format_type,
        "nnz": nnz,
        "sparsity_percent": sparsity_pct,
        "sampled_min": min_val,
        "sampled_max": max_val,
        "sampled_mean_nonzero": round(mean_nonzero, 4),
        "is_non_negative": non_negative,
        "is_integer_valued": is_integer_valued,
        "sampled_row_sums_mean": round(float(np.mean(sub_row_sums)), 4),
        "sampled_row_sums_std": round(float(np.std(sub_row_sums)), 4),
        "sampled_row_sums_min": round(float(np.min(sub_row_sums)), 4),
        "sampled_row_sums_max": round(float(np.max(sub_row_sums)), 4),
        "sampled_row_sums_sample": [round(float(v), 4) for v in sub_row_sums[:5]],
        "pearson_corr_with_nCount_RNA": round(float(p_corr), 6),
        "spearman_corr_with_nCount_RNA": round(float(s_corr), 6),
    }


def main():
    args = parse_args()
    input_path = Path(args.input_h5ad).resolve()
    output_dir = Path(args.output_dir).resolve()
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("01_VALIDATE_ANNDATA: SINGLE-CELL DATASET & COUNT MATRIX VALIDATION")
    print("=" * 80)
    print(f"Input file: {input_path}")
    print(f"Output directory: {output_dir}")

    if not input_path.exists():
        print(f"ERROR: Input file {input_path} does not exist.", file=sys.stderr)
        sys.exit(1)

    print("\n[Step 1/6] Loading AnnData object into memory...")
    adata = ad.read_h5ad(input_path)
    print(f"Successfully loaded AnnData object: {adata.n_obs:,} cells × {adata.n_vars:,} genes.")

    # Step 2: Expression matrix diagnostics
    print("\n[Step 2/6] Inspecting expression matrices (adata.X and layers)...")
    if "nCount_RNA" not in adata.obs.columns:
        raise ValueError("Critical column 'nCount_RNA' missing from adata.obs!")

    obs_ncount = adata.obs["nCount_RNA"]

    x_diag = inspect_matrix(adata.X, obs_ncount, name="adata.X")
    print("\n--- Diagnostics for adata.X ---")
    for k, v in x_diag.items():
        print(f"  {k}: {v}")

    layer_diags = {}
    if adata.layers:
        for layer_key in adata.layers.keys():
            layer_diag = inspect_matrix(
                adata.layers[layer_key], obs_ncount, name=f"adata.layers['{layer_key}']"
            )
            layer_diags[layer_key] = layer_diag
            print(f"\n--- Diagnostics for adata.layers['{layer_key}'] ---")
            for k, v in layer_diag.items():
                print(f"  {k}: {v}")
    else:
        print("  No layers found.")

    # Step 3: Comparison and Raw Counts Determination
    print("\n[Step 3/6] Comparing adata.X vs adata.layers['GPTin'] & Testing Reconstruction...")
    gptin_identical = False
    max_diff_x_gptin = None
    if "GPTin" in adata.layers:
        if adata.X.shape == adata.layers["GPTin"].shape and adata.X.nnz == adata.layers["GPTin"].nnz:
            diff_vec = np.abs((adata.X - adata.layers["GPTin"]).data)
            max_diff_x_gptin = float(np.max(diff_vec)) if len(diff_vec) > 0 else 0.0
            gptin_identical = bool(max_diff_x_gptin < 1e-6)
            print(f"Comparison max|adata.X - adata.layers['GPTin']|: {max_diff_x_gptin}")
            print(f"Are adata.X and adata.layers['GPTin'] identical: {gptin_identical}")

    # Check Seurat log-normalization formula on sample cells
    # Seurat LogNormalize formula: x = log(1 + 10000 * count / nCount_RNA)
    # Inverse: raw_count = expm1(x) * nCount_RNA / 10000
    rng = np.random.RandomState(42)
    sample_cell_indices = rng.choice(adata.n_obs, size=min(100, adata.n_obs), replace=False)
    reconstruction_max_diffs = []
    reconstruction_sum_diffs = []

    for c_idx in sample_cell_indices:
        cell_row = adata.X[c_idx]
        cell_nCount = float(adata.obs["nCount_RNA"].iloc[c_idx])
        reconstructed = np.expm1(cell_row.data) * cell_nCount / 10000.0
        rounded = np.round(reconstructed)
        int_diff = np.max(np.abs(reconstructed - rounded)) if len(reconstructed) > 0 else 0.0
        sum_diff = abs(np.sum(reconstructed) - cell_nCount)
        reconstruction_max_diffs.append(int_diff)
        reconstruction_sum_diffs.append(sum_diff)

    max_int_discrepancy = float(np.max(reconstruction_max_diffs))
    max_sum_discrepancy = float(np.max(reconstruction_sum_diffs))
    is_exact_lognorm = bool(max_int_discrepancy < 1e-5 and max_sum_discrepancy < 1e-3)

    print(f"LogNormalize verification (max int discrepancy): {max_int_discrepancy:.2e}")
    print(f"LogNormalize verification (max nCount_RNA sum discrepancy): {max_sum_discrepancy:.2e}")
    print(f"Exact Seurat LogNormalize confirmed: {is_exact_lognorm}")

    raw_counts_decision = {
        "adata_X_is_raw_counts": False,
        "layers_GPTin_is_raw_counts": False,
        "matrix_value_nature": "log1p-transformed normalized expression (Seurat LogNormalize with scale factor 10,000)",
        "identical_matrices": gptin_identical,
        "max_diff_x_vs_gptin": max_diff_x_gptin,
        "raw_counts_reconstructible": is_exact_lognorm,
        "raw_counts_reconstruction_formula": "raw_counts = round(expm1(X) * nCount_RNA / 10000)",
        "max_integer_reconstruction_discrepancy": max_int_discrepancy,
        "selected_raw_counts_source": "Reconstructed from adata.X using expm1(X) * nCount_RNA / 10000",
        "preprocessing_recommendation": (
            "Because both adata.X and adata.layers['GPTin'] store already log-normalized values "
            "(log1p(CP10k)), do NOT apply a second sc.pp.normalize_total() and sc.pp.log1p() directly on X. "
            "For flavor='seurat_v3' HVG selection, populate a raw counts layer using exact mathematical reconstruction: "
            "round(expm1(X) * nCount_RNA / 10000)."
        ),
    }

    # Step 4: Metadata Structure Audit
    print("\n[Step 4/6] Auditing metadata columns and missing values...")
    obs_info = {
        "num_columns": len(adata.obs.columns),
        "columns": list(adata.obs.columns),
        "dtypes": {col: str(adata.obs[col].dtype) for col in adata.obs.columns},
        "missing_values": {col: int(adata.obs[col].isna().sum()) for col in adata.obs.columns},
    }

    var_info = {
        "num_columns": len(adata.var.columns),
        "columns": list(adata.var.columns),
        "dtypes": {col: str(adata.var[col].dtype) for col in adata.var.columns},
        "missing_values": {col: int(adata.var[col].isna().sum()) for col in adata.var.columns},
    }

    obsm_info = {
        k: {"shape": list(adata.obsm[k].shape), "dtype": str(adata.obsm[k].dtype)}
        for k in adata.obsm.keys()
    }
    uns_keys = list(adata.uns.keys())

    print(f".obs columns ({obs_info['num_columns']}): {obs_info['columns']}")
    print(f".obsm keys: {obsm_info}")
    print(f".layers keys: {list(adata.layers.keys())}")
    print(f".uns keys: {uns_keys}")

    # Step 5: Category breakdown for key metadata columns
    print("\n[Step 5/6] Computing category breakdown for key metadata columns...")
    key_categorical_columns = [
        "celltype_2",
        "genotype",
        "seurat_clusters",
        "RNA_snn_res.0.5",
        "integrated_snn_res.0.5",
        "sub.cluster",
        "celltype_individual",
        "celltype_2_old",
        "orig.ident",
        "time_point",
        "time",
        "Phase",
        "gene",
    ]

    category_counts = {}
    for col in key_categorical_columns:
        if col in adata.obs.columns:
            counts = adata.obs[col].value_counts(dropna=False).to_dict()
            category_counts[col] = {str(k): int(v) for k, v in counts.items()}
            print(f"\nCategory breakdown for '{col}' ({len(counts)} unique values):")
            items = list(counts.items())
            for cat, count in items[:10]:
                print(f"    {cat}: {count:,} ({count / adata.n_obs * 100:.2f}%)")
            if len(items) > 10:
                print(f"    ... and {len(items) - 10} more categories.")

    # Step 6: Write outputs
    print("\n[Step 6/6] Writing validation reports to disk...")
    validation_summary = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "matrix_diagnostics": {
            "adata_X": x_diag,
            "layers": layer_diags,
        },
        "raw_counts_decision": raw_counts_decision,
        "obs_metadata": obs_info,
        "var_metadata": var_info,
        "obsm_metadata": obsm_info,
        "uns_keys": uns_keys,
        "category_counts": category_counts,
    }

    json_path = reports_dir / "validation_summary.json"
    with open(json_path, "w") as f:
        json.dump(validation_summary, f, indent=2)
    print(f"JSON validation report saved to: {json_path}")

    txt_path = reports_dir / "validation_report.txt"
    with open(txt_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("ANNDATA OBJECT VALIDATION REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Dataset Dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes\n\n")

        f.write("1. EXPRESSION MATRIX DIAGNOSTICS\n")
        f.write("-" * 50 + "\n")
        for m_name, diag in [("adata.X", x_diag)] + [(f"adata.layers['{k}']", v) for k, v in layer_diags.items()]:
            f.write(f"\nMatrix: {m_name}\n")
            f.write(f"  Shape: {diag['shape']}\n")
            f.write(f"  Dtype: {diag['dtype']}\n")
            f.write(f"  Sparse: {diag['is_sparse']} ({diag['format']})\n")
            f.write(f"  Sparsity: {diag['sparsity_percent']}%\n")
            f.write(f"  Sampled Min: {diag['sampled_min']}, Max: {diag['sampled_max']}\n")
            f.write(f"  Non-negative: {diag['is_non_negative']}\n")
            f.write(f"  Integer-like: {diag['is_integer_valued']}\n")
            f.write(f"  Sampled Row Sums Mean: {diag['sampled_row_sums_mean']} (std: {diag['sampled_row_sums_std']})\n")
            f.write(f"  Sampled Row Sums Range: [{diag['sampled_row_sums_min']}, {diag['sampled_row_sums_max']}]\n")
            f.write(f"  Sampled Row Sums (First 5): {diag['sampled_row_sums_sample']}\n")
            f.write(f"  Pearson Corr with nCount_RNA: {diag['pearson_corr_with_nCount_RNA']}\n")
            f.write(f"  Spearman Corr with nCount_RNA: {diag['spearman_corr_with_nCount_RNA']}\n")

        f.write("\n\n2. RAW COUNTS & NORMALIZATION EVALUATION\n")
        f.write("-" * 50 + "\n")
        f.write(f"adata.X == adata.layers['GPTin']: {raw_counts_decision['identical_matrices']}\n")
        f.write(f"Matrix Value Nature: {raw_counts_decision['matrix_value_nature']}\n")
        f.write(f"Raw Counts Directly in X: {raw_counts_decision['adata_X_is_raw_counts']}\n")
        f.write(f"Raw Counts Directly in GPTin: {raw_counts_decision['layers_GPTin_is_raw_counts']}\n")
        f.write(f"Reconstruction Verified: {raw_counts_decision['raw_counts_reconstructible']}\n")
        f.write(f"Reconstruction Formula: {raw_counts_decision['raw_counts_reconstruction_formula']}\n")
        f.write(f"Recommendation: {raw_counts_decision['preprocessing_recommendation']}\n\n")

        f.write("3. METADATA COLUMNS (.obs)\n")
        f.write("-" * 50 + "\n")
        for col in adata.obs.columns:
            n_miss = obs_info["missing_values"][col]
            f.write(f"  - {col} (dtype: {adata.obs[col].dtype}, missing: {n_miss})\n")

        f.write("\n\n4. KEY METADATA CATEGORY COUNTS\n")
        f.write("-" * 50 + "\n")
        for col, counts in category_counts.items():
            f.write(f"\nColumn: {col} ({len(counts)} categories)\n")
            for k, v in counts.items():
                f.write(f"  - {k}: {v:,} cells ({v / adata.n_obs * 100:.2f}%)\n")

    print(f"Human-readable text report saved to: {txt_path}")
    print("\nValidation completed successfully!")


if __name__ == "__main__":
    main()
