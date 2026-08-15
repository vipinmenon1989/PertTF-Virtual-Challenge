#!/usr/bin/env python3
"""
02_scanpy_preprocessing.py

Performs standard Scanpy preprocessing on the Perturb-seq AnnData dataset:
  1. Loads the AnnData object.
  2. Preserves the existing Seurat UMAP embedding as obsm['X_umap_existing'].
  3. Mathematically reconstructs the sparse integer count matrix:
       counts = round(expm1(X) * nCount_RNA / 10000)
     Stores it as adata.layers['counts'] and validates cell row sums against obs['nCount_RNA'].
  4. Identifies exactly 5,000 highly variable genes using:
       sc.pp.highly_variable_genes(
           adata,
           n_top_genes=5000,
           flavor="seurat_v3",
           layer="counts",
           batch_key="orig.ident",
           subset=False
       )
  5. Creates a working copy (adata_work) with raw counts, normalizes to 10,000, applies log1p,
     and subsets to the 5,000 HVGs.
  6. Computes PCA with 50 components (reproducible seed, svd_solver='arpack').
  7. Computes neighbor graph with n_neighbors=15, n_pcs=50.
  8. Computes new Scanpy UMAP and stores it as obsm['X_umap_scanpy'].
  9. Validates all embeddings, PCA coordinates, and output properties.
 10. Exports lightweight summaries (CSV/JSON/gzipped tables) without creating or committing .h5ad files.

Usage:
  python 02_scanpy_preprocessing.py --input-h5ad diabetes.h5ad --output-dir results
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
import scanpy as sc
import scipy.sparse as sp


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scanpy preprocessing with 5,000 HVGs, PCA, neighbors, and UMAP."
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
        help="Directory to save summaries and coordinates.",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=5000,
        help="Number of highly variable genes to select (default: 5000).",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of principal components (default: 50).",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for kNN graph (default: 15).",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42).",
    )
    return parser.parse_args()


def reconstruct_raw_counts(adata):
    """
    Reconstruct exact sparse integer counts from Seurat LogNormalize matrix:
      counts = round(expm1(X) * nCount_RNA / 10000)
    Uses CSR structure to compute efficiently without densifying memory.
    """
    print("Reconstructing sparse raw count matrix from log-normalized expression...")
    if not sp.issparse(adata.X):
        raise TypeError("Expected adata.X to be a scipy sparse matrix.")

    csr_x = adata.X.tocsr() if adata.X.format != "csr" else adata.X
    n_obs = adata.n_obs
    nCount = adata.obs["nCount_RNA"].values.astype(np.float64)

    # Compute row indices for each nonzero entry in CSR data
    row_indices = np.repeat(np.arange(n_obs, dtype=np.int32), np.diff(csr_x.indptr))
    cell_scales = (nCount / 10000.0)[row_indices]

    # Reconstruct integer counts
    raw_data = np.round(np.expm1(csr_x.data) * cell_scales)
    # Ensure non-negative
    raw_data = np.clip(raw_data, 0, None)

    counts_csr = sp.csr_matrix(
        (raw_data, csr_x.indices.copy(), csr_x.indptr.copy()),
        shape=csr_x.shape,
        dtype=np.float32,
    )
    return counts_csr


def validate_count_reconstruction(counts_csr, obs_ncount):
    """Validate that reconstructed row sums match nCount_RNA across all cells."""
    row_sums = np.array(counts_csr.sum(axis=1)).flatten()
    ncount_vals = obs_ncount.values.astype(np.float64)
    abs_diff = np.abs(row_sums - ncount_vals)

    max_diff = float(np.max(abs_diff))
    mean_diff = float(np.mean(abs_diff))
    median_diff = float(np.median(abs_diff))
    exact_matches = int(np.sum(abs_diff < 1e-4))
    total_cells = len(ncount_vals)
    pct_exact = round((exact_matches / total_cells) * 100.0, 4)

    validation_stats = {
        "total_cells": total_cells,
        "exact_matches_count": exact_matches,
        "exact_matches_percent": pct_exact,
        "max_absolute_difference": max_diff,
        "mean_absolute_difference": mean_diff,
        "median_absolute_difference": median_diff,
        "is_valid": bool(max_diff < 1.0 and pct_exact > 99.9),
    }
    return validation_stats


def main():
    args = parse_args()
    input_path = Path(args.input_h5ad).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("02_SCANPY_PREPROCESSING: HVG SELECTION, PCA, NEIGHBORS & UMAP")
    print("=" * 80)
    print(f"Input file: {input_path}")
    print(f"Output directory: {output_dir}")
    print(f"Parameters: HVGs={args.n_top_genes}, PCs={args.n_pcs}, Neighbors={args.n_neighbors}, Seed={args.random_state}")

    if not input_path.exists():
        print(f"ERROR: Input file {input_path} does not exist.", file=sys.stderr)
        sys.exit(1)

    # Step 1: Load AnnData
    print("\n[Step 1/8] Loading AnnData object...")
    adata = ad.read_h5ad(input_path)
    print(f"Loaded dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes.")

    # Step 2: Preserve existing UMAP
    print("\n[Step 2/8] Preserving existing UMAP embedding...")
    if "X_umap" in adata.obsm:
        adata.obsm["X_umap_existing"] = adata.obsm["X_umap"].copy()
        print(f"Preserved obsm['X_umap'] as obsm['X_umap_existing'] (shape: {adata.obsm['X_umap_existing'].shape}).")
    else:
        print("WARNING: 'X_umap' was not found in adata.obsm.")

    # Step 3: Reconstruct raw counts & validate
    print("\n[Step 3/8] Reconstructing raw integer count matrix...")
    counts_csr = reconstruct_raw_counts(adata)
    adata.layers["counts"] = counts_csr

    recon_stats = validate_count_reconstruction(counts_csr, adata.obs["nCount_RNA"])
    print("Count Reconstruction Validation:")
    for k, v in recon_stats.items():
        print(f"  {k}: {v}")

    if not recon_stats["is_valid"]:
        raise ValueError(f"Count reconstruction validation failed! Max diff: {recon_stats['max_absolute_difference']}")

    recon_json_path = output_dir / "count_reconstruction_validation.json"
    with open(recon_json_path, "w") as f:
        json.dump(recon_stats, f, indent=2)
    print(f"Saved count reconstruction validation to: {recon_json_path}")

    # Step 4: Highly Variable Genes selection
    print(f"\n[Step 4/8] Selecting {args.n_top_genes} highly variable genes (flavor='seurat_v3', batch_key='orig.ident')...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top_genes,
        flavor="seurat_v3",
        layer="counts",
        batch_key="orig.ident",
        subset=False,
    )

    n_hvg_selected = int(adata.var["highly_variable"].sum())
    print(f"Selected {n_hvg_selected:,} highly variable genes (expected {args.n_top_genes:,}).")
    if n_hvg_selected != args.n_top_genes:
        raise ValueError(f"Expected {args.n_top_genes} HVGs, but found {n_hvg_selected}!")

    # Export HVG summary table
    hvg_cols = [c for c in adata.var.columns if "highly_variable" in c or "variance" in c or "mean" in c or "rank" in c]
    hvg_df = adata.var[hvg_cols].copy()
    hvg_summary_path = output_dir / "hvg_summary.csv"
    hvg_df.to_csv(hvg_summary_path)
    print(f"Saved HVG summary table to: {hvg_summary_path}")

    # Step 5: Create working copy and normalize
    print("\n[Step 5/8] Creating working copy, normalizing to 10,000, and applying log1p...")
    adata_work = adata.copy()
    adata_work.X = adata_work.layers["counts"].copy()

    sc.pp.normalize_total(adata_work, target_sum=10000)
    sc.pp.log1p(adata_work)

    print(f"Subsetting working copy to {args.n_top_genes} HVGs...")
    adata_work = adata_work[:, adata_work.var["highly_variable"]].copy()
    print(f"Working copy shape: {adata_work.n_obs:,} cells × {adata_work.n_vars:,} HVGs.")

    # Step 6: PCA
    print(f"\n[Step 6/8] Running PCA (n_comps={args.n_pcs}, random_state={args.random_state})...")
    sc.pp.pca(
        adata_work,
        n_comps=args.n_pcs,
        zero_center=True,
        svd_solver="arpack",
        random_state=args.random_state,
    )

    pca_coords = adata_work.obsm["X_pca"]
    pca_has_nan = bool(np.isnan(pca_coords).any())
    print(f"PCA shape: {pca_coords.shape}, Has NaN: {pca_has_nan}")
    if pca_has_nan:
        raise ValueError("PCA coordinates contain NaN values!")

    variance_ratio = adata_work.uns["pca"]["variance_ratio"]
    total_var_explained = float(np.sum(variance_ratio))
    print(f"Top {args.n_pcs} PCs explain {total_var_explained * 100:.2f}% of total variance.")

    # Step 7: Neighbors
    print(f"\n[Step 7/8] Computing neighbor graph (n_neighbors={args.n_neighbors}, n_pcs={args.n_pcs})...")
    sc.pp.neighbors(
        adata_work,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        random_state=args.random_state,
    )
    print("Neighbor graph constructed successfully.")

    # Step 8: UMAP
    print(f"\n[Step 8/8] Computing new Scanpy UMAP (random_state={args.random_state})...")
    sc.tl.umap(adata_work, random_state=args.random_state)
    umap_scanpy = adata_work.obsm["X_umap"]
    umap_has_nan = bool(np.isnan(umap_scanpy).any())
    print(f"Scanpy UMAP computed with shape: {umap_scanpy.shape}, Has NaN: {umap_has_nan}")

    if umap_has_nan:
        raise ValueError("Scanpy UMAP coordinates contain NaN values!")

    adata.obsm["X_umap_scanpy"] = umap_scanpy.copy()

    # Step 9: Save lightweight outputs
    print("\n[Step 9/9] Saving lightweight outputs and validation summaries...")
    
    # Save UMAP coordinates + metadata table (compressed CSV)
    coord_df = pd.DataFrame(
        {
            "cell_id": adata.obs.index,
            "UMAP1_scanpy": umap_scanpy[:, 0],
            "UMAP2_scanpy": umap_scanpy[:, 1],
        },
        index=adata.obs.index,
    )
    if "X_umap_existing" in adata.obsm:
        coord_df["UMAP1_existing"] = adata.obsm["X_umap_existing"][:, 0]
        coord_df["UMAP2_existing"] = adata.obsm["X_umap_existing"][:, 1]

    # Include key metadata columns for downstream plotting without reloading heavy h5ad
    key_meta_cols = [
        "celltype_2",
        "genotype",
        "seurat_clusters",
        "RNA_snn_res.0.5",
        "integrated_snn_res.0.5",
        "sub.cluster",
        "celltype_individual",
        "celltype_2_old",
        "orig.ident",
        "time",
        "time_point",
        "Phase",
        "nCount_RNA",
        "nFeature_RNA",
        "percent.mt",
        "gene",
    ]
    for col in key_meta_cols:
        if col in adata.obs.columns:
            coord_df[col] = adata.obs[col].values

    coord_file_path = output_dir / "scanpy_umap_coordinates.csv.gz"
    coord_df.to_csv(coord_file_path, compression="gzip", index=False)
    print(f"Saved UMAP coordinates and metadata to: {coord_file_path}")

    # Save preprocessing summary JSON
    preprocessing_summary = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_top_genes": int(args.n_top_genes),
        "hvg_flavor": "seurat_v3",
        "hvg_batch_key": "orig.ident",
        "pca_n_comps": int(args.n_pcs),
        "pca_variance_explained_ratio_total": round(total_var_explained, 6),
        "pca_variance_ratio_top10": [round(float(x), 6) for x in variance_ratio[:10]],
        "neighbors_n_neighbors": int(args.n_neighbors),
        "neighbors_n_pcs": int(args.n_pcs),
        "umap_shape": list(umap_scanpy.shape),
        "umap_has_nan": umap_has_nan,
        "umap_scanpy_min": [round(float(umap_scanpy[:, 0].min()), 4), round(float(umap_scanpy[:, 1].min()), 4)],
        "umap_scanpy_max": [round(float(umap_scanpy[:, 0].max()), 4), round(float(umap_scanpy[:, 1].max()), 4)],
        "existing_umap_preserved": "X_umap_existing" in adata.obsm,
        "random_state": int(args.random_state),
    }

    summary_json_path = output_dir / "scanpy_preprocessing_summary.json"
    with open(summary_json_path, "w") as f:
        json.dump(preprocessing_summary, f, indent=2)
    print(f"Saved preprocessing summary to: {summary_json_path}")

    print("\n" + "=" * 80)
    print("PREPROCESSING COMPLETED SUCCESSFULLY!")
    print(f"  - Exactly {n_hvg_selected:,} HVGs selected.")
    print(f"  - PCA (50 PCs) completed, total variance explained: {total_var_explained * 100:.2f}%.")
    print(f"  - kNN graph (k=15) constructed.")
    print(f"  - Scanpy UMAP stored as obsm['X_umap_scanpy'] (0 NaNs).")
    print(f"  - Existing UMAP preserved as obsm['X_umap_existing'].")
    print("=" * 80)


if __name__ == "__main__":
    main()
