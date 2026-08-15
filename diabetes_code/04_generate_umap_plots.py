#!/usr/bin/env python3
"""
04_generate_umap_plots.py

Generates publication-quality single-cell UMAP visualizations using new Scanpy UMAP coordinates:
  - seurat_clusters
  - RNA_snn_res.0.5
  - integrated_snn_res.0.5
  - all existing resolution columns containing 'res.'
  - celltype_2
  - genotype
  - sub.cluster
  - celltype_individual
  - celltype_2_old
  - orig.ident
  - Phase

Saves high-resolution PNG (300 DPI) and vector PDF files in results/figures/.
Generates results/plot_manifest.csv summarizing all generated figures and category counts.

Usage:
  python 04_generate_umap_plots.py --coordinates results/scanpy_umap_coordinates.csv.gz --output-dir results
"""

import argparse
import os
import re
import sys
from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description="Generate publication-quality Scanpy UMAP figures.")
    parser.add_argument(
        "--coordinates",
        type=str,
        default="results/scanpy_umap_coordinates.csv.gz",
        help="Path to coordinates & metadata CSV file.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results",
        help="Directory to save figures and plot manifest.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=2.5,
        help="Scatter point size (default: 2.5).",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.6,
        help="Scatter point alpha transparency (default: 0.6).",
    )
    return parser.parse_args()


def set_plotting_style():
    """Configure clean publication-quality styling."""
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["axes.edgecolor"] = "#222222"
    plt.rcParams["axes.linewidth"] = 1.0


def get_color_palette(categories):
    """Generate a distinct, harmonious color palette for any number of categories."""
    n_cats = len(categories)
    if n_cats <= 10:
        base_colors = sns.color_palette("tab10", n_cats)
    elif n_cats <= 20:
        base_colors = sns.color_palette("tab20", n_cats)
    elif n_cats <= 28:
        base_colors = sns.color_palette("tab20b", 20) + sns.color_palette("tab20c", 20)
        base_colors = base_colors[:n_cats]
    else:
        # Use husl for high category counts to maintain distinction
        base_colors = sns.color_palette("husl", n_cats)
    
    color_map = {cat: base_colors[i] for i, cat in enumerate(categories)}
    return color_map


def sanitize_filename(name):
    """Convert column names like 'RNA_snn_res.0.5' to safe filenames 'RNA_snn_res_0_5'."""
    clean = re.sub(r"[^\w\-]", "_", name)
    return clean


def plot_single_umap(df, col, figures_dir, point_size=2.5, alpha=0.6):
    """Generate a single UMAP plot colored by categorical metadata column."""
    x = df["UMAP1_scanpy"].values
    y = df["UMAP2_scanpy"].values
    labels = df[col].astype(str).values
    
    # Sort categories by frequency or natural sort
    val_counts = df[col].value_counts(dropna=False)
    categories = list(val_counts.index.astype(str))
    
    # Try integer sort if all categories are integer digits
    if all(c.isdigit() for c in categories):
        categories = sorted(categories, key=lambda c: int(c))
    
    color_map = get_color_palette(categories)
    
    # Layout configuration
    n_cats = len(categories)
    legend_cols = 1 if n_cats <= 15 else (2 if n_cats <= 30 else 3)
    fig_width = 10.5 if legend_cols == 1 else (12.5 if legend_cols == 2 else 14.5)
    
    fig, ax = plt.subplots(figsize=(fig_width, 7.5))
    
    # Shuffle points so top classes don't completely occlude rare classes
    np.random.seed(42)
    shuffle_idx = np.random.permutation(len(df))
    
    cell_colors = [color_map[l] for l in labels[shuffle_idx]]
    
    scatter = ax.scatter(
        x[shuffle_idx],
        y[shuffle_idx],
        c=cell_colors,
        s=point_size,
        alpha=alpha,
        edgecolors="none",
        rasterized=True,
    )
    
    ax.set_xlabel("Scanpy UMAP 1", fontsize=12, labelpad=8, fontweight="medium")
    ax.set_ylabel("Scanpy UMAP 2", fontsize=12, labelpad=8, fontweight="medium")
    
    col_title = col.replace("_", " ").title() if not ("res." in col or "snn" in col) else col
    ax.set_title(f"Scanpy UMAP: {col} (n = {len(df):,} cells, {n_cats} groups)", fontsize=13, pad=12, fontweight="bold")
    
    # Clean spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle=":", alpha=0.4)
    
    # Create legend proxy artists
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label=f"{cat} ({val_counts.get(cat, val_counts.get(int(cat) if cat.isdigit() else cat, 0)):,})",
               markerfacecolor=color_map[cat], markersize=7)
        for cat in categories
    ]
    
    legend_fontsize = 9 if n_cats <= 20 else (8 if n_cats <= 35 else 7)
    ax.legend(
        handles=legend_elements,
        title=col,
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        frameon=True,
        framealpha=0.95,
        fontsize=legend_fontsize,
        title_fontsize=legend_fontsize + 1,
        ncol=legend_cols,
        borderaxespad=0.0,
    )
    
    plt.tight_layout()
    
    safe_name = sanitize_filename(col)
    png_name = f"umap_scanpy_{safe_name}.png"
    pdf_name = f"umap_scanpy_{safe_name}.pdf"
    
    png_path = figures_dir / png_name
    pdf_path = figures_dir / pdf_name
    
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    
    print(f"Generated UMAP for '{col}': {png_name} ({n_cats} categories)")
    
    return {
        "plot_file": png_name,
        "pdf_file": pdf_name,
        "metadata_column": col,
        "plot_type": "Scanpy_UMAP_Scatter",
        "n_cells": int(len(df)),
        "n_categories": int(n_cats),
        "n_missing": int(df[col].isna().sum()),
    }


def main():
    args = parse_args()
    coord_path = Path(args.coordinates).resolve()
    output_dir = Path(args.output_dir).resolve()
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("04_GENERATE_UMAP_PLOTS: SCANPY UMAP VISUALIZATION SUITE")
    print("=" * 80)
    print(f"Input coordinates: {coord_path}")
    print(f"Figures directory: {figures_dir}")
    
    if not coord_path.exists():
        print(f"ERROR: Coordinate file {coord_path} does not exist.", file=sys.stderr)
        sys.exit(1)
        
    print("\n[Step 1/3] Loading coordinate table and metadata...")
    df = pd.read_csv(coord_path)
    print(f"Loaded dataset: {len(df):,} cells with {len(df.columns)} columns.")
    
    # Check UMAP coordinates
    if "UMAP1_scanpy" not in df.columns or "UMAP2_scanpy" not in df.columns:
        raise ValueError("Missing 'UMAP1_scanpy' or 'UMAP2_scanpy' in coordinate file!")
        
    n_nan_coords = df["UMAP1_scanpy"].isna().sum() + df["UMAP2_scanpy"].isna().sum()
    if n_nan_coords > 0:
        raise ValueError(f"Found {n_nan_coords} NaN coordinates in Scanpy UMAP!")
    
    set_plotting_style()
    
    # Identify target columns to plot
    target_columns = [
        "seurat_clusters",
        "RNA_snn_res.0.5",
        "integrated_snn_res.0.5",
        "celltype_2",
        "genotype",
        "sub.cluster",
        "celltype_individual",
        "celltype_2_old",
        "orig.ident",
        "Phase",
    ]
    
    # Find any additional columns containing 'res.'
    for col in df.columns:
        if "res." in col and col not in target_columns:
            target_columns.append(col)
            
    print(f"\n[Step 2/3] Generating UMAP figures for {len(target_columns)} metadata variables:")
    manifest_records = []
    
    for col in target_columns:
        if col in df.columns:
            rec = plot_single_umap(
                df, col, figures_dir,
                point_size=args.point_size,
                alpha=args.alpha,
            )
            manifest_records.append(rec)
        else:
            print(f"Note: Column '{col}' not found in metadata. Skipping.")
            
    # Step 3: Export plot manifest
    print("\n[Step 3/3] Exporting plot manifest table...")
    manifest_df = pd.DataFrame(manifest_records)
    manifest_path = output_dir / "plot_manifest.csv"
    manifest_df.to_csv(manifest_path, index=False)
    print(f"Saved plot manifest to: {manifest_path}")
    
    print("\n" + "=" * 80)
    print(f"UMAP PLOTTING COMPLETED! Successfully created {len(manifest_records)} figures (PNG & PDF).")
    print("=" * 80)


if __name__ == "__main__":
    main()
