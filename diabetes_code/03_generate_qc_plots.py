#!/usr/bin/env python3
"""
03_generate_qc_plots.py

Generates publication-quality Quality Control (QC) plots for single-cell Perturb-seq analysis:
  1. Scatter plot: nFeature_RNA vs nCount_RNA (log-log scale, colored by percent.mt).
  2. Metric distributions: Histograms and density curves for nCount_RNA, nFeature_RNA, percent.mt.
  3. Sample-level QC: Violin/box plots of QC metrics grouped by orig.ident.
  4. Genotype-level QC: Violin/box plots of QC metrics grouped by genotype.
  5. Statistical Summary: Export qc_summary.csv with counts, min, max, median, mean, and std.

All plots are saved in both PNG (high-DPI) and vector PDF formats in results/figures/.

Usage:
  python 03_generate_qc_plots.py --coordinates results/scanpy_umap_coordinates.csv.gz --output-dir results
"""

import argparse
import os
import sys
from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description="Generate single-cell QC figures and summaries.")
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
        help="Directory to save figures and summary tables.",
    )
    return parser.parse_args()


def set_plotting_style():
    """Configure clean publication-quality styling."""
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["axes.edgecolor"] = "#333333"
    plt.rcParams["axes.linewidth"] = 1.0
    plt.rcParams["axes.grid"] = True
    plt.rcParams["grid.alpha"] = 0.3
    plt.rcParams["grid.linestyle"] = "--"


def plot_nfeature_vs_ncount(df, figures_dir):
    """Plot nFeature_RNA vs nCount_RNA with log-scale axes and percent.mt colormap."""
    fig, ax = plt.subplots(figsize=(8, 6.5))
    
    # Rasterized scatter for crisp vector output with 111k points
    scatter = ax.scatter(
        df["nCount_RNA"],
        df["nFeature_RNA"],
        c=df["percent.mt"],
        cmap="viridis",
        s=4,
        alpha=0.4,
        rasterized=True,
    )
    
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("Mitochondrial Reads (%)", fontsize=11, labelpad=8)
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Total UMI Counts per Cell (nCount_RNA)", fontsize=12, labelpad=8)
    ax.set_ylabel("Detected Genes per Cell (nFeature_RNA)", fontsize=12, labelpad=8)
    ax.set_title(f"Cell Sequencing Depth vs Gene Detection (n = {len(df):,} cells)", fontsize=13, pad=12, fontweight="bold")
    
    # Annotate medians
    med_count = df["nCount_RNA"].median()
    med_feat = df["nFeature_RNA"].median()
    ax.axvline(med_count, color="#d95f02", linestyle=":", linewidth=1.5, label=f"Median nCount: {med_count:,.0f}")
    ax.axhline(med_feat, color="#7570b3", linestyle=":", linewidth=1.5, label=f"Median nFeature: {med_feat:,.0f}")
    ax.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=10)
    
    plt.tight_layout()
    png_path = figures_dir / "qc_nfeature_vs_ncount.png"
    pdf_path = figures_dir / "qc_nfeature_vs_ncount.pdf"
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()
    print(f"Saved: {png_path} & {pdf_path}")


def plot_metric_distributions(df, figures_dir):
    """Plot distribution histograms and KDE for nCount_RNA, nFeature_RNA, and percent.mt."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    
    metrics = [
        ("nCount_RNA", "Total UMI Counts (nCount_RNA)", "#1f77b4", True),
        ("nFeature_RNA", "Detected Genes (nFeature_RNA)", "#2ca02c", False),
        ("percent.mt", "Mitochondrial Reads (%)", "#d62728", False),
    ]
    
    for ax, (col, label, color, log_scale) in zip(axes, metrics):
        data = df[col].dropna()
        if log_scale:
            sns.histplot(data, ax=ax, color=color, kde=True, log_scale=True, bins=50, alpha=0.6)
        else:
            sns.histplot(data, ax=ax, color=color, kde=True, bins=50, alpha=0.6)
        
        median_val = data.median()
        mean_val = data.mean()
        ax.axvline(median_val, color="#e41a1c", linestyle="--", linewidth=1.5, label=f"Median: {median_val:,.1f}")
        ax.axvline(mean_val, color="#377eb8", linestyle=":", linewidth=1.5, label=f"Mean: {mean_val:,.1f}")
        
        ax.set_title(label, fontsize=12, fontweight="bold", pad=10)
        ax.set_xlabel(label, fontsize=10)
        ax.set_ylabel("Number of Cells", fontsize=10)
        ax.legend(frameon=True, fontsize=9)
    
    plt.suptitle("Quality Control Metric Distributions", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    png_path = figures_dir / "qc_metric_distributions.png"
    pdf_path = figures_dir / "qc_metric_distributions.pdf"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {png_path} & {pdf_path}")


def plot_qc_by_sample(df, figures_dir):
    """Plot QC metrics stratified by sample / library origin (orig.ident)."""
    fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
    
    sample_order = df["orig.ident"].value_counts().index.tolist()
    
    # 1. nCount_RNA
    sns.violinplot(
        data=df, x="orig.ident", y="nCount_RNA", order=sample_order,
        ax=axes[0], palette="Blues_r", cut=0, inner="quartile", scale="width"
    )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("nCount_RNA (log scale)", fontsize=11)
    axes[0].set_title("Sequencing Depth by Sample (orig.ident)", fontsize=12, fontweight="bold")
    axes[0].grid(axis="y", linestyle="--", alpha=0.5)
    
    # 2. nFeature_RNA
    sns.violinplot(
        data=df, x="orig.ident", y="nFeature_RNA", order=sample_order,
        ax=axes[1], palette="Greens_r", cut=0, inner="quartile", scale="width"
    )
    axes[1].set_ylabel("nFeature_RNA", fontsize=11)
    axes[1].set_title("Gene Detection by Sample (orig.ident)", fontsize=12, fontweight="bold")
    axes[1].grid(axis="y", linestyle="--", alpha=0.5)
    
    # 3. percent.mt
    sns.violinplot(
        data=df, x="orig.ident", y="percent.mt", order=sample_order,
        ax=axes[2], palette="Reds_r", cut=0, inner="quartile", scale="width"
    )
    axes[2].set_ylabel("percent.mt (%)", fontsize=11)
    axes[2].set_title("Mitochondrial Read Percentage by Sample (orig.ident)", fontsize=12, fontweight="bold")
    axes[2].set_xlabel("Sample Identifier (orig.ident)", fontsize=11, labelpad=8)
    axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, ha="right", fontsize=10)
    axes[2].grid(axis="y", linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    png_path = figures_dir / "qc_by_sample.png"
    pdf_path = figures_dir / "qc_by_sample.pdf"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {png_path} & {pdf_path}")


def plot_qc_by_genotype(df, figures_dir):
    """Plot QC metrics stratified by Perturb-seq target genotype."""
    fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)
    
    genotype_order = df["genotype"].value_counts().index.tolist()
    
    # 1. nCount_RNA
    sns.boxplot(
        data=df, x="genotype", y="nCount_RNA", order=genotype_order,
        ax=axes[0], color="#9ecae1", showfliers=False
    )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("nCount_RNA (log scale)", fontsize=11)
    axes[0].set_title("Sequencing Depth across Perturb-seq Genotypes", fontsize=12, fontweight="bold")
    axes[0].grid(axis="y", linestyle="--", alpha=0.5)
    
    # 2. nFeature_RNA
    sns.boxplot(
        data=df, x="genotype", y="nFeature_RNA", order=genotype_order,
        ax=axes[1], color="#a1d99b", showfliers=False
    )
    axes[1].set_ylabel("nFeature_RNA", fontsize=11)
    axes[1].set_title("Gene Detection across Perturb-seq Genotypes", fontsize=12, fontweight="bold")
    axes[1].grid(axis="y", linestyle="--", alpha=0.5)
    
    # 3. percent.mt
    sns.boxplot(
        data=df, x="genotype", y="percent.mt", order=genotype_order,
        ax=axes[2], color="#fc9272", showfliers=False
    )
    axes[2].set_ylabel("percent.mt (%)", fontsize=11)
    axes[2].set_title("Mitochondrial Percentage across Perturb-seq Genotypes", fontsize=12, fontweight="bold")
    axes[2].set_xlabel("Genotype / Perturbation Target", fontsize=11, labelpad=8)
    axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=60, ha="right", fontsize=9)
    axes[2].grid(axis="y", linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    png_path = figures_dir / "qc_by_genotype.png"
    pdf_path = figures_dir / "qc_by_genotype.pdf"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {png_path} & {pdf_path}")


def generate_qc_summary_table(df, output_dir):
    """Compute summary statistics for QC metrics overall and per sample/genotype."""
    metrics = ["nCount_RNA", "nFeature_RNA", "percent.mt"]
    
    records = []
    
    # Overall summary
    for m in metrics:
        s = df[m].dropna()
        records.append({
            "group_type": "Overall",
            "group_name": "All_Cells",
            "metric": m,
            "count": int(len(s)),
            "min": float(s.min()),
            "max": float(s.max()),
            "median": float(s.median()),
            "mean": round(float(s.mean()), 4),
            "std": round(float(s.std()), 4),
            "q25": float(s.quantile(0.25)),
            "q75": float(s.quantile(0.75)),
        })
    
    # Per sample summary
    for sample, grp in df.groupby("orig.ident"):
        for m in metrics:
            s = grp[m].dropna()
            records.append({
                "group_type": "Sample (orig.ident)",
                "group_name": str(sample),
                "metric": m,
                "count": int(len(s)),
                "min": float(s.min()),
                "max": float(s.max()),
                "median": float(s.median()),
                "mean": round(float(s.mean()), 4),
                "std": round(float(s.std()), 4),
                "q25": float(s.quantile(0.25)),
                "q75": float(s.quantile(0.75)),
            })

    # Per genotype summary
    for geno, grp in df.groupby("genotype"):
        for m in metrics:
            s = grp[m].dropna()
            records.append({
                "group_type": "Genotype",
                "group_name": str(geno),
                "metric": m,
                "count": int(len(s)),
                "min": float(s.min()),
                "max": float(s.max()),
                "median": float(s.median()),
                "mean": round(float(s.mean()), 4),
                "std": round(float(s.std()), 4),
                "q25": float(s.quantile(0.25)),
                "q75": float(s.quantile(0.75)),
            })
    
    summary_df = pd.DataFrame(records)
    csv_path = output_dir / "qc_summary.csv"
    summary_df.to_csv(csv_path, index=False)
    print(f"Saved comprehensive QC summary table to: {csv_path}")
    return summary_df


def main():
    args = parse_args()
    coord_path = Path(args.coordinates).resolve()
    output_dir = Path(args.output_dir).resolve()
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("03_GENERATE_QC_PLOTS: QUALITY CONTROL VISUALIZATIONS & SUMMARY")
    print("=" * 80)
    print(f"Input file: {coord_path}")
    print(f"Figures directory: {figures_dir}")
    
    if not coord_path.exists():
        print(f"ERROR: Coordinate file {coord_path} does not exist.", file=sys.stderr)
        sys.exit(1)
        
    print("\n[Step 1/5] Loading coordinate and metadata table...")
    df = pd.read_csv(coord_path)
    print(f"Loaded metadata for {len(df):,} cells with columns: {list(df.columns)}")
    
    # Verify QC variables
    required_cols = ["nCount_RNA", "nFeature_RNA", "percent.mt", "orig.ident", "genotype"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Required QC column '{col}' missing from {coord_path}!")
        n_nan = df[col].isna().sum()
        if n_nan > 0:
            print(f"WARNING: Column '{col}' contains {n_nan} NaN values.")
    
    set_plotting_style()
    
    print("\n[Step 2/5] Generating nFeature vs nCount scatter plot...")
    plot_nfeature_vs_ncount(df, figures_dir)
    
    print("\n[Step 3/5] Generating QC metric distribution histograms...")
    plot_metric_distributions(df, figures_dir)
    
    print("\n[Step 4/5] Generating QC violin plots by sample (orig.ident)...")
    plot_qc_by_sample(df, figures_dir)
    
    print("\n[Step 5/5] Generating QC box plots by genotype...")
    plot_qc_by_genotype(df, figures_dir)
    
    print("\nGenerating statistical QC summary table...")
    generate_qc_summary_table(df, output_dir)
    
    print("\n" + "=" * 80)
    print("QC PLOTTING COMPLETED SUCCESSFULLY!")
    print("=" * 80)


if __name__ == "__main__":
    main()
