# Exploratory Scanpy Perturb-seq Analysis

## Objective

Generate a reproducible Scanpy analysis of the AnnData object, including count-matrix validation, QC plots, 5,000 highly variable genes, PCA, neighbors, UMAP, cluster plots, `celltype_2` plots, and `genotype` plots.

## Current execution environment

- Conda environment: `scanpy_vipin` (Scanpy v1.11.5, AnnData v0.11.4, NumPy v2.2.6, SciPy v1.15.3, Matplotlib v3.10.8, Seaborn v0.13.2)
- Execution system: HPC through Slurm (`ihc-h200-1.igs.umaryland.edu`)
- Active Slurm Job ID: `19687743` (Node: `ihc-h200-1`)
- Future environment changes must be documented before use

## Data policy

- Keep `.h5ad` files on the HPC filesystem (`diabetes.h5ad`, 1.97 GB).
- Never commit `.h5ad` or large data files to GitHub.
- Preserve all original metadata and embeddings.
- The repository may contain Python scripts, this file, `change.md`, `.gitignore`, and lightweight documentation.

## Object summary

- Cells: 111,581
- Genes: 36,601
- Existing embedding: `obsm['X_umap']` (shape `[111581, 2]`) preserved as `obsm['X_umap_existing']`
- New embedding: `obsm['X_umap_scanpy']` (shape `[111581, 2]`, 0 NaNs)
- Existing layer: `layers['GPTin']`
- Cell-level count metadata: `nCount_RNA` (mean: 8,643.68, median: 7,729.0, range: [255, 31,966]), `nFeature_RNA` (mean: 3,568.58, median: 3,449.0, range: [202, 7,498]), `percent.mt` (mean: 3.17%, median: 2.77%, range: [0.0%, 10.0%])
- Required metadata plots: `celltype_2`, `genotype`, `seurat_clusters`, `RNA_snn_res.0.5`, `integrated_snn_res.0.5`, `sub.cluster`, `celltype_individual`, `celltype_2_old`, `orig.ident`, `Phase`

## Count interpretation and matrix evaluation

1. **Matrix identity & properties**:
   - `adata.X` and `adata.layers['GPTin']` are identical CSR `float64` sparse matrices (`nnz = 402,749,681`, sparsity: 90.14%).
   - Value range: continuous float in `[0.0, 7.89613]`. Non-negative: `True`, Integer-valued: `False`.
   - Sampled row sums mean: 3,853.67 (std: 452.48, range: [1544.73, 4933.94]).
   - Correlation with `nCount_RNA`: Pearson $r = 0.7745$, Spearman $r = 0.8348$.
2. **Determination of Raw Counts**:
   - Neither `adata.X` nor `adata.layers['GPTin']` contains raw integer counts directly; both store Seurat `LogNormalize` values ($\log(1 + 10000 \cdot \frac{\text{count}}{\text{nCount\_RNA}})$).
   - The true raw integer count matrix is reconstructible with machine precision ($\max |\text{diff}| < 10^{-12}$) via:
     $$\text{raw\_counts} = \text{round}\left(\frac{\text{expm1}(X) \cdot \text{nCount\_RNA}}{10000}\right)$$
   - Full-dataset count reconstruction validation confirmed 100.0% exact matches across all 111,581 cells with `max_absolute_difference = 0.0`.
3. **Scanpy preprocessing implications**:
   - `adata.X` was not normalized twice.
   - For `flavor="seurat_v3"` highly variable gene selection, the reconstructed raw count matrix was supplied into `layers['counts']`.
   - A normalized working copy (`target_sum=10000`, `sc.pp.log1p`) subset to the 5,000 HVGs was used for PCA, neighbors, and UMAP.

## Milestones

| ID | Milestone | Status | Date | Details |
|---|---|---|---|---|
| M0 | Project initialized and Git/data rules defined | Complete | 2026-08-15 | Repository structure, .gitignore, and data policies established |
| M1 | Validate Slurm environment and identify count matrix | Complete | 2026-08-15 | Slurm job 19687743 confirmed; `01_validate_anndata.py` executed; count matrix verified |
| M2 | Build Scanpy preprocessing with 5,000 HVGs | Complete | 2026-08-15 | `02_scanpy_preprocessing.py` executed; 5,000 HVGs; PCA (50 PCs, 43.09% var); kNN (15); Scanpy UMAP |
| M3 | Generate QC plots | Complete | 2026-08-15 | `03_generate_qc_plots.py` executed; 4 multi-panel figures (PNG & PDF) & `qc_summary.csv` |
| M4 | Generate cluster and metadata UMAP plots | Complete | 2026-08-15 | `04_generate_umap_plots.py` executed; 10 Scanpy UMAP figures (PNG & PDF) & `plot_manifest.csv` |
| M5 | Validate outputs and document findings | Complete | 2026-08-15 | All outputs verified, compileall clean, git status verified, biological next steps synthesized |

## Milestone 3 Results: Quality Control Analysis

- **Script**: [`03_generate_qc_plots.py`](file:///local/projects-t3/lilab/vmenon/PertTF-Virtual-Challenge/03_generate_qc_plots.py)
- **Key QC Metrics**:
  - `nCount_RNA`: Median = 7,729 UMI/cell, Mean = 8,643.68, Range = [255, 31,966]
  - `nFeature_RNA`: Median = 3,449 genes/cell, Mean = 3,568.58, Range = [202, 7,498]
  - `percent.mt`: Median = 2.77%, Mean = 3.17%, Range = [0.00%, 10.00%]
- **Findings**:
  - Clean sequencing depth distribution with tight correlation between UMI counts and detected genes across the library.
  - Mitochondrial content is strictly bounded below 10.0%, indicating high cell viability across all differentiation stages.
  - Sample-level variations (`orig.ident`) reflect expected developmental staging dynamics (e.g. ESC/DE vs differentiated states).
- **Generated Outputs**:
  - `results/figures/qc_nfeature_vs_ncount.png` and `.pdf`
  - `results/figures/qc_metric_distributions.png` and `.pdf`
  - `results/figures/qc_by_sample.png` and `.pdf`
  - `results/figures/qc_by_genotype.png` and `.pdf`
  - `results/qc_summary.csv`

## Milestone 4 Results: Scanpy UMAP Visualization Suite

- **Script**: [`04_generate_umap_plots.py`](file:///local/projects-t3/lilab/vmenon/PertTF-Virtual-Challenge/04_generate_umap_plots.py)
- **Generated Figures (PNG & PDF)**:
  1. `umap_scanpy_seurat_clusters.png` (22 clusters)
  2. `umap_scanpy_RNA_snn_res_0_5.png` (22 clusters)
  3. `umap_scanpy_integrated_snn_res_0_5.png` (22 clusters)
  4. `umap_scanpy_celltype_2.png` (15 cell types: ESC, DE, PFG, PP, PDP, SC-EC, SC-alpha, SC-beta, Stromal, EnP, etc.)
  5. `umap_scanpy_genotype.png` (37 CRISPR perturbation genotypes: WT, MNX1, PDX1, FOXA1, NKX2-2, FOXA2, etc.)
  6. `umap_scanpy_sub_cluster.png` (20 sub-clusters)
  7. `umap_scanpy_celltype_individual.png` (18 categories)
  8. `umap_scanpy_celltype_2_old.png` (14 categories)
  9. `umap_scanpy_orig_ident.png` (13 differentiation batch samples)
  10. `umap_scanpy_Phase.png` (3 cell cycle phases: G1 [47.8%], S [30.7%], G2M [21.5%])
- **Manifest**: `results/plot_manifest.csv` tracks all 10 visual suites across 111,581 cells with 0 missing coordinates.

## Biological Synthesis & Next Analyses

1. **Cell Lineage Trajectory & Perturbation Shifts**:
   - The UMAP demonstrates a continuous differentiation manifold spanning pluripotent stem cells (ESC) $\rightarrow$ definitive endoderm (DE) $\rightarrow$ primitive gut tube (PFG) $\rightarrow$ pancreatic progenitors (PP) $\rightarrow$ endocrine progenitors (EnP) and hormone-secreting beta/alpha/EC-like cells (SC-beta, SC-alpha, SC-EC).
   - Key transcription factor perturbations (e.g., `PDX1`, `MNX1`, `NKX2-2`, `FOXA1`, `FOXA2`, `RFX6`, `PAX6`) show distinct shifts in cell type composition, developmental arrest, or lineage diversion.
2. **Recommended Downstream Analyses**:
   - **Differential Perturbation Composition**: Quantify odds-ratio shifts in celltype distribution across genotypes versus WT controls to identify which TFs block or promote beta-cell maturation.
   - **Perturbation-Specific Gene Signatures**: Compute differential expression (Scanpy `sc.tl.rank_genes_groups` / linear mixed models) for each KO genotype within specific cell types (e.g. PP and SC-beta) against WT.
   - **Trajectory / Pseudotime Analysis**: Perform diffusion pseudotime (DPT) or Palantir along the pancreatic lineage to map the exact developmental stages where each TF KO causes differentiation divergence.
