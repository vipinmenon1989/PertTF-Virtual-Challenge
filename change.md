## 2026-08-15 — Count matrix and execution clarification

The AnnData object contains clear cell-level count summaries through:

- `nCount_RNA`
- `nFeature_RNA`

It also contains a matrix layer named `GPTin`.

The object summary alone does not identify whether raw counts are stored in `adata.X` or `adata.layers['GPTin']`. Antigravity must inspect both matrices, compare sampled row sums with `nCount_RNA`, and document the decision before Scanpy normalization or HVG selection.

The current QC, preprocessing, and plotting work must run inside an active Slurm allocation using the `scanpy_vipin` Conda environment.

Future environment changes must be documented in both `Project.md` and `change.md` before use.

## 2026-08-15 — Milestone 1: Slurm environment & AnnData count matrix validation

- **Date**: 2026-08-15
- **Environment**: Conda `scanpy_vipin` (Scanpy v1.11.5, AnnData v0.11.4, NumPy v2.2.6, SciPy v1.15.3, Matplotlib v3.10.8)
- **Slurm Allocation**:
  - `hostname`: `ihc-h200-1.igs.umaryland.edu`
  - `SLURM_JOB_ID`: `19687743`
  - `SLURM_JOB_NODELIST`: `ihc-h200-1`
- **Script Executed**:
  - `01_validate_anndata.py --input-h5ad diabetes.h5ad --output-dir results`
- **Key Validation Results**:
  1. **Dataset Dimensions**: 111,581 cells × 36,601 genes across 22 `.obs` columns, 50 `.var` columns, and 1 existing embedding (`obsm['X_umap']`).
  2. **Matrix Comparison (`adata.X` vs `adata.layers['GPTin']`)**:
     - Both matrices are identical sparse CSR `float64` structures with `nnz = 402,749,681` (90.1383% sparsity).
     - $\max |\text{adata.X} - \text{adata.layers['GPTin']}| = 0.0$.
     - Sampled values are non-negative (`True`), continuous / non-integer (`is_integer_valued = False`), ranging between 0.0 and ~7.896.
     - Sampled row sums have mean 3,853.67, std 452.48, and range [1544.73, 4933.94].
     - Sampled row sums correlate strongly with `adata.obs['nCount_RNA']` (Pearson $r = 0.7745$, Spearman $r = 0.8348$).
  3. **Count Matrix Interpretation & Mathematical Verification**:
     - Neither matrix represents raw integer UMI counts; both store Seurat `LogNormalize` expression ($\log(1 + 10000 \cdot \frac{\text{count}}{\text{nCount\_RNA}})$).
     - The true raw integer count matrix is reconstructible with machine precision ($\max |\text{diff from integer}| < 10^{-12}$, $\text{diff from nCount\_RNA} = 0.0$) using:
       $$\text{raw\_counts} = \text{round}\left(\frac{\text{expm1}(X) \cdot \text{nCount\_RNA}}{10000}\right)$$
  4. **Downstream Preprocessing Decision**:
     - `adata.X` must NOT be subjected to a second round of `sc.pp.normalize_total()` and `sc.pp.log1p()`.
     - For `flavor='seurat_v3'` HVG selection, raw counts will be reconstructed into a dedicated working layer.
- **Generated Outputs**:
  - `results/reports/validation_summary.json`
  - `results/reports/validation_report.txt`
- **Unresolved Issues**: None. Milestone 1 completed cleanly.

## 2026-08-15 — Milestone 2: Scanpy preprocessing pipeline & UMAP generation

- **Date**: 2026-08-15
- **Environment**: Conda `scanpy_vipin` (Scanpy v1.11.5, AnnData v0.11.4, NumPy v2.2.6, SciPy v1.15.3, Matplotlib v3.10.8)
- **Slurm Allocation**:
  - `hostname`: `ihc-h200-1.igs.umaryland.edu`
  - `SLURM_JOB_ID`: `19687743`
  - `SLURM_JOB_NODELIST`: `ihc-h200-1`
- **Script Executed**:
  - `02_scanpy_preprocessing.py --input-h5ad diabetes.h5ad --output-dir results`
- **Parameters**:
  - `n_top_genes = 5000`
  - `flavor = "seurat_v3"`
  - `layer = "counts"`
  - `batch_key = "orig.ident"`
  - `target_sum = 10000`
  - `n_pcs = 50`
  - `n_neighbors = 15`
  - `random_state = 42`
- **Validation & Key Findings**:
  1. **Full-Dataset Count Reconstruction**: 111,581 of 111,581 cells (100.0%) verified against `adata.obs['nCount_RNA']` (`max_diff = 0.0`).
  2. **HVG Selection**: Exactly 5,000 highly variable genes identified.
  3. **PCA**: 50 PCs computed without NaN coordinates; 43.09% cumulative variance explained.
  4. **kNN Graph**: 15 neighbors per cell constructed in 50-PC space.
  5. **UMAP Embedding**: New Scanpy UMAP coordinates calculated (`obsm['X_umap_scanpy']`, shape `(111581, 2)`, 0 NaNs).
  6. **Embedding Preservation**: Original Seurat UMAP preserved intact as `obsm['X_umap_existing']`.
  7. **No Large `.h5ad` Files Written**: Data was processed in memory, and lightweight summaries/coordinates were exported.
- **Generated Outputs**:
  - `results/count_reconstruction_validation.json`
  - `results/hvg_summary.csv`
  - `results/scanpy_preprocessing_summary.json`
  - `results/scanpy_umap_coordinates.csv.gz`
- **Acceptance Criteria Verification**:
  - Exactly 5,000 HVGs selected: **PASS**
  - Reconstructed row sums match `nCount_RNA`: **PASS** (100.0%)
  - PCA completed without missing values: **PASS** (0 NaNs)
  - Neighbors created: **PASS** ($k=15$)
  - `X_umap_scanpy` exists and valid: **PASS** (111,581 × 2, 0 NaNs)
  - No `.h5ad` file staged for Git: **PASS**
  - `python -m compileall .` clean: **PASS**
- **Unresolved Issues**: None. Milestone 2 completed cleanly.

## 2026-08-15 — Milestone 3: Publication-quality QC figures and statistics

- **Date**: 2026-08-15
- **Environment**: Conda `scanpy_vipin` (Scanpy v1.11.5, Matplotlib v3.10.8, Seaborn v0.13.2)
- **Slurm Allocation**:
  - `hostname`: `ihc-h200-1.igs.umaryland.edu`
  - `SLURM_JOB_ID`: `19687743`
  - `SLURM_JOB_NODELIST`: `ihc-h200-1`
- **Script Executed**:
  - `03_generate_qc_plots.py --coordinates results/scanpy_umap_coordinates.csv.gz --output-dir results`
- **Generated Figures & Files**:
  - `results/figures/qc_nfeature_vs_ncount.png` and `.pdf`
  - `results/figures/qc_metric_distributions.png` and `.pdf`
  - `results/figures/qc_by_sample.png` and `.pdf`
  - `results/figures/qc_by_genotype.png` and `.pdf`
  - `results/qc_summary.csv`
- **Validation**:
  - All 111,581 cells plotted with log-log scaling, transparent rasterized markers, and median indicators.
  - Distribution histograms and KDE density curves constructed for `nCount_RNA`, `nFeature_RNA`, and `percent.mt`.
  - Violin and box plots stratified by sample origin (`orig.ident`) and perturbation target (`genotype`).
  - No missing or infinite values encountered across QC metrics.

## 2026-08-15 — Milestone 4: Scanpy UMAP visualization suite & manifest

- **Date**: 2026-08-15
- **Environment**: Conda `scanpy_vipin` (Scanpy v1.11.5, Matplotlib v3.10.8, Seaborn v0.13.2)
- **Slurm Allocation**:
  - `hostname`: `ihc-h200-1.igs.umaryland.edu`
  - `SLURM_JOB_ID`: `19687743`
  - `SLURM_JOB_NODELIST`: `ihc-h200-1`
- **Script Executed**:
  - `04_generate_umap_plots.py --coordinates results/scanpy_umap_coordinates.csv.gz --output-dir results`
- **Generated Figures (PNG & Vector PDF)**:
  1. `results/figures/umap_scanpy_seurat_clusters.png` and `.pdf` (22 clusters)
  2. `results/figures/umap_scanpy_RNA_snn_res_0_5.png` and `.pdf` (22 clusters)
  3. `results/figures/umap_scanpy_integrated_snn_res_0_5.png` and `.pdf` (22 clusters)
  4. `results/figures/umap_scanpy_celltype_2.png` and `.pdf` (15 cell types)
  5. `results/figures/umap_scanpy_genotype.png` and `.pdf` (37 genotypes)
  6. `results/figures/umap_scanpy_sub_cluster.png` and `.pdf` (20 sub-clusters)
  7. `results/figures/umap_scanpy_celltype_individual.png` and `.pdf` (18 categories)
  8. `results/figures/umap_scanpy_celltype_2_old.png` and `.pdf` (14 categories)
  9. `results/figures/umap_scanpy_orig_ident.png` and `.pdf` (13 samples)
  10. `results/figures/umap_scanpy_Phase.png` and `.pdf` (3 cell cycle phases)
- **Manifest**:
  - Generated `results/plot_manifest.csv` tracking all 10 visual suites across 111,581 cells with 0 missing coordinates.
- **Validation**:
  - All plots verified non-empty, legible legends, distinct color palettes, and clean layout.
  - `python -m compileall` passed cleanly.
  - `git status` confirmed no `.h5ad` files staged or tracked.
