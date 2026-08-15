import scanpy as sc

# Load your H5AD file
adata = sc.read_h5ad("diabetes.h5ad")

# Basic summary
print(adata)

# Dimensions
print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)

# Metadata columns
print("\nCell metadata columns (.obs):")
print(adata.obs.columns.tolist())

print("\nGene metadata columns (.var):")
print(adata.var.columns.tolist())

# Preview metadata
print("\nCell metadata:")
print(adata.obs.head())

print("\nGene metadata:")
print(adata.var.head())

# Additional stored information
print("\nUnstructured metadata (.uns):", list(adata.uns.keys()))
print("Embeddings (.obsm):", list(adata.obsm.keys()))
print("Layers:", list(adata.layers.keys()))
