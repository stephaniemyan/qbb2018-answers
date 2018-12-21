#!/usr/bin/env python3.6

"""
Usage: annotate.py

Include filepath to dataset in the "adata =" field. The dataset should be a binary file containing a 
cell x gene matrix, where each cell is matched to all of its aligned reads.

The program uses Scanpy's tools to:
1. Filter data (with the method used by Zheng et al. 2017)
2. Make PCA plots of the cells
3. Identify clusters with Louvain
4. Create UMAP and t-SNE plots of clusters
5. Identify genes distinguishing each cluster
6. Plot differential expression of genes of interest on t-SNE
7. Annotate clusters by cell type and re-plot t-SNE
"""

import sys
import matplotlib
matplotlib.use("Agg")

import scanpy.api as sc
sc.settings.autoshow = False

# Read 10x dataset
adata = sc.read_10x_h5("neuron_10k_v3_filtered_feature_bc_matrix.h5")
# Make variable names (in this case the genes) unique
adata.var_names_make_unique()


# Pre-filtering PCA
sc.tl.pca(adata, svd_solver='auto')
sc.pl.pca(adata, save="pre-filter.png", title="Pre-filtering PCA plot of scRNA-seq")


# Filter data using the Zheng 2017 recipe
sc.pp.recipe_zheng17(adata, n_top_genes=1000, log=True, plot=False, copy=False)


# Post-filtering PCA
sc.tl.pca(adata, svd_solver='auto')
sc.pl.pca(adata, save="post-filter.png", title="Post-filtering PCA plot of scRNA-seq")


# Use Louvain method (in Scanpy) to identify clusters in the data
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=None)
sc.tl.louvain(adata)

# Create UMAP plot of clusters
sc.tl.umap(adata)
sc.pl.umap(adata, color="louvain", save=".png", title="UMAP plot of neuron cell type clusters", legend_loc="right margin")

# Create t-SNE plot of clusters
sc.tl.tsne(adata)
sc.pl.tsne(adata, color="louvain", save=".png", title="t-SNE plot of neuron cell type clusters", legend_loc="right margin")


# Identify the genes that distinguish each cluster identified by the Louvain method
# Logistic regression approach
sc.tl.rank_genes_groups(adata, groupby="louvain", method="logreg")
sc.pl.rank_genes_groups(adata, save="_logreg.png")
# t-test approach
sc.tl.rank_genes_groups(adata, groupby="louvain", method="t-test")
sc.pl.rank_genes_groups(adata, save="_t-test.png")

# Plot dotplot showing expression of genes in Louvain clusters
sc.pl.rank_genes_groups_dotplot(adata, groupby="louvain", save=".png")


# Identify types of neurons by looking at expression of marker genes
# Peter showed me this way to figure out of my gene of interest is in the dataset:
genes = adata.var.gene_ids.index.get_values()
candidates = ['Sox17', 'Gad1', 'Gad2', 'Gad2', 'Lhx6', 'Hcls1', 'Hes1', 'Pax6', 'Chat']
present = [x for x in candidates if x in genes]
print(present)

# Plot t-SNE showing differential expression of genes of interest
sc.pl.tsne(adata, color=["Hcls1", "Sox17", "Gad1", "Htr3a", "Lhx6", "Sst", "Pax6", "Hes5"], \
            save="_markers.png")

# Plot t-SNE labeling groups by cell type
# Change louvain cluster labels to names of the cell types I can identify
cluster_names = ['0', 'SST interneurons', '2', 'Glial cells', '4', '5', '6', 'Radial glia 1', '8', '9', \
                'Radial glia 2', '11', '12', '13', 'VIP interneurons', '15', '16', \
                'PV interneurons', '18', 'Oligodendrocytes', '20', '21', '22']
adata.rename_categories('louvain', cluster_names)
# Re-plot t-SNE
sc.pl.tsne(adata, color="louvain", save="_annotated.png", title="Labeled t-SNE plot of neuron cell type clusters", \
            legend_loc="on data", legend_fontsize=6)
