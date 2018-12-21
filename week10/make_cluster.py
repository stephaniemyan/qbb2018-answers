#!/usr/bin/env python3

"""
Usage: make_cluster.py <data.txt>

<data.txt> - A table of genes with corresponding FPKMs at different stages of cell differentiation

This code takes a dataset of gene expression in differentiating cells, and gives information about
the differentiation sequence and gene expression in different stages.

PART 1
Use Scipy to perform linkage analysis on gene expression, to get the order of differentiation.

PART 2
Use results from Part 1 to plot differentiation sequence on a dendrogram.

PART 3
Plot a heatmap of gene expression during all stages.

PART 4
Use k-means clustering to identify clusters of similarly expressed genes in two stages of your choice.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from sklearn.cluster import KMeans

f = open(sys.argv[1])

"""
PART 1

Do linkage analysis on gene expression data.
"""

df = pd.read_table(f, header=0, index_col=0)

# Pull out cell types and gene names from dataframe labels
cell_types = df.columns.values
gene_names = df.index.values

# Convert dataframe values into an array for linkage analysis
array = df.values

# Do linkage analysis on array (for gene order) with Ward's method
Z = linkage(array, "ward")
# leaves_list gives the order of items from left to right in the linkage dendrogram
gene_order = leaves_list(Z)
# Reorder gene names according to the leaves_list order, to label the dendrogram accurately
gene_names_ordered = gene_names[gene_order]

# Do linkage analysis on transposed array (for cell type order)
ZT = linkage(array.T, "ward")
# leaves_list gives the order of items from left to right in the linkage dendrogram
cell_order = leaves_list(ZT)
# Reorder cell types according to the leaves_list order, to label the dendrogram accurately
cell_types_ordered = cell_types[cell_order]


"""
PART 2

Plot a dendrogram displaying the order of cell type differentiation.
"""

fig, ax = plt.subplots()
# Use linkage data from the ZT array to get the order by cell type, rather than gene
dendrogram(
    ZT,
    show_leaf_counts=False,
    leaf_font_size=10,
    show_contracted=True,
    labels = cell_types_ordered)
ax.set_xlabel("Cell Types")
ax.set_ylabel("Distance")
ax.set_title("Dendrogram of cell type differentiation")
fig.savefig("dendrogram.png")
plt.close(fig)


"""
PART 3

Plot a heatmap of gene expression.
"""

# Reorganize (i.e. cluster) the array of FPKM values according to the gene order from linkage analysis
data = array[gene_order, :]
# Reorganize the reorganized array so that it also corresponds to the cell type order from linkage analysis
data = data[:, cell_order]

# Normalize gene expression data to keep heatmap colors within a range
X = (data-np.average(data,axis=0))/np.std(data,axis=0)
# Set max gene expression value to be the max of the normalized data
m = np.max(np.abs(X))

fig, ax = plt.subplots(figsize=(8, 6))
# Use Peter's code to create a heatmap and colorbar
im = ax.pcolor(X, cmap="RdBu", vmin=-1*m, vmax=m)
cbar = fig.colorbar(im, ax=ax)
ax.grid(False)
ax.set_xticks(
	np.arange(0.5, X.shape[1]+0.5),
	)
ax.set_xticklabels(cell_types_ordered, rotation=50)
ax.set_yticks([])
ax.set_title("Heatmap of gene expression by cell type")
fig.savefig("heatmap.png")
plt.close(fig)


"""
PART 4

Predict the # of expression clusters by eye from the heatmap.

Plot expression data for CFU and poly cells, and use kmeans clustering to predict clusters of similar
gene expression.
"""

# Pull out only CFU and Poly gene expression columns from the dataset
df2 = df.loc[:, ("CFU", "poly")]

# Convert dataframe values into an array for linkage analysis
array = df2.values

# Use KMeans function from sklearn to cluster data
kmeans = KMeans(n_clusters=5)
# I don't know what this part does, sorry Peter
kmeans = kmeans.fit(array)
# Use KMeans to predict clusters
labels = kmeans.predict(array)
# Get coordinates for centers of clusters
centroids = kmeans.cluster_centers_

# Plot clusters and centers
fig, ax = plt.subplots(figsize=(8, 6))
# Plot CFU and Poly FPKM values on x and y axes, and color dots by KMeans clustering
ax.scatter(array[:,0], array[:,1], c=labels, s=20)
# Plot predicted cluster centers
ax.scatter(centroids[:, 0], centroids[:, 1], c="black", s=200, alpha = 0.5)
ax.set_title("K-means clustering of CFU and poly gene expression (k=5)")
fig.savefig("kmeans.png")
plt.close(fig)