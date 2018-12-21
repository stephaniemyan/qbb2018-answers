#!/usr/bin/env python3

"""
Usage: diff_exp.py <data.txt>

<data.txt> - A table of genes with corresponding FPKMs at different stages

This finds genes that are differentially expressed at significant levels between two sets of stages.
It returns these genes of interest, as well as genes that are co-expressed with them.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans

f = open(sys.argv[1])
df = pd.read_table(f, header=0, index_col=0)


"""
PART 1

Use t-test to find genes that are differentially expressed at significant levels between two sets of
stages (early and late in differentiation).
"""

# Based on the dendrogram, I'm saying that unk and CFU are the two earliest stages, and mys and mid are the two latest stages
# If that's wrong I apologize in advance (sorry Peter)
# Pull out unk, CFU, mys, and mid gene expression columns; drop rows with NaN values
df = df.loc[:, ("unk", "CFU", "mys", "mid")].dropna()
early = ["unk", "CFU"]
late = ["mys", "mid"]
array = df.values

# Use Scipy t-test to calculate pvalues for different expression between early and late cell stages
stat, pvalue = ttest_ind(df[early], df[late], axis=1)
# Append pvalues to dataframe
df["pvalue"] = pvalue
# The Qi et al. paper used ratio change in addition to t-test to identify genes with significant differences in expression
# If I have time I'll do that too, but for now I'm just going to use the t-test

# Create new dataframe column with only pvalues < 0.05
df.loc[df['pvalue'] < 0.05, 'pvalue(significant)'] = df['pvalue']
# Make new dataframe with only genes that are significantly differentially expressed
df_sigs = df.dropna().loc[:,"pvalue(significant)"]
df_sigs.to_csv(sys.stdout, sep="\t", header=True)


"""
PART 2

Use kmeans clustering to find genes co-expressed with the differentially expressed genes that were
just identified.

This is a really weird method in which I assume that if a gene's kmeans "label" is the same as the 
"label" of my gene of interest, then they must be co-expressed. I have no idea if this is biologically
significant in any way.
"""

# Pull out all 500 genes from the original dataframe
gene_list = df.index.values
# Pull out differentially expressed genes identified in previous step
sig_gene_list = df_sigs.index.tolist()
# Get index values of all the diferentially expressed genes
sig_gene_indexes = []
for gene in sig_gene_list:
    index = df.index.get_loc(gene)
    sig_gene_indexes.append(index)
 
# Find kmeans cluster coordinates for all genes
kmeans = KMeans(n_clusters=10)
kmeans.fit(array)
labels = kmeans.predict(array)

# Iterate through the list of differentially expressed genes
# zip lets me get their index #s at the same time
for gene, index in zip(sig_gene_list, sig_gene_indexes):
    print("\n"+ "Genes co-expressed with " + gene + ":")
    # Get kmeans "label" for the gene of interest
    label_of_interest = labels[index]
    # If there are other genes with the same label as my gene of interest, print them out
    for i, label in enumerate(labels):
        if label == label_of_interest:
            print(gene_list[i])