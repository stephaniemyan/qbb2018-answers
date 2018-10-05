#!/usr/bin/env python3

"""
Usage: ./make_pca.py <eigenvec file>

<eigenvec file> Eigenvector output file from PLINK

Makes a PCA plot from the PLINK analysis output.
"""

import sys
import matplotlib.pyplot as plt
plt.style.use('ggplot')

# Open eigenvector file
f = open(sys.argv[1])

# Parse file to extract PCA1 and PCA2 values
pca1_list = []
pca2_list = []
for line in f:
    fields = line.rstrip("\r\n").split()
    pca1 = float(fields[2])
    pca2 = float(fields[3])
    pca1_list.append(pca1)
    pca2_list.append(pca2)
    
# Plot in a scatterplot
fig, ax = plt.subplots()
ax.scatter(pca1_list, pca2_list, alpha=0.3, s=10, color="mediumvioletred")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("PCA")
plt.tight_layout()
fig.savefig("pca.png")
plt.close(fig)