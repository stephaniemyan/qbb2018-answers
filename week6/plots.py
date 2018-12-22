#!/usr/bin/env python3

"""
Usage: ./plots.py <gained> <lost> <G1E file> <ER4 file> <features>

<gained> - Table of CTCF binding sites (start and end positions) gained from G1E to ER4 
differentiation
<lost> - Table of CTCF binding sites (start and end positions) lost from G1E to ER4 differentiation
<G1E file> - Table of CTCF binding sites in G1E cells
<ER4 file> - Table of CTCF binding sites in ER4 cells
<features> - Table of genome features (exon, intron, promoter) and their start/end positions

This script analyzes ChIP-seq peaks identified through MACS.

It outputs two plots:
1. Barplot showing # of CTCF binding sites gained and lost between states
2. Number of CTCF binding sites associated with each feature type
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

"""
PART 1

Count number of CTCF sites lost and gained between G1E and ER4.
"""

gained = open(sys.argv[1])
lost = open(sys.argv[2])

# Count number of lines in tje gained and lost files
gain_count = 0
for line in gained:
    gain_count +=1
loss_count = 0
for line in lost:
    loss_count += 1

# Get y and x values for barplot
y_vals = [gain_count, loss_count]
x_vals_1 = np.arange(2)


"""
PART 2

Count how many CTCF sites in G1E and ER4 overlap annotated genome features.

This method is unable to handle the following cases:

- Overlapping features. If an exon and intron overlap, the second annotated feature will overwrite
the part of the first that it overlaps.
- CTCF sites that only partially overlap features. If the CTCF site starts outside of a feature, it
will be counted as "other."
- CTCF sites that overlap multiple features. Every site will be counted in the feature that has its
start position.
"""

G1E_file = open(sys.argv[3])
ER4_file = open(sys.argv[4])
features_file = open(sys.argv[5])

# Make dictionary of all nucleotide positions that features cover, and which features they correspond to
positions = {}
for line in features_file:
    fields = line.rstrip("\r\n").split()
    # Get start and end position of each feature
    start = int(fields[1])
    end = int(fields[2])
    feat_type = fields[3]
    # For every nucleotide within the feature, annotate it with the feature type
    for position in range(start, end):
        positions[position] = feat_type

# Iterate through the list of CTCF sites and keep track of which features they overlap with
def count_features(filename):
    # Create dictionary where the keys are the feature types and the values are a running count of CTCF sites
    features = {"other":0}
    for line in filename:
        fields = line.rstrip("\r\n").split()
        start = int(fields[1])
        # With this method, every CTCF site will be categorized according to whether its start site is within a feature
        if start in positions:
            feat_type = positions[start]
            if feat_type in features:
                features[feat_type] += 1
            else:
                features[feat_type] = 1
        # If the start site doesn't fall within an annotated feature, count it as "other"
        else:
            features["other"] += 1
    return(features)

# Count features 
G1E_features = count_features(G1E_file)
ER4_features = count_features(ER4_file)

# Get x and y values for barplot
y_vals_G1E = [G1E_features["intron"], G1E_features["exon"], G1E_features["promoter"], \
            G1E_features["other"]]
y_vals_ER4 = [ER4_features["intron"], ER4_features["exon"], ER4_features["promoter"], \
            ER4_features["other"]]
x_vals_2 = np.arange(4)


"""
PART 3

Make plots.
"""

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(13,5))

# Plot barplot of gained/lost CTCF sites
ax1.bar(x_vals_1, y_vals, width=0.5, color="royalblue")
ax1.set_xticks(x_vals_1)
ax1.set_xticklabels(["Gained", "Lost"])
ax1.set_xlabel("Change during differentiation")
ax1.set_ylabel("Number of CTCF sites")
ax1.set_title("Change in CTCF binding sites during G1E to ER4 differentiation")

# Plot stacked barplot of CTCF sites overlapping with features
p1 = ax2.bar(x_vals_2, y_vals_G1E, color="rebeccapurple")
p2 = ax2.bar(x_vals_2, y_vals_ER4, bottom=y_vals_G1E, color="pink")
ax2.set_xticks(x_vals_2)
ax2.set_xticklabels(["Introns", "Exons", "Promoters", "Other"])
ax2.legend((p1[0], p2[0]), ("G1E", "ER4"))
ax2.set_xlabel("Type of region")
ax2.set_ylabel("Number of CTCF sites")
ax2.set_title("Types of features bound by CTCF")

plt.tight_layout()
fig.savefig("plots.png")
plt.close(fig)