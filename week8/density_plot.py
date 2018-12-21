#!/usr/bin/env python3

"""
Usage: density_plot.py <bedtools_file.txt>

<bedtools_file.txt> - A table containing the start and end positions of ChIP-seq peaks, and the 
start position of the motif associated with each peak

This script plots a histogram of the relative start positions of identified motifs within ChIP-seq
peaks.
"""

import sys
import matplotlib.pyplot as plt
plt.style.use('ggplot')

f = open(sys.argv[1])

positions = []

for line in f:
    fields = line.rstrip("\t").split()
    # Pull out ChIP-seq peak start and end position, and motif start position
    seq_start = int(fields[1])
    seq_end = int(fields[2])
    motif_start = int(fields[13])
    
    # Calculate length of the ChIP-seq peak
    seq_length = seq_end - seq_start
    # Calculate distance between motif start position and the ChIP-seq peak start position
    position = motif_start - seq_start
    # Calculate relative motif start position
    rel_pos = position/seq_length
    positions.append(rel_pos)
    
# Plot relative start positions in a histogram
fig, ax = plt.subplots(figsize=(8,5))
ax.hist(positions, color="royalblue", bins=30)
ax.set_title("Distribution of relative motif locations in ChIP-seq peaks")
ax.set_ylabel("Number of motifs")
ax.set_xlabel("Relative location within ChIP-seq peak")
plt.savefig("position_freq.png")
plt.close(fig)