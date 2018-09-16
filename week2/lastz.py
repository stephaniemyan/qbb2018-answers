#!/usr/bin/env python3

"""
Usage: ./lastz.py <LASTZ file>

<LASTZ file> A LASTZ general output file where the first three fields are the reference 
sequence's name, zstart, and zend.

Produces a dotplot of each contig's alignment to the reference genome. Every contig has a 
separate position on the x-axis (to account for duplicate coverage of one reference sequence 
site), and the script graphs position in contig sequence vs. position in reference genome.

For more aesthetic results, pre-sort the LASTZ file so that the contigs are ordered by start 
position in the reference sequence:
sort -k 2 -n <LASTZ file>

The default naming of the PNG works best for files named in this format: 
<lastz_[name]_sort.out>
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

# Open file and extract file name
f = open(sys.argv[1])
filename = (sys.argv[1].split("_"))[1]

# Plot dotplot: contig vs. positions on reference sequence
fig, ax = plt.subplots(figsize=(25, 10))

# Set counter for current position in x-axis (i.e. current contig)
x_pos = 0

for line in f:
    
    # Skip first line of file with column headers
    if line.startswith("#name1"):
        continue
    
    # Extract start position, end position, and length of contig in the reference sequence
    fields = line.rstrip("\r\n").split("\t")
    start = int(fields[1])
    end = int(fields[2])
    length = end - start
    
    # x range is the next spot on the x-axis, long enough for the whole contig
    x = np.linspace(x_pos,x_pos+length)  
    # y range is between this contig's start and end positions on the reference genome
    y = np.linspace(start,end)
    ax.plot(x, y)
    
    # Increment the x position counter
    x_pos += length
    
# I set the y-axis to end at 100000 for aesthetic purposes
# But this cuts out a handful of contigs at very far positions in the reference sequence
ax.set_ylim(0,100000)    

ax.set_xlabel("Contigs")
ax.set_ylabel("Reference sequence position")
ax.set_title("Contigs from " + filename + " aligned to reference sequence")
fig.savefig(filename + "_contigs.png")
plt.close(fig)