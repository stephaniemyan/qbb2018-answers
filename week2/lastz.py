#!/usr/bin/env python3

"""
Usage: ./lastz.py <velvet.fa>

sort -k 2 -n lastz_velvet.out
"""

import sys
import os
import fasta
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

f = open(sys.argv[1])
filename = (sys.argv[1].split("_"))[1]

fig, ax = plt.subplots(figsize=(25, 10))

x_pos = 0
for line in f:
    if line.startswith("#name1"):
        continue
    fields = line.rstrip("\r\n").split("\t")
    start = int(fields[1])
    end = int(fields[2])
    length = end - start
    y = np.linspace(start,end)
    x = np.linspace(x_pos,x_pos+length)
    ax.plot(x, y)
    x_pos += length
    
ax.set_xlabel("Contigs")
ax.set_ylabel("Reference sequence position")
ax.set_title("Contigs from " + filename + " aligned to reference sequence")
ax.set_ylim(0,100000)
fig.savefig(filename + "_contigs.png")
plt.close(fig)