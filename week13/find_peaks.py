#!/usr/bin/env python2.7

"""
Usage: find_peaks.py <CTCF file>

<CTCF file> - A table with the start and end positions of CTCF binding sites

This script takes hifive data and incorporates it with ChIP-seq data to find regions of the genome 
where your protein of interest binds and creates interactions between DNA.

In this case, it looks at a specific region of mouse chr17 to finds all of the hifive interaction 
enrichments greater than 1 that have at least one CTCF peak at both ends.
"""

import sys
import numpy as np
import hifive
import pandas as pd

# This is a file output from hifive with interactions between regions of a genome
hic = hifive.HiC('hifive.hcp')


"""
PART 1

Make a 2D enrichment matrix from the hifive file. All of this code was written by Mike Sauria.
"""

# Get data into numpy 3D array
data = hic.cis_heatmap(chrom='chr17', start=15000000, stop=17500000, binsize=10000, \
datatype='fend', arraytype='full')

# Make square enrichment matrix
data[:, :, 1] *= np.sum(data[:, :, 0]) / np.sum(data[:, :, 1])

# Finds bins where expected value is > 0 and only preserve those bins in the data matrix
where = np.where(data[:, :, 1] > 0)
data[where[0], where[1], 0] /= data[where[0], where[1], 1]
data = data[:, :, 0]


"""
PART 2

Pull out all CTCF binding sites and convert their locations into bin positions.
"""

bins = []
for line in open(sys.argv[1]):
    
    # Skip first line
    if "length" in line:
        continue
    # Pull out start and end positions
    fields = line.rstrip("\r\n").split()
    start = int(fields[1])
    end = int(fields[2])
    
    # Only look at CTCF occupancy sites within the range of the enrichment matrix
    if (fields[0] == "chr17") and (start >= 15000000) and (end <= 17500000):
        # Find midpoint of CTCF site
        midpoint = (end + start)/2
        # Calculate the bin that the midpoint is in. The enrichment matrix bins begin at position 15000000, and the bin size is 10000
        this_bin = (midpoint - 15000000)/10000
        bins.append(this_bin)

# Take out redundant bins
unique_bins = np.unique(bins)


"""
PART 3

Find CTCF interaction sites that have enrichments >= 1. Convert the bin positions of these sites back
into positions in the genome, and sort them by enrichment value.
"""

starts = []
ends = []
enrichments = []

# Iterate through the enrichment matrix, finding the enrichment values at every combination of bins in the bins list
# If the enrichment value is greater than 1, record the value and its start and end bins
# tbh I don't really know why I'm doing this. But it works
for i in range(len(unique_bins)):
    for j in range(i, len(unique_bins)):
        enrichment = data[unique_bins[i],unique_bins[j]]
        if enrichment >= 1:
            # Convert start and end bins to genomic positions
            start_pos = (unique_bins[i]*10000) + 15000000
            end_pos = (unique_bins[j]*10000) + 15000000
            starts.append(start_pos)
            ends.append(end_pos)
            enrichments.append(enrichment)

# Create dataframe with CTCF start and end positions and corresponding enrichment values
df = pd.DataFrame(
    {"start" : starts,
     "end" : ends,
     "enrichment_value" : enrichments
    })

# Reorder columns
col_order = ["start", "end", "enrichment_value"]
df = df.reindex(columns=col_order)

# Sort dataframe from highest to lowest enrichment value
df = df.sort_values(by="enrichment_value", ascending=False)
df.to_csv(sys.stdout, sep="\t", index=False)