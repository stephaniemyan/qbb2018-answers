#!/usr/bin/env python3

# day3-exercise-3-v2.py <ref genome> <chromosome> <position>
# Command line input 1: Filepath to reference genome
# Command line input 2: Chromosome name where your gene is
# Command line input 3: Position of your gene on chromosome

import sys

# Process command line inputs
f = open(sys.argv[1])
chrom = sys.argv[2]
position = int(sys.argv[3])

# To compare this_dist to; you want to start out really big
shortest_dist_coding = 999999999999
shortest_dist_ncode = 999999999999

# Set your looping this_dist variable to 0
this_dist = 0

# Get rid of headers
for line in f:
    fields = line.strip().split("\t")
    if len(fields) < 2:
        continue
        
# Grab lines that are on the chromosome and are genes
# Define gene start and end positions
    elif fields[0] == chrom and fields[2] == "gene":
        gene_start = int(fields[3])
        gene_end = int(fields[4])
        
# Figure out distance to gene from position of interest
        if position < gene_start:
            this_dist = gene_start - position
        elif position > gene_end:
            this_dist = position - gene_end
            
# For protein-coding genes:
# Compare this_dist to your shortest distance so far
# If it's shortest, output last column with info as shortest_info_coding
        if this_dist < shortest_dist_coding and "protein_coding" in fields[-1]:
            shortest_dist_coding = this_dist
            shortest_info_coding = fields[-1]
            
# For all other genes:
# Compare this_dist to your shortest distance so far
# If it's shortest, output last column with info as shortest_info_ncode
        elif this_dist < shortest_dist_ncode:
            shortest_dist_ncode = this_dist
            shortest_info_ncode = fields[-1]

# For protein-coding genes:
# Split items up by "; " to get the gene_name part
# Print out gene name and corresponding length
coding_split = shortest_info_coding.strip().split("; ")
for i, item in enumerate(coding_split):
    if "gene_name" in item:
        coding_name = (item.strip().split())[1]
    i += 1
print("Nearest protein-coding gene: " + coding_name + ", " + str(shortest_dist_coding) + " bases away")

# For all other genes:
# Split items up by "; " to get the gene_name part
# Print out gene name and corresponding length
ncode_split = shortest_info_ncode.strip().split("; ")
for i, item in enumerate(ncode_split):
    if "gene_name" in item:
        ncode_name = (item.strip().split())[1]
    i += 1
print("Nearest non-protein coding gene: " + ncode_name + ", " + str(shortest_dist_ncode) + " bases away")






