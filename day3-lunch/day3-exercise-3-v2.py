#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    f = open(sys.argv[1])
else:
    f = sys.stdin
       
position = 21378950
shortest_dist_coding = 999999999999
shortest_dist_ncode = 999999999999
this_dist = 0
shortest_info = ""

# this finds all genes in 3R
for line in f:
    fields = line.strip().split("\t")
    if len(fields) < 2:
        continue 
    elif fields[0] == "3R" and fields[2] == "gene":
        gene_start = int(fields[3])
        gene_end = int(fields[4])
        if position < gene_start:
            this_dist = gene_start - position
        elif position > gene_end:
            this_dist = position - gene_end
        if this_dist < shortest_dist_coding and "protein_coding" in fields[-1]:
            shortest_dist_coding = this_dist
            shortest_info_coding = fields[-1]
        elif this_dist < shortest_dist_ncode:
            shortest_dist_ncode = this_dist
            shortest_info_ncode = fields[-1]
        
coding_split = shortest_info_coding.strip().split("; ")
ncode_split = shortest_info_ncode.strip().split("; ")

for i, item in enumerate(coding_split):
    if "gene_name" in item:
        coding_name = (item.strip().split())[1]
    i += 1
    
for i, item in enumerate(ncode_split):
    if "gene_name" in item:
        ncode_name = (item.strip().split())[1]
    i += 1
print("Nearest protein-coding gene: " + coding_name + ", " + str(shortest_dist_coding) + " bases away")
print("Nearest non-protein coding gene: " + ncode_name + ", " + str(shortest_dist_ncode) + " bases away")


