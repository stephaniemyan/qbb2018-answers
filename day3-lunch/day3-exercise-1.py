#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    f = open(sys.argv[1])
else:
    f = sys.stdin

# grab all lines that say column 3 is "gene"
# split up the last column by "; "
# see if gene_biotype and protein_coding are both in the same index item
# if so, add to count
# print count

count = 0
for line in f:
    fields = line.strip().split("\t")
    if len(fields) < 2:
        continue 
    elif fields[2] == "gene":
        fields_last_col = line.strip().split("; ")
        for i, item in enumerate(fields_last_col):
            if ("gene_biotype" in fields_last_col[i] and "protein_coding" in fields_last_col[i]):
                count += 1
            i += 1
print(count)