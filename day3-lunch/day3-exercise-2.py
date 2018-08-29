#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    f = open(sys.argv[1])
else:
    f = sys.stdin

# This generates a list with all the biotype names
# grab all lines that say column 3 is "gene"
# split up the last column by "; "
# see if gene_biotype is in the index item
# if so, take the gene_biotype item and make a list called 'last_col_split'
# last_col_split looks like ['gene_biotype', 'whatever', 'gene_biotype' etc]

biotypes = []
for line in f:
    fields = line.strip().split("\t")
    if len(fields) < 2:
        continue 
    elif fields[2] == "gene":
        fields_last_col = line.strip().split("; ")
        for i, item in enumerate(fields_last_col):
            if "gene_biotype" in fields_last_col[i]:
                last_col_split = fields_last_col[i].strip().split()
                for j, word in enumerate(last_col_split):
                    if word != "gene_biotype":
                        biotypes.append(word)
                    j += 1
            i += 1


# This counts unique items in the biotypes list, using a dictionary
# dictionary key is biotype name
# value is the count

dict = {}
for biotype_name in biotypes:
    if biotype_name in dict:
        dict[biotype_name] += 1
    else:
        dict[biotype_name] = 1
print(dict)