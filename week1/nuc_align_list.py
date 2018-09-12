#!/usr/bin/env python3

"""
Usage: ./nuc_align.py <BLAST output> <MAFFT output>

<BLAST output> Must be formatted as a FASTA file. Provides nucleotide alignment
<MAFFT output> Provides amino acid alignment

PART 1
For each aligned sequence, this script considers both the AA alignment and the original DNA 
sequence from BLAST. Wherever there is a gap in the AA alignment, it inserts three gaps in 
the DNA sequence.

This creates two lists of lists. The outer list contains each individual AA or DNA alignment.
Each individual alignment is a list of codons or AAs.

PART 2
The script counts the number of synonymous and nonsynonymous mutations at each codon position. 
Synonymous - codon changes, AA stays the same; nonsynonymous - codon changes, AA also changes.

This creates two lists. One lists the numbers of synonymous changes at each codon position.
The other lists the numbers of nonsynonymous changes.

PART 3
Run z-test

PART 4
Plot dN/dS vs. codon position. Highlight sites under positive selection in red.
"""

import sys
import fasta
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from statsmodels.stats import weightstats as stests

# Use FASTAReader to read BLAST and MAFFT output files
blast = fasta.FASTAReader(open(sys.argv[1]))
mafft = fasta.FASTAReader(open(sys.argv[2]))



# PART 1
# For every MAFFT AA alignment and its corresponding nucleotide alignment:
# Wherever there is a gap in the AA alignment, insert 3 nucleotide gaps (dashes ---) to the nucleotide alignment

# Create lists to add gapped DNA and AA alignments to
all_nuc_aligns = []
all_aa_aligns = []

# zip iterates through the BLAST and MAFFT files simultaneously
# At a given time, you are working with one specific AA alignment and its corresponding nucleotide alignment
for (dna_id, dna), (aa_id, aa) in zip(blast, mafft):
    
    # Create lists to add gaps/aligned nucleotides and AAs to
    new_nucs = []
    new_aas = []
    
    # j iterates separately through the nucleotide alignment. Iterating parallel with the AA alignment will skip over codons where the AA is a gap
    # j is outside the for loop so that it resets at 0 for each new alignment
    j = 0
    
    for i in range(0, len(aa)):
        
        # Define a as the current amino acid and nuc as the current codon
        a = aa[i]
        nuc = dna[j*3:(j+1)*3]
        # Add current amino acid to list of AAs
        new_aas.append(a)
        
        # If AA is a gap, add gap to nucleotide alignment
        if a == "-":
            new_nucs.append("---")
        # Otherwise, add corresponding codon to nucleotide alignment
        # Increment j to move on to next codon in the next loop
        else:
            new_nucs.append(nuc)
            j += 1
    
    # Add newly gapped AA and DNA alignments to their master lists
    all_nuc_aligns.append(new_nucs)
    all_aa_aligns.append(new_aas)



# PART 2
# Calculate the number of synonymous and nonsynonymous mutations at each position in the alignment

# Pull out the query AA and DNA sequences from their master lists
# (They should be the first item of each list)
query_aas = all_aa_aligns[0]
query_codons = all_nuc_aligns[0]

# Create lists to count the synonymous and nonsynonymous mutations in
list_of_dS = [0] * len(query_aas)
list_of_dN = [0] * len(query_aas)

# Count keeps track of how many positions there are across all the alignments. Is this even necessary? Probably not
count = 0
    
# Iterate through each alignment in the master list
for codons_list, aas_list in zip(all_nuc_aligns[1:], all_aa_aligns[1:]):
    
    # Iterate through every AA position in the alignment
    for i in range(0, len(codons_list)):
        
        # If the query codon is an indel, don't use it
        #if query_codons[i] = "---":
        #    continue
        # If the current codon does not match the corresponding query codon...
        elif codons_list[i] != query_codons[i]:
            # ...but the current AA matches the query AA, add to the count of synonymous mutations at this position
            if aas_list[i] == query_aas[i]:
                list_of_dS[i] += 1
            # ...and the current AA does not match the query AA, add to the count of nonsynonymous mutations at this position
            else:
                list_of_dN[i] += 1
        
        # Increment the other count (for all positions across all alignments)
        count += 1
        
# Here are a bunch of print statements for troubleshooting/checking things
# print(count)
# print(list_of_dS)
# print(list_of_dN)
print("Number of codons = " + str(len(query_codons)))
print("Number of AAs = " + str(len(query_aas)))
print("Number of nuc alignments = " + str(len(all_nuc_aligns)))
print("Number of aa alignments = " + str(len(all_aa_aligns)))
print("Length of dN list = " + str(len(list_of_dN)))
print("Length of dS list = " + str(len(list_of_dS)))
print("Number of nonsynonymous = " + str(sum(list_of_dN)))
print("Number of synonymous = " + str(sum(list_of_dS)))
print("Number of synonymous + nonsynonymous = " + str(sum(list_of_dS) + sum(list_of_dN)))



# PART 3
# Z-test. I don't even know

result = stests.ztest([list_of_dS[0], list_of_dN[0]])
#print(result)
pvalue = result[1]
#print(pvalue)



# PART 4
# Plot dN/dS ratios vs. codon position

# Create a list of dN/dS ratios by dividing the corresponding dN and dS numbers
list_of_ratios = []
for i in range(len(list_of_dS)):
    # I'm adding 1 to the denominator avoid dividing by 0. It's pretty sus but w/e
    ratio = (list_of_dN[i])/(list_of_dS[i] + 1)
    list_of_ratios.append(ratio)
    
# This is the version for subtracting dN - dS. I probably won't need it but idk what James Taylor wants from me right now
#list_of_difference = []
    #difference = (list_of_dN[i] - list_of_dS[i])
    #list_of_difference.append(difference)
#print(list_of_ratios)

x = range(0,len(list_of_ratios))
y = list_of_ratios

fig1, ax1 = plt.subplots()
ax1.scatter(x,y, alpha=0.3, s=2, color="blue")
ax1.set_xlabel("Codon position")
ax1.set_ylabel("dN/dS ratio")
# ax1.set_xscale('log')
# ax1.set_yscale('log')
ax1.set_title("dN vs dS")
fig1.savefig("dN_dS.png")
plt.close(fig1)




