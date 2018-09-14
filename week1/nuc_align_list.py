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
Calculate dN - dS for each position, and calculate the z value of that difference to 
determine if it's significantly different from the other dN - dS values. The null hypothesis
is dN = dS, or no selection; positions with significant z values are undergoing selection.

PART 4
Plot dN/dS vs. codon position. Highlight sites under positive selection in red.
"""

import sys
import fasta
from math import sqrt
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
# Create a list to count the # of changes, synonymous or nonsynonymous, at each position
count = [0] * len(query_aas)

# count_indels keeps track of all the --- codons in the query sequence. This is just for fun
count_indels = 0

# Iterate through each alignment in the master list
for codons_list, aas_list in zip(all_nuc_aligns[1:], all_aa_aligns[1:]):
    
    # Iterate through every AA position in the alignment
    for i in range(0, len(codons_list)):
        
        # If the query codon is an indel, don't use it, and increment the indels count
        if query_codons[i] == "---":
            count_indels += 1
        # If the current codon does not match the corresponding query codon...
        elif codons_list[i] != query_codons[i]:
            count[i] += 1
            # ...but the current AA matches the query AA, add to the count of synonymous mutations at this position
            if aas_list[i] == query_aas[i]:
                list_of_dS[i] += 1
            # ...and the current AA does not match the query AA, add to the count of nonsynonymous mutations at this position
            else:
                list_of_dN[i] += 1
        
# Here are a bunch of print statements for troubleshooting/checking things
print("Total indels in all alignments = " + str(count_indels))
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
# Z-test each dN - dS value to find out if it is significantly different from expectation
# Null hypothesis: dN - dS = 0, or no selection. Calculate z for each position, and determine if the z is significant at p < 0.001
# If a site is undergoing significant positive selection, we expect dN > dS (mutations that change the AA are encouraged), so dN - dS should be negative

# Create a list of differences for dN - dS
list_of_difference = []
for i in range(len(list_of_dS)):
    difference = (list_of_dN[i] - list_of_dS[i])
    list_of_difference.append(difference)

# Create a list of dN/dS ratios by dividing the corresponding dN and dS numbers
# This is what is actually being plotted
list_of_ratios = []
for i in range(len(list_of_dS)):
    # I'm adding 1 to the denominator avoid dividing by 0. It's pretty sus but w/e
    ratio = (list_of_dN[i])/(list_of_dS[i] + 1)
    list_of_ratios.append(ratio)

# Calculate z score for every dN - dS difference. These z scores will tell you how far from the mean (in std devs) a particular dN - dS difference is
# z = (x - mu) / sigma
# x: dN - dS value at that position
# mu: "Population mean" - for some reason in this case it's the null hypothesis, 0
# sigma: Standard error - stdev of all the dN - dS values, divided by the # of samples (i.e. # of changes at that position)

# For plotting purposes, make dictionaries to store the points that undergo significant positive selection and the points that don't
pos_selection = {}
other = {}

std = np.std(list_of_difference)
for i, difference in enumerate(list_of_difference):
    
    # Skip positions that have no changes
    if count[i] == 0:
        continue
        
    else:
        # Calculate stderror and z value
        stderror = std / sqrt(count[i])
        z = (difference - 0)/stderror
        
        # p < 0.001, so significant z values are z < -3.29
        # If this dN - dS is significant, add it to the dictionary of positive selection positions
        if z < -3.29:
            pos_selection[i] = np.log(list_of_ratios[i])
        # Otherwise, add it to the list of other positions
        else:
            other[i] = np.log(list_of_ratios[i])



# PART 4
# Plot dN/dS ratios vs. codon position

fig, ax = plt.subplots(figsize=(20, 8))
ax.scatter( pos_selection.keys(), pos_selection.values(), 
            alpha=1, s=6, color="mediumvioletred", label="Positive selection"
            )
ax.scatter( other.keys(), other.values(), 
            alpha=1, s=6, color="royalblue", label="Negative or no selection"
            )
ax.set_xlabel("Codon position")
ax.set_ylabel("log(dN/dS ratio)")
ax.set_title("Ratio of nonsynonymous to synonymous changes at each codon")
ax.legend(bbox_to_anchor=(1.01,0.52), loc=2, borderaxespad=0.)
plt.tight_layout()
fig.savefig("dN_dS.png")
plt.close(fig)


