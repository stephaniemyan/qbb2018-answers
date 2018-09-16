#!/usr/bin/env python3

"""
Usage: ./count_contigs.py <contigs.fa>

<contigs.fa> FASTA-formatted output of a contig assembler

Takes a contig assembly output file and returns total # of contigs; avg, min, and max contig
lengths; and n50.
"""

import sys
import fasta

# Read contigs file
contigs_file = fasta.FASTAReader(open(sys.argv[1]))

# Calculate # of contigs and min, max, and avg contig length
# Set counter for # of contigs and a list to record contig lengths
num_contigs = 0
contig_lengths = []

for ident, sequence in contigs_file:
    
    # Increment contig counter for each new contig
    num_contigs += 1
    # Add contig length to list
    contig_lengths.append(len(sequence))

# Sort list of lengths (from largest to smallest, for the n50 calculation)
contig_lengths.sort(reverse=True)
print("Total # of contigs = " + str(num_contigs))
print("Avg contig length = " + str(sum(contig_lengths)/num_contigs))
print("Min contig length = " + str(contig_lengths[-1]))
print("Max contig length = " + str(contig_lengths[0]))


# Calculate n50
# At least 50% of the genome is covered by contigs of minimum length n50
# If the x longest contigs are the minimum requirement to cover 50% of the genome, n50 is the length of the xth contig

# Calculate 50% of length of the genome, set counter for current position
half_length = sum(contig_lengths)/2
counter = 0

for i, length in enumerate(contig_lengths):
    
    # If the current position is greater than the 50% position, print out the previous contig length
    if counter >= half_length:
        print("1/2 total sequence length = " + str(half_length))
        print("Current position in sequence = " + str(counter))
        print("n50 = " + str(contig_lengths[i-1]))
        break
    # Otherwise, add current contig length to the position count
    else:
        counter += length