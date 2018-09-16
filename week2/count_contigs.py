#!/usr/bin/env python3

"""
Usage: ./blast_to_fasta.py <contigs.fa>
"""

import sys
import fasta

contigs_fa = fasta.FASTAReader(open(sys.argv[1]))

def do_summary(fasta_file):
    num_contigs = 0
    contig_lengths = []
    
    for ident, sequence in fasta_file:
        num_contigs += 1
        contig_lengths.append(len(sequence))
        
    contig_lengths.sort(reverse=True)
    print("Total # of contigs = " + str(num_contigs))
    print("Avg contig length = " + str(sum(contig_lengths)/num_contigs))
    print("Min contig length = " + str(contig_lengths[-1]))
    print("Max contig length = " + str(contig_lengths[0]))
    
    half_length = sum(contig_lengths)/2
    counter = 0
    for i, length in enumerate(contig_lengths):
        if counter >= half_length:
            print("1/2 total sequence length = " + str(half_length))
            print("Current position in sequence = " + str(counter))
            print("n50 = " + str(contig_lengths[i-1]))
            break
        else:
            counter += length

do_summary(contigs_fa)