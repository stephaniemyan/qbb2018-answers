#!/usr/bin/env python3

# kmer_matcher.py takes a known FASTA sequence and aligns kmers from a query sequence

# USE THIS SYNTAX IN COMMAND LINE:
# kmer_matcher.py <target.fa> <query.fa> <k>
# Command line input 1: Target sequence as a FASTA file
# Command line input 2: Query sequence as a FASTA file
# Command line input 3: Desired k value as integer



import sys
import fasta

# Set k value according to command line input
k = int(sys.argv[3])



# PART 1
# Take target file and create dictionary of kmers in target sequence, with their corresponding positions
# key: k-mer nucleotides; value: start location in fasta file

# Read target file, run through FASTA reader
target_file = open(sys.argv[1])
target = fasta.FASTAReader(target_file)

# Generate target kmers from target sequence
# Make dictionary of kmers and their positions
# Key: kmer sequence; value: kmer start position
kmer_dict = {}
for t_ident, t_sequence in target:
    for i in range(0, len(t_sequence) - k):
        kmer = t_sequence[i:i+k]
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [i]
        else:
            kmer_dict[kmer].append(i)



# PART 2
# Look up kmers from query file in dictionary of target kmers

# Read query file, run through FASTA reader
query_file = open(sys.argv[2])
query = fasta.FASTAReader(query_file)

# Generate query kmers from query sequence
# Look up query kmers in target dictionary
# If query kmer is in dictionary, print out relevant info
# t_ident is the target sequence name, and comes from the fasta.py function
for q_ident, q_sequence in query:
    for i in range(0, len(q_sequence) - k):
        query_kmer = q_sequence[i:i+k]
        if query_kmer in kmer_dict:
            print(
                "Target sequence name: " + t_ident + "\n"
                "Target sequence start(s): " + ", ".join(str(x) for x in kmer_dict[query_kmer]) + "\n"
                "Query sequence start: " + str(i)+ "\n"
                "k-mer: " + query_kmer + "\n"
            )