#!/usr/bin/env python3

"""
Usage: ./blast_to_fasta.py < <blast_output>

<blast_output> Filepath to BLAST output file. This should be a tabular file with ONLY 
two fields: the subject seq-id (sseqid) and the aligned part of the subject sequence 
(sseq).

Here is the command line syntax to get this output:
  blastn -db nr -query week1_query.fa -evalue .0001 -num_alignments 1000 -ungapped -outfmt \
  "6 sseqid sseq" -out blast.out -remote

- -db nr: Searches only the non-redundant BLAST database
- -evalue .0001: Limits results to only E-values < .0001
- -num_alignments 1000: Limits results to only top 1000 alignments
- -ungapped: Makes sure all alignments are ungapped
- -outfmt "6 sseqid sseq": Formats output as a tabular file (6), with subject seq-id (sseqid)
                           and subject sequence (sseq)
- -remote: Searches against the NCBI database

The script takes this BLAST output and reformats it as a FASTA file where each header is the
subject seq-id.

It outputs directly to the terminal and you will have to pipe it to a file yourself.
"""

import sys

# Split each line by tab and strip out line breaks
# Print >subject sequence ID, line break, and subject sequence
# You will have to pipe the output to a file destination yourself
for line in sys.stdin:
    split_line = line.rstrip("\r\n").split("\t")
    print(">" + split_line[0] + "\n" + split_line[1] + "\n")
    
    
# To translate the new BLAST.fa file to amino acids, use transeq from emboss in command line:
# transeq <BLAST.fa> <output name>

# To align amino acid sequences, use MAFFT:
# mafft <translated_blast_file.fa> > <output name>