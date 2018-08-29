#!/usr/bin/env python3

# uniprot_map.py prints out a table with FlyBase IDs and corresponding UniProt IDs
# You will have to pipe output into a .txt file yourself

# USE THIS SYNTAX IN COMMAND LINE:
# uniprot_map.py <.txt file>
# Command line input 1: .txt file that includes some kind of table matching FlyBase IDs and UniProt IDs

import sys

# Open txt file
f = open(sys.argv[1])

# 'Grep' lines that include "DROME"
# Select for lines with "FBGn", which indicates that they include a FlyBase ID
# Print everything out

for line in f:
    if "DROME" in line:
        fields = line.strip().split()
        if fields[-1].startswith("FBgn"):
            print(fields[3], fields[2])