#!/usr/bin/env python3

# ident_map.py maps FlyBase IDs from a c_tab file to corresponding UniProt IDs

# USE THIS SYNTAX IN COMMAND LINE:
# ident_map.py <.txt file> <.c_tab file> [-m]
# Command line input 1: .txt file with UniProt-FlyBase correspondences
# Command line input 2: c_tab file
# Command line input 3: (OPTIONAL) -m will print "No match" if c_tab FlyBase ID has no corresponding UniProt ID

import sys

# Open files
map_file = open(sys.argv[1])
ctab_file = open(sys.argv[2])

# Convert .txt file into dict with FlyBase-UniProt correspondences

dict = {}

for line in map_file:
    fields = line.strip().split()
    key = fields[0]
    val = fields[1]
    dict[key] = val
    
# Match FlyBase ID in c_tab file to corresponding UniProt ID in the dict
# Print out FlyBase IDs that have corresponding UniProt IDs
# If user includes -m input, print "No match" if c_tab FlyBase ID has no corresponding UniProt ID
    
for line in ctab_file:
    fields = line.strip().split("\t")
    fly = fields[8]
    if fly in dict:
        unip = dict[fly]
        print(fly + "\t" + unip)
    elif len(sys.argv) > 3 and sys.argv[3] == "-m":
        print(fly + "\t" + "No match")