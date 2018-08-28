#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    f = open(sys.argv[1])
else:
    f = sys.stdin
    # command < file

for line in f:
    if line.startswith("SRR"):
        for i, line in enumerate(f):
            if i < 10:
                cut_line = line.strip().split("\t")
                print("Read 1: " + cut_line[2])
                i += 1