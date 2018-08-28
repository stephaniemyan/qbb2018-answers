#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    f = open(sys.argv[1])
else:
    f = sys.stdin
    # command < file

count = 0
for line in f:
    if line.startswith("SRR"):
        cut_line = line.strip().split("\t")
        if int(cut_line[1]) & (1 << 4) > 0:
            count += 1
print(count)