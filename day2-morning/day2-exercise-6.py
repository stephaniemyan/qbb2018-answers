#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    f = open(sys.argv[1])
else:
    f = sys.stdin
    # command < file

count = 0
for line in f:
    if line.startswith("SRR") and "2L" in line:
        cut_line = line.strip().split("\t")
        if int(cut_line[3]) >= 10000 and int(cut_line[3]) <= 20000:
            count += 1
print(count)