#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    f = open(sys.argv[1])
else:
    f = sys.stdin
    # command < file

count = 0
total = 0

for line in f:
    if line.startswith("SRR"):
        count += 1
        cut_line = line.strip().split("\t")
        total += int(cut_line[4])
        #print(total)

print("Number of lines is: " + str(count))
print("Total of all MAPQ scores is: " + str(total))

avg_mapq = total / count
print("Avg MAPQ score is: " + str(avg_mapq))