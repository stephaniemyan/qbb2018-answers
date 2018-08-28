#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    f = open(sys.argv[1])
else:
    f = sys.stdin
    # command < file

count = 0
for line in f:
    if "NM:i:0" in line:
        count += 1
print(count)