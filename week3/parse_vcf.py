#!/usr/bin/env python3

import sys

f = open(sys.argv[1])

for line in f:
    if line.startswith("##"):
        continue
    else:
        print(line)