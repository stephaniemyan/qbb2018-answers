#!/usr/bin/env python3

"""
Usage: ./make_af.py <VCF file>

<VCF file> Output of freebayes with extensive filtering. The only information in the INFO 
column should be the AF value(s).

Makes a histogram of allele frequency values from a freebayes VCF output file.
"""

import sys
import matplotlib.pyplot as plt
plt.style.use('ggplot')

# Open file
f = open(sys.argv[1])

# Parse file to extract AF values
af_list = []
for line in f:
    # Skip header lines
    if line.startswith("#"):
        continue
    fields = line.rstrip("\r\n").split("\t")
    # the 8th field is INFO. Split it by the = and take out the values, which are on the right side of the = sign
    af = (fields[7].split("="))[1]
    if "," in af:
        afs = af.split(",")
        # Account for multiple AF values in one line by appending all of them to the list
        for i in range(len(afs)):
            af_list.append(float(afs[i]))
    else:
        af_list.append(float(af))

# Make histogram
fig, ax = plt.subplots()
ax.hist(af_list, color="royalblue", bins=100)
ax.set_xlabel("Allele frequency")
ax.set_ylabel("Number of variants")
ax.set_title("Distribution of allele frequencies among individual variants")
plt.tight_layout()
fig.savefig("af.png")
plt.close(fig)