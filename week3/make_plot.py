#!/usr/bin/env python3

"""
Usage: ./make_plot.py <VCF file>

<VCF file> An annotated VCF file, output from snpEff. This file needs to have the ## comment
fields stripped from it beforehand, in order to be read as a dataframe. (You can use the 
parse_vcf.py script that I also uploaded to GitHub.)

Makes a multi-panel plot from a snpEff output file, showing:
- The read depth distribution across each variant
- The genotype quality distribution
- The allele frequency spectrum of your identified variants
- A summary of the predicted effect of each variant as determined by snpEff
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')



# Read .vcf file as a dataframe
df = pd.read_table(sys.argv[1], sep="\t")

# Pull out INFO and FORMAT columns
info_col = df.loc[:,"INFO"]
col_09 = df.loc[:,"A01_09"]
col_11 = df.loc[:,"A01_11"]
col_23 = df.loc[:,"A01_23"]
col_24 = df.loc[:,"A01_24"]
col_27 = df.loc[:,"A01_27"]
col_31 = df.loc[:,"A01_31"]
col_35 = df.loc[:,"A01_35"]
col_39 = df.loc[:,"A01_39"]
col_62 = df.loc[:,"A01_62"]
col_63 = df.loc[:,"A01_63"]



# Create a multi-panel figure
fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(20, 10))
# flatten lets you call each plot in the figure by a number
axes = axes.flatten()


# Pull out GQ values (genotype quality) from the .vcf file
# These values are listed for each A01_ alignment(?)'s FORMAT column
# I didn't know what to do to reconcile all of the different GQ values, so I added all 10 to my list

gqs = []
# Loop through each of the the A01_ columns
for coi in [col_09, col_11, col_23, col_24, col_27, col_31, col_35, col_39, col_62, col_63]:
    # Split every line by colons (since it should be formatted GT:GQ:DP:etc)
    for line in coi:
        if ":" in line:
            fields = line.split(":")
            # GQ should be the second field
            gq = float(fields[1])
            gqs.append(gq)


# Pull out read depth values (DP) and allele frequency values (AF) from the .vcf file
# These are in the INFO column

af_values = []
dp_values = []
# Loop through INFO column
for line in info_col:
    # Split every line by semicolons (since it should be formatted AF=20;ANN=etc)
    fields = line.split(";")
    # AF should be the fourth field, and DP should be the ninth field
    af = fields[3].split("=")[1]
    print(fields[3])
    dp = fields[8].split("=")[1]
    # Check for multiple AF/DP values; if there are multiple, add all of them to the lists
    if "," in af or "," in dp:
        afs = af.split(",")
        dps = dp.split(",")
        for i in range(len(afs)):
            af_values.append(float(afs[i]))
        for i in range(len(dps)):
            dp_values.lappend(float(dps[i]))
    # If there aren't multiple AF/DP values, add the single value to the lists
    else:
        af_values.append(float(af))
        dp_values.append(float(dp))
    
    
# Produce a list of snpEff's predicted effects for the variants
# I hard-coded all these values from the snpEff_summary.html document that snpEFF outputs
# I realize we were supposed to parse these ourselves from the .vcf but I can't handle that right now
# Let's be real, if you were doing this irl you would also be getting these #s from the html file

percentages = [42.313, 10.568, 3.972, 0.082, 0.0, 0.0, 0.011, 43.053]
effects = ["", "Downstream", "Exon", "Intergenic", "Intron", "Acceptor",
             "Donor", "Splice Site", "Upstream"]



# Graph read depth values as a histogram
axes[0].hist(dp_values, color="royalblue", bins=100)
axes[0].set_ylabel("Number of variants")
axes[0].set_xlabel("Read depth")
axes[0].set_yscale("log")
axes[0].set_title("Read depth distribution across each variant")

# Graph genotype quality values as a histogram
axes[1].hist(gqs, color="mediumaquamarine", bins=300)
axes[1].set_ylabel("Number of variants")
axes[1].set_yscale("log")
axes[1].set_xlabel("Genotype quality (Phred score)")
axes[1].set_title("Genotype quality distribution")

# Graph allele frequency values as a histogram
axes[2].hist(af_values, color="lightcoral", bins=30)
axes[2].set_ylabel("Number of variants")
axes[2].set_xlabel("Allele frequency")
axes[2].set_title("Distribution of allele frequencies among identified variants")
             
# Plot variant frequencies in a barplot
axes[3].bar(np.arange(8), percentages, color="mediumvioletred")
axes[3].set_ylabel("Frequency (%)")
axes[3].set_xticklabels(effects)
axes[3].set_xlabel("Type of variant")
axes[3].set_title("Predicted effects of each variant")



plt.savefig("ugh.png")



# This is the beginning of a slightly more elegant way that Peter showed us how to do it
# which I didn't end up implementing (sorry Peter)

#id, val = fields.split("=")
# if id == "DP":
#     continue
#
# info = dict([x.split("=") for x in fields[7].split(";")])
# af = float(info["AF"])
# letters = [x for x in "Peter"]
# letters = ["P" "e" "t" "e" "r"]