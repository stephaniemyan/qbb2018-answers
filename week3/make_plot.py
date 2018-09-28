#!/usr/bin/env python3

"""
Usage: ./make_plot.py <VCF file 1> <VCF file 2>

<VCF file 1> An annotated VCF file, output from snpEff
<VCF file 2> An annotated VCF file, output from snpEff, with GQ values. I had to do this because
my original VCF file didn't include genotype quality scores, so I had to hard-code this extra
file in

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


# Read .vcf files as dataframes
df = pd.read_table(sys.argv[1], sep="\t")
df2 = pd.read_table(sys.argv[2], sep="\t")

# Pull out INFO and FORMAT columns
info_col = df.loc[:,"INFO"]
# col_09 = df2.loc[:,"A01_09"]
# col_11 = df2.loc[:,"A01_11"]
# col_23 = df2.loc[:,"A01_23"]
# col_24 = df2.loc[:,"A01_24"]
# col_27 = df2.loc[:,"A01_27"]
# col_31 = df2.loc[:,"A01_31"]
# col_35 = df2.loc[:,"A01_35"]
# col_39 = df2.loc[:,"A01_39"]
# col_62 = df2.loc[:,"A01_62"]
# col_63 = df2.loc[:,"A01_63"]

col_09 = df2.loc[:,"A1_09"]
col_11 = df2.loc[:,"A1_11"]
col_23 = df2.loc[:,"A1_23"]
col_24 = df2.loc[:,"A1_24"]
col_27 = df2.loc[:,"A1_27"]
col_31 = df2.loc[:,"A1_31"]
col_35 = df2.loc[:,"A1_35"]
col_39 = df2.loc[:,"A1_39"]
col_62 = df2.loc[:,"A1_62"]
col_63 = df2.loc[:,"A1_63"]


# Create a multi-panel figure
fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(20, 10))
# flatten lets you call each plot in the figure by a number
axes = axes.flatten()


# Produce lists of DP values (read depth) and GQ values (genotype quality) from the .vcf file
# These values are listed for each A01_ alignment(?)'s FORMAT column
# I didn't know what to do to reconcile all of them, so I added all 10 DP/GQ values to my lists

# Create DP and GQ lists
dps = []
gqs = []

# Loop through each of the the A01_ columns
for coi in [col_09, col_11, col_23, col_24, col_27, col_31, col_35, col_39, col_62, col_63]:
    
    # Split every line by colons (since it should be formatted DP:GQ:etc)
    for line in coi:
        if ":" in line:
            fields = line.split(":")
            # GQ should be the first field, and DP should be the second field
            gq = float(fields[1])
            dp = float(fields[2])
            dps.append(dp)
            gqs.append(gq)
            

# Graph read depth values as a histogram
axes[0].hist(dps, color="royalblue", bins=100)
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


# Produce a list of allele frequency values (AF), which is in the INFO column
af_values = []

id, val = fields.split("=")

for line in info_col:
    fields = line.split(";")
    id, val = fields.split("=")
    if id == "DP": etc
    af = fields[3].split("=")[1]
    if "," in af:
        afs = af.split(",")
        for i in range(len(afs)):
            af_values.append(float(afs[i]))
    else:
        af_values.append(float(af))

# Graph allele frequency values as a histogram
axes[2].hist(af_values, color="lightcoral", bins=30)
axes[2].set_ylabel("Number of variants")
axes[2].set_xlabel("Allele frequency")
axes[2].set_title("Distribution of allele frequencies among identified variants")


info = dict([x.split("=") for x in fields[7].split(";")])
af = float(info["AF"])
letters = [x for x in "Peter"]
letters = ["P" "e" "t" "e" "r"]
# I hard-coded all these values from the snpEff_summary.html document that snpEFF outputs
# I realize we were supposed to parse these ourselves from the .vcf but I can't handle that right now
# Let's be real, if you were doing this irl you would also be getting these #s from the html file

y_values = [42.313, 10.568, 3.972, 0.082, 0.0, 0.0, 0.011, 43.053]
x_values = ["", "Downstream", "Exon", "Intergenic", "Intron", "Acceptor",
             "Donor", "Splice Site", "Upstream"]
             
# Plot variant frequencies in a barplot
axes[3].bar(np.arange(8), y_values, color="mediumvioletred")
axes[3].set_ylabel("Frequency (%)")
axes[3].set_xticklabels(x_values)
axes[3].set_xlabel("Type of variant")
axes[3].set_title("Predicted effects of each variant")


plt.savefig("ugh.png")