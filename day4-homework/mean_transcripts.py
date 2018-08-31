#!/usr/bin/env python3

"""
Usage: ./mean_transcripts.py <gene_name> <samples.csv> <ctab_dir>

<gene_name> Name of gene of interest
<samples.csv> - Path to samples file
<ctab_dir> - Path to directory with .ctab files

Create a timecourse of FPKMs for a given gene transcript and plot vs. stage
Will plot genders separately
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

gene_name = sys.argv[1]
samples_file = sys.argv[2]
ctab_file = sys.argv[3]


# mean_transcripts returns the mean of all FPKMs for a gene of interest, within one sample
# mean_transcripts takes gender as an argument to output different results for male/female samples
def mean_transcripts(gender):
    
    df = pd.read_csv(samples_file)
    soi = df.loc[:,"sex"] == str(gender)
    df = df.loc[soi,:]

    all_fpkms = []
    mean_fpkms = []
    # Read .ctab file for sample
    for index, sample, sex, stage in df.itertuples():
        filename = os.path.join(ctab_file, sample, "t_data.ctab")
        ctab_df = pd.read_table(filename)
    
    # Identify rows that have the correct gene name
    # Extract FPKMs for those rows
    # Take mean of FPKMs, and add it to a mean_fpkms - a list that keeps track of them for all samples
        roi = ctab_df.loc[:,"gene_name"] == gene_name
        all_fpkms.append(ctab_df.loc[roi,"FPKM"])
        mean_fpkms.append(np.mean(all_fpkms))
    return mean_fpkms
    
    
# Run mean_transcripts to generate list of FPKMs for male and female samples
mean_f = mean_transcripts("female")
mean_m = mean_transcripts("male")


# Plot mean_transcripts results
fig, ax = plt.subplots()
ax.plot(mean_f, label="Female")
ax.plot(mean_m, label="Male")

ax.set_title(str(gene_name))
ax.set_xlabel("developmental stage")
ax.set_ylabel("mean FPKM of transcript")

my_xticks = ["9", "10","11","12","13","14A","14B","14C","14D"]
ax.set_xticklabels(labels = my_xticks)
plt.xticks(rotation=90)

plt.legend(bbox_to_anchor=(1.07,0.55), loc=2, borderaxespad=0.)

plt.tight_layout()
fig.savefig(str(gene_name) + "_mean_transcripts.png", bbox_inches='tight')
plt.close(fig)

    