#!/usr/bin/env python3

"""
Usage: ./mean_transcripts.py <samples.csv> <ctab_dir> <gene_name1> <gene_name2> ... <gene_namen>

<samples.csv> - Path to samples file
<ctab_dir> - Path to directory with .ctab files
<gene_name 1 ... n> Arbitrary # of genes of interest

Create timecourses of FPKMs for n genes of interest and plot each vs. stage
WIll plot genders separately
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Create list of genes of interest by looping over command line arguments
gene_names = []
for item in sys.argv[3:]:
    gene_names.append(item)
samples_file = sys.argv[1]
ctab_file = sys.argv[2]


# mean_transcripts returns the mean of all FPKMs for a gene of interest, within one sample
# mean_transcripts takes gender and gene name as arguments, to output different results for male/female samples and different genes
def mean_transcripts(gender, gene):
    
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
        roi = ctab_df.loc[:,"gene_name"] == gene
        all_fpkms.append(ctab_df.loc[roi,"FPKM"])
        mean_fpkms.append(np.mean(all_fpkms))
    return mean_fpkms
    
# Run mean_transcripts and plot for every name in list of gene names
for name in gene_names:
    mean_f = mean_transcripts("female", name)
    mean_m = mean_transcripts("male", name)
    
    fig, ax = plt.subplots()
    ax.plot(mean_f, label="Female")
    ax.plot(mean_m, label="Male")

    ax.set_title(str(name))
    ax.set_xlabel("developmental stage")
    ax.set_ylabel("mean FPKM of transcript")

    my_xticks = ["9", "10","11","12","13","14A","14B","14C","14D"]
    ax.set_xticklabels(labels = my_xticks)
    plt.xticks(rotation=90)

    plt.legend(bbox_to_anchor=(1.07,0.55), loc=2, borderaxespad=0.)

    plt.tight_layout()
    fig.savefig(str(name) + "_mean_transcripts.png", bbox_inches='tight')
plt.close(fig)
    
    
    
    
    
    