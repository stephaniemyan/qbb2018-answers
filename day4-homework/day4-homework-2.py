#!/usr/bin/env python3

"""
Usage: ./day4-homework-2.py <t_name> <samples.csv> <replicates.csv> <ctab_dir>

<t_name> Name of transcript of interest
<samples.csv> - Path to samples file
<replicates.csv> Path to replicates file
<ctab_dir> - Path to directory with .ctab files

Create a timecourse of FPKMs for a given gene transcript and plot vs. stage
Plot genders separately, and will also plot replicate samples separately
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


transcript = sys.argv[1]
samples_file = sys.argv[2]
replicates_file = sys.argv[3]
ctab_file = sys.argv[4]


# timecourse_<> extracts FPKM column corresponding to a specific transcript
# timecourse_<> takes gender as an argument to output different results for male/female samples
def timecourse_samp(gender):
    
    df = pd.read_csv(samples_file)
    soi = df.loc[:,"sex"] == str(gender)
    df = df.loc[soi,:]

    fpkms = []
    # Read .ctab file for sample
    for index, sample, sex, stage in df.itertuples():
        filename = os.path.join(ctab_file, sample, "t_data.ctab")
        ctab_df = pd.read_table(filename, index_col="t_name")
        
    # Extract FPKMs corresponding to transcript of interest
        fpkms.append(ctab_df.loc[transcript,"FPKM"])
    return fpkms

# Create separate lists of FPKMs for replicate files
def timecourse_rep(gender):
    
    df = pd.read_csv(replicates_file)
    soi = df.loc[:,"sex"] == str(gender)
    df = df.loc[soi,:]

    fpkms = []
    for index, sample, sex, stage in df.itertuples():
        filename = os.path.join(ctab_file, sample, "t_data.ctab")
        ctab_df = pd.read_table(filename, index_col="t_name")
        fpkms.append(ctab_df.loc[transcript,"FPKM"])
    return fpkms


# Run timecourse to generate list of FPKMs for male and female samples
fpkms_m = timecourse_samp("male")
fpkms_f = timecourse_samp("female")

reps_m = timecourse_rep("male")
reps_f = timecourse_rep("female")   


# Plot timecourse results
fig, ax = plt.subplots()
ax.plot(fpkms_f, label="Female", color="orange")
ax.plot(fpkms_m, label="Male", color="blue")
# Replicates are only plotted along x-ticks 2-7 (i.e. stages 14A-D)
ax.plot([4,5,6,7], reps_f, label="Replicate Female", color="red")
ax.plot([4,5,6,7], reps_m, label="Replicate Male", color="green")

ax.set_title("FBtr0331261")
ax.set_xlabel("developmental stage")
ax.set_ylabel("mRNA abundance (FPKM)")

my_xticks = ["9", "10","11","12","13","14A","14B","14C","14D"]
ax.set_xticklabels(labels = my_xticks)
plt.xticks(rotation=90)

plt.legend(bbox_to_anchor=(1.07,0.55), loc=2, borderaxespad=0.)

plt.tight_layout()
fig.savefig("timecourse.png", bbox_inches='tight')
plt.close(fig)


