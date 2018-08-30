#!/usr/bin/env python3

"""
Usage: ./day4-homework-2.py <t_name> <samples.csv> <ctab_dir>

samples - name of sample files
ctab_dir - location of directory with files
Create a timecourse of a given transcript (FBtr0331261) for females
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

transcript = sys.argv[1]
samples_file = sys.argv[2]
ctab_file = sys.argv[3]
replicates_file = sys.argv[4]

    
def timecourse_samp(gender):
    
    df = pd.read_csv(samples_file)
    soi = df.loc[:,"sex"] == str(gender)
    df = df.loc[soi,:]

    fpkms = []
    for index, sample, sex, stage in df.itertuples():
        filename = os.path.join(ctab_file, sample, "t_data.ctab")
        ctab_df = pd.read_table(filename, index_col="t_name")
        fpkms.append(ctab_df.loc[transcript,"FPKM"])
    return fpkms
    
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

fpkms_m = timecourse_samp("male")
fpkms_f = timecourse_samp("female")

reps_m = timecourse_rep("male")
reps_f = timecourse_rep("female")   


fig, ax = plt.subplots()
ax.plot(fpkms_f, label="Female", color="orange")
ax.plot(fpkms_m, label="Male", color="blue")
ax.plot([4,5,6,7], reps_f, label="Replicate Female", color="red")
ax.plot([4,5,6,7], reps_m, label="Replicate Male", color="green")

ax.set_title("FBtr0331261")
ax.set_xlabel("developmental stage")
ax.set_ylabel("mRNA abundance (FPKM)")

my_xticks = ["9", "10","11","12","13","14A","14B","14C","14D"]
ax.set_xticklabels(labels = my_xticks)
plt.xticks(rotation=90)
#ax.set_yticklabels(labels = np.arange(0,120,20))

plt.legend(bbox_to_anchor=(1.07,0.55), loc=2, borderaxespad=0.)

#plt.show()

plt.tight_layout()
fig.savefig("timecourse.png", bbox_inches='tight')
plt.close(fig)




