#!/usr/bin/env python3

"""
Usage: ./day4-homework-1.py <samples.csv> <ctab_dir>

<t_name> Transcript of interest
<samples.csv> Path to samples file
<ctab_dir> Path to directory with .ctab files

Output all FPKMs from all samples in samples.csv
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

# Read file with list of all samples
df = pd.read_csv(sys.argv[1])

# Grab all FPKMs from the c_tab file corresponding to your current sample
# Make a dictionary entry
# Key: Sample sex/stage; value = list of FPKMS
all_dict = {}
for index, sample, sex, stage in df.itertuples():
    filename = os.path.join(sys.argv[2], sample, "t_data.ctab")
    ctab_fpkm = pd.read_table(filename, index_col="t_name" ).loc[:,"FPKM"]
    all_dict[str(sex) + "_" +str(stage)] = ctab_fpkm
df = pd.DataFrame(all_dict)
df.to_csv(sys.stdout)