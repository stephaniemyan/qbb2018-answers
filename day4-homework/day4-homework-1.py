#!/usr/bin/env python3

"""
Usage: ./day4-homework-1.py <t_name> <samples.csv> <ctab_dir>

samples.csv - path for samples file
ctab_dir - location of directory with corresponding ctab files

Output all FPKMs from all samples in samples.csv, ordered by transcript name
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(sys.argv[1])
# This is the samples.csv file

all_dict = {}
for index, sample, sex, stage in df.itertuples():
    filename = os.path.join(sys.argv[2], sample, "t_data.ctab")
    ctab_fpkm = pd.read_table(filename, index_col="t_name" ).loc[:,"FPKM"]
    all_dict[str(sex)+str(stage)] = ctab_fpkm
df = pd.DataFrame(all_dict)
df.to_csv(sys.stdout)