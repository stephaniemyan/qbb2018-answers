#!/usr/bin/env python3

"""
Usage: ./merge_fpkms.sh <threshold> <ctab_file1> <ctab_file2> ... <ctab_filen>

Only report transcript when total FPKM > threshold
"""

import sys
import os
import pandas as pd

num_files = len(sys.argv) - 2
list_of_files = sys.argv[2:len(sys.argv)]
threshold = float(sys.argv[1])

dict = {}
for file in list_of_files:
    name = file.split(os.sep)[-2]
    fpkms = pd.read_csv(file, sep="\t", index_col="t_name").loc[:,"FPKM"]
    dict[name] = fpkms
fpkms_all = pd.DataFrame(dict)
#print(fpkms_all)

sum_column = fpkms_all.sum(axis=1)
fpkms_sums = fpkms_all.assign(sums=sum_column)
over_thresh_bool = fpkms_sums.loc[:, "sums"] > threshold

print(fpkms_sums.loc[over_thresh_bool,:])