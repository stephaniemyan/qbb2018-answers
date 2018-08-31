#!/usr/bin/env python3

"""
Usage: ./ma_plot.py <ctab_file1> <ctab_file2>

Create an MA plot comparing FPKMs of the same transcript in two different .c_tab files
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ctab1 = sys.argv[1]
ctab2 = sys.argv[2]

name1 = str(ctab1.split(os.sep)[-2])
name2 = str(ctab2.split(os.sep)[-2])

df1 = pd.read_csv(ctab1, sep="\t", index_col="t_name")
df2 = pd.read_csv(ctab2, sep="\t", index_col="t_name")

fpkms1 = df1.loc[:, "FPKM"]
fpkms2 = df2.loc[:, "FPKM"]

# log_df1 = np.log(fpkms1 + 1)
# log_df2 = np.log(fpkms2 + 1)

m = np.log2((fpkms1+1)/(fpkms2+1))
a = 0.5 * np.log2((fpkms1+1) * (fpkms2+ 1))

fig, (ax) = plt.subplots()
ax.scatter(a, m, alpha = 0.1, s=2, color="purple")
ax.set_title("Comparative gene expression in samples SRR072893 and SRR072915")
ax.set_xlabel("avg intensity (log2(FPKM1/FPKM2))")
ax.set_ylabel("intensity ratio 1/2(log2(FPKM1*FPKM2))")
fig.savefig("ma_plot.png")
plt.close()















