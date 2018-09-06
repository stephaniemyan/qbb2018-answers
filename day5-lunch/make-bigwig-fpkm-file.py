#!/usr/bin/env python3

"""
Usage: ./make-bigwig-fpkm-file.py <ctab_file> <tab_file 1> ... <tab_file n>

<ctab_file> Filepath to .ctab file with FPKMs
<tab_file 1 ... n> Filepaths to .tab file with bigWigAverageOverBed outputs

This script creates a tab-delimited file with FPKMs and their corresponding bigWig
means ('average over just covered bases', as defined by the bigWig usage page).
These means indicate how much of that histone modification is present in the promoter
of a particular transcript.

The script takes the FPKM and bigWig info from a .ctab file and the bigWig .tab
output files, and organizes it all into one file.
"""

import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm


# Read .ctab file, ordered by transcript name
ctab_df = pd.read_table(sys.argv[1], index_col="t_name")

# Make a new dataframe with just FPKMs from the ctab file
df_new = pd.DataFrame({"FPKM" : ctab_df.loc[:,"FPKM"]})


# For every bigWig output file provided in the command line, add a column with the means of the histone modification coverage to the new dataframe

for item in sys.argv[2:len(sys.argv)]:
    
    # Extract histone modification name from bigWig file name
    hisname = item.split(".")[-2]
    
    # Convert bigWig file into a dataframe, manually adding headers
    his_df = (pd.read_csv(item,
             sep="\t",
             header=None,
             index_col="t_name",
             names = ["t_name","size","covered","sum","mean0","mean"]))
    
    # Sort bigWig dataframe by transcript name, so it matches the order of the FPKMs
    his_df = his_df.sort_values("t_name")
    
    # Add column with histone coverage means to the FPKM dataframe
    df_new = df_new.assign(name = his_df.loc[:,"mean"])
    # Rename new column with the name of the histone modication
    df_new = df_new.rename(columns={"name": str(hisname)})

# You will have to pipe the output to a .tab file yourself
df_new.to_csv(sys.stdout, sep="\t")