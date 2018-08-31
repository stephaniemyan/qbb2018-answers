#!/usr/bin/env python3

"""
Usage: ./day5-exercise-4.py <ctab_file> <tab_file1-5

<ctab_file> Filepath to ctab file

"""

import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm

#MAKE OUR NEW TABLE

ctab_df = pd.read_table(sys.argv[1], index_col="t_name")
df_new = pd.DataFrame({"FPKM" : ctab_df.loc[:,"FPKM"]})

for item in sys.argv[2:len(sys.argv)]:
    
    hisname = item.split(".")[-2]
    his_df = (pd.read_csv(item,
             sep="\t",
             header=None,
             index_col="t_name",
             names = ["t_name","size","covered","sum","mean0","mean"]))
    his_df = his_df.sort_values("t_name")
    #his_df = pd.DataFrame({"hello" : his_df.loc[:,"mean"]})
    df_new = df_new.assign(name = his_df.loc[:,"mean"])
    df_new = df_new.rename(columns={"name": str(hisname)})

df_new.to_csv(sys.stdout, sep="\t")