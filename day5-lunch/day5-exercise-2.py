#!/usr/bin/env python3

"""
Usage: ./day5-exercise-2.py <ctab_file>

<ctab_file> Filepath to ctab file

"""

import sys
import pandas as pd

# Select for specific columns
df = pd.read_table(sys.argv[1])

starts_col = df.loc[:,"start"] - 500
ends_col = df.loc[:,"start"] + 500

new_df = pd.DataFrame(
    {"chr" : df.loc[:,"chr"],
     "start" : starts_col,
     "end" : ends_col,
     "t_name" : df.loc[:,"t_name"]
    })
    
new_df.to_csv(sys.stdout, sep="\t")

#new_df = pd.DataFrame(all_dict)


#tarts_col = df.loc[:,"start"].assign(start=sum_column)
#ends_col = df.loc[:,"start"].assign(start=sum_column)

#df.to_csv(sys.stdout, sep="\t")


# d = {"first" : ctab1,
#      "second" : ctab2 }
#
# df = pd.DataFrame(d)
# df.to_csv(sys.stdout, sep="\t")