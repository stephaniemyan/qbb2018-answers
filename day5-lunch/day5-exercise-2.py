#!/usr/bin/env python3

"""
Usage: ./day5-exercise-2.py <ctab_file>

<ctab_file> Filepath to .ctab file

This script takes a .ctab file of your choice, and approximates the promoter region
for each transcript by returning the region +/- 500bp from the transcription start
site.

It outputs a .bed file with columns: chromosome, promoter start, promoter end, 
transcript name.

"""

import sys
import pandas as pd

# Open file
df = pd.read_table(sys.argv[1])


# Read .ctab file to find start site for each transcript
# Create list of promoter start sites (start - 500) and end sites (start + 500)
starts_col = []
ends_col = []

# Take out strand (+/-), start site, and end site columns from dataframe
strands = df.loc[:,"strand"]
starts = df.loc[:,"start"]
ends = df.loc[:,"end"]

for i, item in enumerate(strands):
    # If transcript is on reverse strand, it actually begins at the "end" site
    if item == "-":
        if ends[i] <= 499:
            # bigWigAverageOverBed doesn't take negative start site #s, so we pre-emptively set all negatives to 0
            starts_col.append(0)
            ends_col.append(0)
        else:
            starts_col.append(ends[i] + 500)
            ends_col.append(ends[i] - 500)
    # The positive strand starts at the "start" site
    else:
        if starts[i] <= 499:
            # bigWigAverageOverBed doesn't take negative start site #s, so we pre-emptively set all negatives to 0
            starts_col.append(0)
            ends_col.append(0)
        else:
            starts_col.append(starts[i] - 500)
            ends_col.append(starts[i] + 500)


# Create new dataframe with chromosome, promoter start/end, and transcript name
new_df = pd.DataFrame(
    {"chr" : df.loc[:,"chr"],
     "start" : starts_col,
     "end" : ends_col,
     "t_name" : df.loc[:,"t_name"]
    })

# Output dataframe without index #s or header, for ease of use with bigWig
# You will have to pipe the output to a .bed file yourself
new_df.to_csv(sys.stdout, sep="\t", index=False, header=False)