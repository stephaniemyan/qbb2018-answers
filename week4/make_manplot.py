#!/usr/bin/env python3

"""
Usage: ./make_manplot.py <.qassoc file 1> ... <.qassoc file n>

<.qassoc file> Any number of .qassoc files output from plink.

Produce a Manhattan plot, showing association of SNPs with a specific phenotype, from a plink 
.qassoc file. Highlight SNPs with p-values less than 10^-5. The chromosomes will be plotted 
out of order because their names are Roman numerals and pandas sorts them alphabetically.

I'd like to thank Elad Joseph on StackOverflow for the dataframe idea, and Rebekka for working
with me. I'd also like to thank Peter because I used a bunch of his dataframe plotting code 
from the week 5 review.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



for f in sys.argv[1:]:
    
    
    # Extract name of the experimental treatment from the file name
    treatment = f.split(".")[1]
    # Read the file as a dataframe
    df = pd.read_table(f, delim_whitespace=True)
    
    
    """
    PART 1
    
    Create a column of index values that lets you plot each chromosome consecutively on the
    same plot.
    
    Peter's pandas method indexes each SNP one after the other, without accounting for the 
    distance between SNPs. That's a lot easier, and it looks pretty much the same as this
    method, but oh well, we coded this already...
    
    Witn this indexing, we space SNPs relative to each other on the chromosome. Also, it 
    indexes consecutive chromosomes right after each other - if chrI is 100bp long and chrII 
    is 200bp long, the indexing for chrII will start at 100 and count up from there. This 
    doesn't account for chromosomes where the SNPs in the VCF file don't cover the entire 
    length of the chromosome, but that's probably not super important.
    """
  
    # Pull out the column that shows each SNP's location on the chromosome
    bp = df.loc[:,"BP"].tolist()

    # positions is a list of each SNP's new index position
    positions = []
    # Count keeps track of which position you're currently at
    count = 0

    # Loop through the list of BP positions
    for i, pos in enumerate(bp):
        # At the very first BP position, create a 0 index value in the positions list
        if i == 0:
            positions.append(0)
        # If your current BP value is less than the previous one, this indicates that a new chromosome has started
        # So you want this chr's indexing to start right after the previous one ended
        # Which means your next index position is just 1 more than the previous one
        elif bp[i] < bp[i-1]:
            positions.append(count + 1)
            count += 1
        # Otherwise add the location for your new SNP
        # The difference in position between the previous SNP and this one should = how far your new SNP is from the previous position
        else:
            difference = bp[i] - bp[i-1]
            positions.append(count + difference)
            count += difference
    
    
    """
    PART 2
    
    Reorganize the dataframe so it has the information you need:
    - -log(p-values)
    - significant p-values (p < 10^-5)
    - grouping by chromosome (for plotting purposes)
    """
        
    # Add a column with -log10 of the p-values
    pvalue = -np.log10(df.loc[:,"P"])
    df = df.assign(pvalue=pvalue, position=positions)
    
    # Add a column with all the significant p-values (p < 10^-5)
    # I'm still not entirely sure how this syntax works, but it works. Thanks Rebekka
    df.loc[df['pvalue'] > 5, 'sigs'] = df['pvalue']
    
    # Group by chromosome name. I don't know why we need to say sort=False, but things are out of order otherwise
    groups = df.groupby('CHR', sort=False)


    """
    PART 3
    
    Make a Manhattan plot from the dataframe by iterating through the groups (i.e. through 
    each individual chromosome), plotting the SNPs and significant SNPs.
    """
    
    fig, ax = plt.subplots(figsize=(20,10))
    
    # Pick two colors to label alternate chromosomes and their significant p-values. Thanks Peter
    colors = ['lightblue', 'thistle']
    highlights = ['dodgerblue', 'purple']
    
    # Iterate through all the chromosomes and plot each one separately
    # While iterating, keep running lists of chromosome names and x-label positions as you encounter them
    x_labels = []
    x_labels_pos = []
    
    for i, (name, group) in enumerate(groups):
        ax.scatter(group.position, group.pvalue, s=1, color=colors[i%2])
        ax.scatter(group.position, group.sigs, s=1, color=highlights[i%2])
        x_labels.append(name)
        # Set the position of the x label to be the center of this chromosome's x-range
        x_labels_pos.append(np.median(group.position))
    
    ax.set_xticklabels(x_labels, rotation=50, fontsize=8)
    ax.set_xticks(x_labels_pos)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log(p-value)")
    ax.set_title(treatment)
    fig.savefig(treatment + "_manplot.png")
    plt.close()