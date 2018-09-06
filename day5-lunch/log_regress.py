#!/usr/bin/env python3

"""
Usage: ./log_regress.py <tab_file>

<ctab_file> Filepath to .tab file with transcript IDs,
and their associated histone modifications (by bigWig) and FPKMs

This script takes a file with gene transcripts, their FPKMs, and bigWigAverageOverBed's
predicted histone modifications over the promoter region (for as many types of histone
modification as you like). It fits a multiple linear regression model to the data, 
trying to use histone modifications to predict a transcript's FPKM.

To make the data nicer, here we use log(FPKM + 1) instead of FPKM to fit the linear 
regression model. This results in some strange artifacts.

The script then compares the model's projected FPKM values to the actual FPKMs in the 
.tab file, and plots the residuals in histogram format.
"""

import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from sklearn import datasets

# Read .tab file
filename = pd.read_table(sys.argv[1])


# Fit linear regression model to .tab data
# Y: Log of FPKM values + 1 (to account for zeroes)
# X: Corresponding avg signal for each of 5 histone markers
Y = np.log(1 + filename.loc[:,"FPKM"])
X = filename.iloc[:,2:]
# Not entirely sure why the constant needs to be added
X = sm.add_constant(X)
# Fit to multiple linear regression model
model = sm.OLS(Y,X)
results = model.fit()
results.params
results.tvalues
# Print summary to output R^2, p-value, coefficients, etc
print(results.summary())


# Make histogram

# new_y - FPKM values predicted by the model for each transcript
# In the histogram, we want to plot the difference between the predicted FPKM values
# and the actual values from the .tab file (i.e. the 'residuals')
new_y = results.fittedvalues
r = new_y - Y

# For the purpose of aesthetic, the xlim is set smaller than this histogram's default range
fig, ax = plt.subplots(figsize=(8, 6))
ax.hist(r, bins=500)
ax.set_ylabel("Frequency")
ax.set_xlabel("Residuals")
ax.set_title("Linear regression residuals")
fig.savefig("log_residuals.png")
plt.close(fig)