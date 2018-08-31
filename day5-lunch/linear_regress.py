#!/usr/bin/env python3

"""
Usage: ./day5-exercise-4.py <ctab_file>

<ctab_file> Filepath to ctab file

"""

import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from sklearn import datasets

filename = pd.read_table(sys.argv[1])

Y = filename.loc[:,"FPKM"]
X = filename.iloc[:,2:]
X = sm.add_constant(X)
model = sm.OLS(Y,X)
results = model.fit()
results.params
results.tvalues

new_y = results.fittedvalues
r = new_y - Y

#HISTOGRAM

fig, ax = plt.subplots(figsize=(8, 6))
ax.hist(r, bins=1000)
ax.set_xlim(-1000,1000)
ax.set_ylabel("Frequency")
ax.set_xlabel("Residuals")
ax.set_title("Linear regression residuals")
fig.savefig("multiple_regression.png")
plt.close(fig)