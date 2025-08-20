#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 21:11:04 2021

@author: duaghk

# vaccc data generation.

"""

import os
import pandas as pd

maindir = "/Users/duaghk/data/strumpi/vaccc_data/raw_data"

os.listdir(maindir)

data1 = os.path.join(maindir, "VACCC01_strumpi_input.csv")
data5 = os.path.join(maindir, "VACCC05_strumpi_input.csv")

def make_input(df, outdir):
    for allele, df in df.groupby("allele"):
        peplist = df["peptide"].tolist()
        with open(os.path.join(outdir, f"{allele}.csv"), 'w') as f:
            f.write("\n".join(peplist))
    return

df1 = pd.read_csv(data1)
make_input(df1, "/Users/duaghk/data/strumpi/vaccc_data/input/vaccc-01")

df5 = pd.read_csv(data5)
make_input(df5, "/Users/duaghk/data/strumpi/vaccc_data/input/vaccc-05")


