#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 11:27:04 2021

@author: duaghk

build_template validation

check input peptide in peptide chain.

"""

# import library
import os
import pandas as pd
from tqdm import tqdm

# set global variable.
global aa_dict
aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
           'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


# define functions
def get_pdb_peptide(pdb_fn):
    with open(pdb_fn) as f:
        lines = [x.strip() for x in f.readlines() if len(x) > 21 and x[21] == "P"]
    pep_aa_dict = {}
    for line in lines:
        aa = aa_dict[line[17:20]]
        pos = line[22:26].strip()
        pep_aa_dict[pos] = aa
    tmpseries = pd.Series(pep_aa_dict)
    tmpseries = tmpseries.sort_index()
    returnpep = "".join(tmpseries.tolist())
    return returnpep


def parse_all(maindir, allelelist):
    returndf = pd.DataFrame()
    for allele in tqdm(allelelist):
        allele_dir = os.path.join(maindir, allele)
        ap_list = os.listdir(allele_dir)
        for ap in ap_list:
            allele, peptide = ap.split("_")
            ap_dir = os.path.join(allele_dir, ap)
            tgtfile = os.path.join(ap_dir, "builded_template.pdb")
            builded_peptide = get_pdb_peptide(tgtfile)
            tmpseries = pd.Series([allele, peptide, builded_peptide],
                                  index=["allele", "queried_peptide", "builded_peptide"])
            returndf = returndf.append(tmpseries, ignore_index=True)
    returndf["queried_matched_builded"] = [row["queried_peptide"] == row["builded_peptide"]
                                           for i, row in returndf.iterrows()]
    return returndf


# set variable.
maindir = "/Users/duaghk/data/strumpi/new_strumpi_output/output/"
allelelist = ["A2402", "B5101"]
matched_df = parse_all(maindir, allelelist)
sum(matched_df["queried_matched_builded"])  # 936.
for allele in allelelist:
    tmpdf = matched_df[matched_df["allele"].isin([allele])]
    printline = (
        f"{allele} number of queried peptide: {len(tmpdf)}, "
        f"number of queried peptide == builded peptide: {sum(tmpdf['queried_matched_builded'])}"
        )
    print(printline)

# %% code test.
tmpdict = {"2":"Q", "4":"S", "1":"P", "3":"R"}

tmpdict
tmpseries = pd.Series(tmpdict)
tmpseries = tmpseries.sort_index()
tmpseries.tolist()


