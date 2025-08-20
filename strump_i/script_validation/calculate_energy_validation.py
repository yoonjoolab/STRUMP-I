#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 17:33:25 2021

@author: duaghk

# energy calculation validation.

"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from math import floor
from matplotlib import pyplot as plt
from string import ascii_uppercase as abt


# define functions.
def check_skiplines(fn):
    count = 0
    with open(fn) as f:
        while True:
            line = f.readline()
            if line.startswith("Pdb"):
                break
            count += 1
    return count


def parse_data(fn):
    skipnum = check_skiplines(fn)
    df = pd.read_csv(fn, sep="\t", skiprows=skipnum)
    df["Pdb"] = [x.split("/")[-2] for x in df["Pdb"]]
    df = df[df["Group1"].isin(["A"])]
    return df


def chain_pos_dict_gen(tgt_n: int):
    chain = abt[15:15+tgt_n]
    pos = list([x+1 for x in range(tgt_n)])
    returndict = {chain[i]: pos[i] for i in range(len(pos))}
    return returndict


def parse_data_with_chain(maindir, allele_list):
    returndf = pd.DataFrame()
    for allele in tqdm(allele_list):
        allele_dir = os.path.join(maindir, allele)
        ap_list = os.listdir(allele_dir)
        ap_list = [x for x in ap_list if allele in x]
        for ap in ap_list:
            ap_dir = os.path.join(allele_dir, ap)
            # find tgt files.
            tgt_fn = [x for x in os.listdir(ap_dir)
                      if x.startswith("Interaction") and "modified" in x]
            if len(tgt_fn):
                tgt_fn = tgt_fn[0]
            else:
                continue
            tmpdf = parse_data(os.path.join(ap_dir, tgt_fn))
            returndf = returndf.append(tmpdf)
    returndf = returndf.reset_index(drop=True)
    return returndf


def parse_data_with_peptide(maindir, allele_list):
    returndf = pd.DataFrame()
    for allele in tqdm(allele_list):
        allele_dir = os.path.join(maindir, allele)
        ap_list = os.listdir(allele_dir)
        ap_list = [x for x in ap_list if allele in x]
        for ap in ap_list:
            ap_dir = os.path.join(allele_dir, ap)
            # find tgt files.
            tgt_fn = [x for x in os.listdir(ap_dir)
                      if x.startswith("Interaction") and "Repair" in x]
            if len(tgt_fn):
                tgt_fn = tgt_fn[0]
            else:
                continue
            tmpdf = parse_data(os.path.join(ap_dir, tgt_fn))
            returndf = returndf.append(tmpdf)
    returndf = returndf.reset_index(drop=True)
    return returndf


def aggregate_data(df):
    returndf = pd.DataFrame()
    infocol = ["Pdb", "Group1", "Group2"]
    elsecol = df.columns.drop(infocol).tolist()
    for ap, df in tqdm(df.groupby("Pdb")):
        infoseries = df[infocol].iloc[1]
        df = df[elsecol]
        tmpseries = df.sum()
        tmpseries = infoseries.append(tmpseries)
        returndf = returndf.append(tmpseries, ignore_index=True)
    return returndf


maindir = "/Users/duaghk/data/strumpi/new_strumpi_output/repair_test/new"
allele_list = ["A0201"]

chain_wise_df = parse_data_with_chain(maindir, allele_list)
chain_aggregated = aggregate_data(chain_wise_df)
peptide_df = parse_data_with_peptide(maindir, allele_list)

infocol = ["Pdb", "Group1", "Group2"]
elsecol = chain_wise_df.columns.drop(infocol).tolist()

merged_df = pd.merge(peptide_df, chain_aggregated, on="Pdb")

fig, axes = plt.subplots(nrows=5, ncols=6, figsize=(20, 15))
for i in range(len(elsecol)):
    colname = elsecol[i]
    ncol = i % 6
    nrow = floor(i/6)
    axes[nrow][ncol].scatter(merged_df[f"{colname}_x"],
                             merged_df[f"{colname}_y"],
                             s=3, alpha=0.5)
    axes[nrow][ncol].set_xlabel("Peptide")
    axes[nrow][ncol].set_ylabel("Position-wise")
    axes[nrow][ncol].set_title(colname)
plt.tight_layout()
plt.show()

# %% analysis peptide chain.(in PDB)
# make 2d daataframe.
'''
rows = peptide
cols = chain(position wise)
'''

def parse_chain(tgt_fn):
    tmpdict = {}
    aadict = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
        }
    with open(tgt_fn) as f_read:
        lines = f_read.readlines()
    lines = [x for x in lines if x.startswith("ATOM") and x[21] != "A"]
    for line in lines:
        aa = aadict[line[17:20]]
        chain = line[21]
        tmpdict[chain] = aa
    returnseries = pd.Series(tmpdict)
    return returnseries


def parse_peptide_chain(maindir, allele_list):
    returndf = pd.DataFrame()
    for allele in tqdm(allele_list):
        allele_dir = os.path.join(maindir, allele)
        ap_list = os.listdir(allele_dir)
        ap_list = [x for x in ap_list if allele in x]
        for ap in ap_list:
            ap_dir = os.path.join(allele_dir, ap)
            # find tgt files.
            tgt_fn = "builded_template_minim_peptide_chain_modified.pdb"
            tmpseries = parse_chain(os.path.join(ap_dir, tgt_fn))
            peptides = ap.split("_")[1]
            tmpseries.name = peptides
            returndf = returndf.append(tmpseries)
    returndf = returndf
    return returndf


pep_seq_chain = parse_peptide_chain(maindir, allele_list)
pep_seq_chain["parsed_peptide"] = pep_seq_chain.sum(axis=1)
pep_seq_chain["check_equal"] = [row.name == row["parsed_peptide"] for i, row in pep_seq_chain.iterrows()]
sum(pep_seq_chain["check_equal"])



