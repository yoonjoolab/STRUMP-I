#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 15:23:24 2021

@author: duaghk
"""

import os
import pandas as pd
from tqdm import tqdm
global aa_dict
aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
           'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def parse_pdb_to_df(templatefile):
    returndf = pd.DataFrame()
    lines = open(templatefile).readlines()
    seriesname = ["linestart", "aa", "line_middle", "chain", "pos", "alt", "line_end", "group_key"]
    for line in lines:
        if len(line) < 50:
            continue
        line_start = line[0:17]
        aa = line[17:20]
        line_middle = line[20]
        chain = line[21]
        pos = int(line[22:26].strip())
        alt = line[26]
        groupkey = line[22:27]
        line_end = line[27:]
        tmpseries = pd.Series(data=[line_start, aa, line_middle, chain,
                                    pos, alt, line_end, groupkey],
                              index=seriesname)
        returndf = returndf.append(tmpseries, ignore_index=True)
    returndf = returndf[seriesname]
    return returndf


def parse_indiv(indiv_fn):
    returndict = {}
    with open(indiv_fn) as f:
        line = f.readline()
    tokens = line.strip().replace(";", "").split(",")
    for token in tokens:
        chain = token[1]
        pos = token[2:-1]
        altered = token[-1]
        if chain not in returndict.keys():
            returndict[chain] = {}
        returndict[chain][pos] = altered
    return returndict


def check_matched(pdbdf, indiv_dict):
    for i, row in pdbdf.iterrows():
        chain = row["chain"]
        pos = int(row["pos"])
        aa = aa_dict[row["aa"]]
        if pos in indiv_dict[chain].keys():
            yield aa != indiv_dict[chain][str(pos)]


maindir = "/Users/duaghk/data/strumpi/ref_template_change/A0101_test/A0101"
# ap = "A0101_ASSLPTTMNY"


build_check_dict = {}
for ap in os.listdir(maindir):
    ap_dir = os.path.join(maindir, ap)
    indiv_fn = os.path.join(ap_dir, "individual_list.txt")
    pdb_fn = os.path.join(ap_dir, "builded_template_minim_best_Repair.pdb")
    pdbdf = parse_pdb_to_df(pdb_fn)
    indiv_dict = parse_indiv(indiv_fn)
    build_check_dict[ap] = sum(check_matched(pdbdf, indiv_dict))

sum(build_check_dict.values())
row = pdbdf.iloc[343]



