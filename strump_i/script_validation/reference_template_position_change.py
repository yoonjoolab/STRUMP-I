#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 23:18:36 2021

@author: duaghk

STRUMP-I template change
alter AA position change.

"""

import os
import pandas as pd
from tqdm import tqdm


def parse_pdb_to_df(templatefile):
    returndf = pd.DataFrame()
    lines = open(templatefile).readlines()
    seriesname = ["linestart", "aa", "line_middle", "chain", "pos", "alt", "line_end", "group_key"]
    for line in lines:
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


def reorder_pdb_position(pdbdf):
    pos_count = 0
    hla_chain = pdbdf[pdbdf["chain"].isin(["A"])]
    pep_chain = pdbdf[pdbdf["chain"].isin(["P"])]
    changed_hla = []
    for key, df in hla_chain.groupby("group_key"):
        pos_count += 1
        df["pos"] = pos_count
        changed_hla.append(df)
    changed_hla = pd.concat(changed_hla)

    changed_pdb = pd.concat([changed_hla, pep_chain])
    changed_pdb = changed_pdb.drop(columns="group_key")
    changed_pdb["pos"] = [str(int(x)) for x in changed_pdb["pos"]]
    changed_pdb["pos"] = [x.rjust(4, " ") for x in changed_pdb["pos"]]
    changed_pdb["pos"] = [x.ljust(5, " ") for x in changed_pdb["pos"]]
    changed_pdb = changed_pdb.drop(columns="alt")
    return changed_pdb


def process_ref(template_fn, outdir):
    base_n = os.path.basename(template_fn)
    write_pdb = os.path.join(outdir, base_n)
    pdbdf = parse_pdb_to_df(template_fn)
    changed = reorder_pdb_position(pdbdf)
    with open(write_pdb, 'w') as f_write:
        for i in changed.index:
            line = "".join(list(map(str, changed.loc[i])))
            f_write.write(line)


pdbdir = "./reference/renumbered_pdb"
pdbdir2 = "./reference/renumbered_pdb2"
os.mkdir(pdbdir2)

for pdb in tqdm(os.listdir(pdbdir)):
    process_ref(os.path.join(pdbdir, pdb), pdbdir2)


def check_ordering(pdb_fn):
    pdbdf = parse_pdb_to_df(pdb_fn)
    start_pos = 1
    for pos in pdbdf["pos"][pdbdf["chain"].isin(["A"])]:
        if int(pos)-start_pos > 2 or int(pos)-start_pos < 0:
            print(pos)
        start_pos = pos






    