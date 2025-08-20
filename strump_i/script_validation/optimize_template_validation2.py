#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 11:30:25 2021

@author: duaghk

# validation energy calculation

"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from matplotlib import pyplot as plt
from string import ascii_uppercase as abt
from sklearn.metrics import roc_auc_score


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


def aggregate_data(maindir, allele_list):
    returndf = pd.DataFrame()
    for allele in tqdm(allele_list):
        allele_dir = os.path.join(maindir, allele)
        ap_list = os.listdir(allele_dir)
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


def aggregate_data2(maindir, allele_list):
    returndf = pd.DataFrame()
    for allele in tqdm(allele_list):
        allele_dir = os.path.join(maindir, allele)
        ap_list = os.listdir(allele_dir)
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


def chain_pos_dict_gen(tgt_n: int):
    chain = abt[15:15+tgt_n]
    pos = list([x+1 for x in range(tgt_n)])
    returndict = {chain[i]: pos[i] for i in range(len(pos))}
    return returndict


def get_minimized_log(tgtdir, allele_list):
    returndf = pd.DataFrame()
    for allele in allele_list:
        allele_dir = os.path.join(tgtdir, allele)
        ap_list = os.listdir(allele_dir)
        ap_list = [x for x in ap_list if allele in x]
        for ap in tqdm(ap_list):
            ap_dir = os.path.join(allele_dir, ap)
            log_fn = os.path.join(ap_dir, f"{ap}.log")
            with open(log_fn) as f_read:
                lines = f_read.readlines()
            tgt_line = [x for x in lines if "Final RMS Gradient" in x][0]
            # print(tgt_line)
            final_rms = tgt_line.strip().split(" ")[-1]
            tmpseries = pd.Series([ap, final_rms],
                                  index=["ap", "optimized_RMS"])
            returndf = returndf.append(tmpseries, ignore_index=True)
    return returndf


# %% minimized check
tgt_dir = "/Users/duaghk/data/strumpi/new_strumpi_output/repair_test/new/"
allele_list = ["A0201"]
minimized = get_minimized_log(tgt_dir, allele_list)
minimized["optimized_RMS"] = minimized["optimized_RMS"].astype("float")
# histogram.
plt.axvline(x=1, label="Queried RMS", color="gray", alpha=0.5)
plt.hist(minimized["optimized_RMS"],
         bins = np.linspace(minimized["optimized_RMS"].min(),
                            minimized["optimized_RMS"].max(),
                            50),
         density=True, alpha=0.5)
plt.legend()
plt.xlabel("Final RMS")
plt.title("Optimized RMS")
plt.show()

# %% repaire test
maindir = "/Users/duaghk/data/strumpi/new_strumpi_output/output/"
allele_list = ["A0201", "A2402", "B5101"]
full_inter = aggregate_data(maindir, allele_list)
chain_pos_dict = chain_pos_dict_gen(9)
ap_inter = aggregate_data2(maindir, allele_list)

