#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 10:26:12 2021

@author: duaghk

optimize template validation.

validation procedure
1. check minimize energy comparison.
1-1. compare input RMS vs minimized final RMS

2. check repairing
2.1 compare first energy and last energy.
2.2 last energy change check.
2.3 iteration number check.

"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from matplotlib import pyplot as plt

# define functions
def parse_optimize_log(fn):
    # get target lines.
    with open(fn) as f_read:
        lines = f_read.readlines()
    # get optimized RMS
    tgt_str = "Final RMS Gradient"
    tgt_line = [x for x in lines if tgt_str in x]
    final_grad = float(tgt_line[0].split(" ")[-1].strip())
    # get iter time.
    tgt_str = "Attempt"
    tgt_lines = [x for x in lines if tgt_str in x]
    iter_nums = [int(x.split(" ")[-1].replace(")", "")) for x in tgt_lines]
    iter_nums = max(iter_nums)
    # get repairing energy.
    tgt_str = "Improvement"
    tgt_lines = [x.split("INFO - ")[1].split(" (")[0] for x in lines if tgt_str in x]
    e_list = [float(x.split("-> ")[1]) for x in tgt_lines]
    first_e = e_list[0]
    final_e = e_list[-1]
    total_imp = round((first_e - final_e), 4)
    last_imp = round((e_list[-2] - final_e), 4)
    returndict = {}
    returndict["optimized_RMS"] = final_grad
    returndict["iteration"] = iter_nums
    returndict["primary_e"] = first_e
    returndict["final_e"] = final_e
    returndict["total_improvement"] = total_imp
    returndict["last_improvement"] = last_imp
    return returndict


def parse_optimize_log2(fn):
    # get target lines.
    with open(fn) as f_read:
        lines = f_read.readlines()
    # get optimized RMS
    tgt_str = "Final RMS Gradient"
    tgt_line = [x for x in lines if tgt_str in x]
    final_grad = float(tgt_line[0].split(" ")[-1].strip())
    # get iter time.
    tgt_str = "Attempt"
    tgt_lines = [x for x in lines if tgt_str in x]
    iter_nums = [int(x.split(" /")[0].split(" ")[-1].replace(")", "")) for x in tgt_lines]
    iter_nums = max(iter_nums)
    # get repairing energy.
    tgt_str = "Improvement"
    tgt_lines = [x.split("INFO - ")[1].split(" (")[0] for x in lines if tgt_str in x]
    e_list = [float(x.split("-> ")[1]) for x in tgt_lines]
    first_e = e_list[0]
    final_e = e_list[-1]
    total_imp = round((first_e - final_e), 4)
    last_imp = round((e_list[-2] - final_e), 4)
    returndict = {}
    returndict["optimized_RMS"] = final_grad
    returndict["iteration"] = iter_nums
    returndict["primary_e"] = first_e
    returndict["final_e"] = final_e
    returndict["total_improvement"] = total_imp
    returndict["last_improvement"] = last_imp
    return returndict


# set variables
maindir = "/Users/duaghk/data/strumpi/new_strumpi_output/sampled_full/"
allelelist = ["A0201"]
prev_optimized_df = pd.DataFrame()
for allele in allelelist:
    allele_dir = os.path.join(maindir, "previous", allele)
    ap_list = os.listdir(allele_dir)
    for ap in ap_list:
        tmp_fn = os.path.join(allele_dir, ap, f"{ap}.log")
        tmp_dict = parse_optimize_log(tmp_fn)
        tmp_series = pd.Series(tmp_dict, name=ap)
        prev_optimized_df = prev_optimized_df.append(tmp_series)

new_optimized_df = pd.DataFrame()
for allele in allelelist:
    allele_dir = os.path.join(maindir, "new", allele)
    ap_list = os.listdir(allele_dir)
    for ap in ap_list:
        tmp_fn = os.path.join(allele_dir, ap, f"{ap}.log")
        tmp_dict = parse_optimize_log2(tmp_fn)
        tmp_series = pd.Series(tmp_dict, name=ap)
        new_optimized_df = new_optimized_df.append(tmp_series)

merged = pd.merge(new_optimized_df, prev_optimized_df,
                  left_index=True, right_index=True)

# %% plot drawing

# plt.figure(figsize=(8, 6))
tgt_col = "optimized_RMS"
plt.scatter(merged[f"{tgt_col}_x"], merged[f"{tgt_col}_y"],
            s=5, c="blue", alpha=0.5)
plt.axhline(y=1, label="queried cutoff", color="gray")
plt.axvline(x=1, color="gray")
plt.xlabel("New script")
plt.ylabel("Previous script")
plt.legend()
plt.title("Optimized RMS")
plt.show()

# plt.figure(figsize=(8,6))
tgt_col = "iteration"
plt.scatter(merged[f"{tgt_col}_x"], merged[f"{tgt_col}_y"],
            s=5, c="blue", alpha=0.5)
plt.xlabel("New script")
plt.ylabel("Previous script")
plt.title("Iteration comparision")
plt.show()

tgt_col = "primary_e"
plt.scatter(merged[f"{tgt_col}_x"], merged[f"{tgt_col}_y"],
            s=5, c="blue", alpha=0.5)
plt.axline((merged[f"{tgt_col}_x"].min(), merged[f"{tgt_col}_x"].min()),
           (merged[f"{tgt_col}_x"].max(), merged[f"{tgt_col}_x"].max()), alpha=0.2, color="gray")
plt.xlabel("New script")
plt.ylabel("Previous script")
plt.title("Primary energy")
plt.show()

tgt_col = "final_e"
plt.scatter(merged[f"{tgt_col}_x"], merged[f"{tgt_col}_y"],
            s=5, c="blue", alpha=0.5)
plt.axline((merged[f"{tgt_col}_x"].min(), merged[f"{tgt_col}_x"].min()),
           (merged[f"{tgt_col}_x"].max(), merged[f"{tgt_col}_x"].max()), alpha=0.2, color="gray")
plt.xlabel("New script")
plt.ylabel("Previous script")
plt.title("Last energy")
plt.show()


tgt_col = "total_improvement"
plt.scatter(merged[f"{tgt_col}_x"], merged[f"{tgt_col}_y"],
            s=5, c="blue", alpha=0.5)
plt.axline((merged[f"{tgt_col}_x"].min(), merged[f"{tgt_col}_x"].min()),
           (merged[f"{tgt_col}_x"].max(), merged[f"{tgt_col}_x"].max()), alpha=0.2, color="gray")
plt.xlabel("New script")
plt.ylabel("Previous script")
plt.title("Energy improvement")
plt.show()

merged["energy_gap"] = merged["total_improvement_x"] - merged["total_improvement_y"]

plt.hist(merged["energy_gap"])

print(merged["energy_gap"] )


# %% 500 iteration results comparison
maindir = "/Users/duaghk/data/strumpi/new_strumpi_output/repair_test"
tgt_list = ["previous", "new"]

def parse_results(datadir):
    returndf = pd.DataFrame()
    tgt_fn = os.path.join(datadir, "A0201", "repair_results.csv")
    with open(tgt_fn) as f_read:
        lines = f_read.read().split("/data1")[1:]
        for line in tqdm(lines):
            tmpseries = pd.Series(line.split(","),
                                  index = ["pdb", "start_e", "final_e", "iter_n"])
            returndf = returndf.append(tmpseries, ignore_index=True)
    # pdb change.
    returndf["pdb"] = [os.path.basename(os.path.dirname(x)) for x in returndf["pdb"]]
    # get energy gap
    returndf["final_e"] = returndf["final_e"].astype(float)
    returndf["start_e"] = returndf["start_e"].astype(float)
    returndf["gap_e"] = returndf["final_e"] - returndf["start_e"]
    return returndf

prev_df = parse_results(os.path.join(maindir, tgt_list[0]))
new_df = parse_results(os.path.join(maindir, tgt_list[1]))

merged = pd.merge(prev_df, new_df, on="pdb")
merged["improved"] = [True if row["gap_e_x"] > row["gap_e_y"] else False for i, row in merged.iterrows()]

tgt_col = "gap_e"
plt.figure(figsize=(8,6))
plt.axline((merged[f"{tgt_col}_x"].min(), merged[f"{tgt_col}_x"].min()), 
           (merged[f"{tgt_col}_x"].max(), merged[f"{tgt_col}_x"].max()), 
           color="gray", alpha=0.5)
plt.scatter(merged[f"{tgt_col}_x"], merged[f"{tgt_col}_y"],
            s=5, c="blue", alpha=0.5, label="not_imporved")
plt.xlabel("Previous")
plt.ylabel("New")
plt.title("Energy gap")
plt.show()

tgt_col = "start_e"
plt.figure(figsize=(8,6))
plt.axline((merged[f"{tgt_col}_x"].min(), merged[f"{tgt_col}_x"].min()), 
           (merged[f"{tgt_col}_x"].max(), merged[f"{tgt_col}_x"].max()), 
           color="gray", alpha=0.5)
plt.scatter(merged[f"{tgt_col}_x"], merged[f"{tgt_col}_y"],
            s=5, c="blue", alpha=0.5, label="not_imporved")
plt.xlabel("Previous")
plt.ylabel("New")
plt.title("Start energy")
plt.show()

