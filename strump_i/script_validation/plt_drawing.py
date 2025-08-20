#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 13:24:32 2021

@author: duaghk
"""

# %% import library
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from math import floor
from random import sample
from matplotlib import pyplot as plt
from string import ascii_uppercase as abt


# define function
def cal_bins(tgt_series, num=100):
    bins = np.linspace(tgt_series.min(), tgt_series.max(), num)
    return bins


def chain_pos_dict_gen(tgt_n: int):
    chain = abt[15:15+tgt_n]
    pos = list([x+1 for x in range(tgt_n)])
    returndict = {chain[i]: pos[i] for i in range(len(pos))}
    return returndict


# %% set default values
maindir = "/Users/duaghk/data/strumpi/new_data/runned_output"
savedir = os.path.dirname(maindir)
allele_list = os.listdir(maindir)

# data load
chain_wise_df = pd.read_csv(os.path.join(savedir, "chain_wise_data_with_stability.csv"),
                            header=0)
chain_aggregated = pd.read_csv(os.path.join(savedir, "chain_wise_aggregated_data.csv"),
                               header=0)
peptide_df = pd.read_csv(os.path.join(savedir, "peptide_data_with_stability.csv"),
                         header=0)

chain_pos_dict = chain_pos_dict_gen(11)
chain_wise_df["Group2"] = chain_wise_df["Group2"].replace(chain_pos_dict)
chain_wise_df["pos_len"] = [len(x.split("_")[1]) for x in chain_wise_df["Pdb"]]

# allele, position plot.
entropy_feature = ['Interaction Energy', 'Backbone Hbond', 'Sidechain Hbond', 'entropy sidechain', 'entropy mainchain']
hydrophobic_feature = ['Interaction Energy', 'Van der Waals', 'Solvation Hydrophobic', 'Van der Waals clashes']
hydrophilic_feature = ['Interaction Energy', 'Solvation Polar', 'Electrostatics', 'electrostatic kon']
stability_feature = ['Interaction Energy', 'StabilityGroup1', 'StabilityGroup2']

feature_dict = {"entropy": entropy_feature, "hydrophobic": hydrophobic_feature,
                "hydrophilic": hydrophilic_feature, "stability": stability_feature}


for allele in tqdm(allele_list):
    allele_df = chain_wise_df[[True if allele in x else False for x in chain_wise_df["Pdb"]]]
    for p_len in allele_df["pos_len"].unique():
        for feat_n, feat_list in feature_dict.items():
            nrows = len(feat_list)
            ncols = p_len+1
            tmpdf = allele_df[allele_df["pos_len"].isin([p_len])]
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                     figsize=(ncols*4, nrows*3))
            for nrow in range(len(feat_list)):
                feat = feat_list[nrow]
                for ncol in range(p_len):
                    pos = ncol+1
                    tmpdf2 = tmpdf[tmpdf["Group2"].isin([pos])]
                    bins = cal_bins(tmpdf2[feat])
                    axes[nrow][ncol].hist(tmpdf2[feat][tmpdf2["Quality"].isin([0])],
                                          color="red", alpha=0.5, label="Non-binder",
                                          bins=bins, density=True)
                    axes[nrow][ncol].hist(tmpdf2[feat][tmpdf2["Quality"].isin([1])],
                                          color="blue", alpha=0.5, label="Binder",
                                          bins=bins, density=True)
                    if ncol == 0:
                        axes[nrow][ncol].set_ylabel(feat)
                    if nrow == 0:
                        axes[nrow][ncol].set_title(pos)
                    axes[nrow][ncol].legend()
                bins = cal_bins(tmpdf[feat])
                axes[nrow][ncols-1].hist(tmpdf[feat][tmpdf["Quality"].isin([0])],
                                         color="red", alpha=0.5, label="Non-binder",
                                         bins=bins, density=True)
                axes[nrow][ncols-1].hist(tmpdf[feat][tmpdf["Quality"].isin([1])],
                                         color="blue", alpha=0.5, label="Binder",
                                         bins=bins, density=True)
                if nrow == 0:
                    axes[nrow][ncols-1].set_title("Total")
            fig.suptitle(f"Length {p_len} {feat_n} histogram color by Answers", fontsize=15)
            plt.tight_layout()
            figname = f"{allele}_{p_len}_{feat_n}.png"
            plt.savefig(f"{os.path.join(savedir, 'figure', figname)}")
            plt.close()
