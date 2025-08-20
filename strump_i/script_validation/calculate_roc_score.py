#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 08:28:49 2021

@author: duaghk

Calculate ROC score

"""

# %% import library
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot as plt
from seaborn import clustermap
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


def make_len_seq_dict(p_len):
    lendict = list(range(1, p_len+1))
    forward = lendict[:4]
    forward_dict = {i: i for i in forward}
    backward = lendict[-4:]
    backlist = list(range(-4, 0, 1))
    backward_dict = {backward[i]: i for i in backlist}
    forward_dict.update(backward_dict)
    return forward_dict


def df_sampling_and_iteration(sampled, featurelist, iter_n):
    returndf = pd.DataFrame()
    for pos in sampled["Group2"].unique():
        sampled_pos = sampled[sampled["Group2"].isin([pos])]
        tmpdict = {x: roc_auc_score(sampled_pos["Quality"],
                                    sampled_pos[x]) for x in featurelist}
        tmpdict.update({"pos": pos, "iter_n": iter_n})
        tmpseries = pd.Series(tmpdict)
        returndf = returndf.append(tmpseries, ignore_index=True)
    return returndf


def iter_roc(df, n_sampling, featurelist):
    returndf = []
    iter_df = [
        df_sampling_and_iteration(df.groupby(["Group2",
                                              "Quality"]).sample(n=n_sampling,
                                                                 replace=True),
                                  featurelist,
                                  i) for i in range(30)]
    iter_df = pd.concat(iter_df).reset_index(drop=True)
    tmpdf = iter_df.groupby("pos").mean().drop(columns="iter_n")
    for pos in tmpdf.index:
        tmpdf2 = pd.DataFrame(tmpdf.loc[pos])
        tmpdf2.columns = ["roc_score"]
        tmpdf2["feature"] = tmpdf2.index
        tmpdf2["pos"] = pos
        returndf.append(tmpdf2)
    returndf = pd.concat(returndf)
    return returndf


def parsing_df_to_cal_auroc(df, featurelist):
    # allele, p_len, feature, pos, roc(iter 30 average), n_sampling
    allele = [x.split("_")[0] for x in df["Pdb"]][0]
    p_len = df["pos_len"].unique().tolist()[0]
    # check df Quality.
    qual_count = df["Quality"].value_counts()/p_len
    qual_count = qual_count.astype(int)
    n_sampling = qual_count.min()
    if len(qual_count) == 1 or n_sampling < 30:
        returndf = None
    else:
        pos_dict = make_len_seq_dict(p_len)
        df = df[df["Group2"].isin(pos_dict.keys())]
        df["Group2"] = df["Group2"].replace(pos_dict)
        returndf = iter_roc(df, n_sampling, featurelist)
        returndf["allele"] = allele
        returndf["p_len"] = p_len
        returndf["n_sampling"] = n_sampling
        returndf = returndf.reset_index(drop=True)
        # column reordering
        returndf = returndf[["allele", "pos", "p_len",
                             "feature", "n_sampling", "roc_score"]]
    return returndf


# set values
# %% set default values
maindir = "/Users/duaghk/data/strumpi/new_data/runned_output"
savedir = os.path.dirname(maindir)
# allele_list = os.listdir(maindir)

# data load
chain_wise_df = pd.read_csv(os.path.join(savedir, "chain_wise_data_with_stability_211006.csv"),
                            header=0)
chain_aggregated = pd.read_csv(os.path.join(savedir, "chain_wise_aggregated_data_211006.csv"),
                               header=0)
peptide_df = pd.read_csv(os.path.join(savedir, "peptide_data_with_stability_211006.csv"),
                         header=0)
allele_list = [x.split("_")[0] for x in peptide_df["Pdb"]]
allele_list = list(set(allele_list))
# data filtering.
wrong_pdb = peptide_df[peptide_df["StabilityGroup1"] > 100]
chain_wise_df = chain_wise_df[~chain_wise_df["Pdb"].isin(wrong_pdb["Pdb"].tolist())]

# Group2 (peptide chain) change to position number
chain_pos_dict = chain_pos_dict_gen(11)
chain_wise_df["Group2"] = chain_wise_df["Group2"].replace(chain_pos_dict)
chain_wise_df["pos_len"] = [len(x.split("_")[1]) for x in chain_wise_df["Pdb"]]

features = [
    'Interaction Energy', 'Backbone Hbond', 'Sidechain Hbond',
    'entropy sidechain', 'entropy mainchain', 'Interaction Energy',
    'Van der Waals', 'Solvation Hydrophobic', 'Van der Waals clashes',
    'Interaction Energy', 'Solvation Polar', 'Electrostatics',
    'electrostatic kon', 'Interaction Energy', 'StabilityGroup1',
    'StabilityGroup2'
    ]


# # code test
# tmpdf = chain_wise_df[[True if "A2402" in x else False for x in chain_wise_df["Pdb"]]]
# tmpdf = tmpdf[tmpdf["pos_len"].isin([10])]
# df = tmpdf

roc_df = []
for allele in tqdm(allele_list):
    tmpdf = chain_wise_df[[True if allele in x else False for x in chain_wise_df["Pdb"]]]
    for pos in tmpdf["pos_len"].unique().tolist():
        tmpdf2 = tmpdf[tmpdf["pos_len"].isin([pos])]
        returns = parsing_df_to_cal_auroc(tmpdf2, features)
        roc_df.append(returns)

roc_df = pd.concat(roc_df, ignore_index=True)
# reversed ROC change
roc_df["roc_score"] = roc_df["roc_score"]-0.5
roc_df["roc_score"] = abs(roc_df["roc_score"])

roc_df_back = roc_df.copy()
roc_df = roc_df[roc_df["p_len"].isin([9,10])]
# heatmap
heatmap_df = pd.DataFrame()
for allele, tmpdf in roc_df.groupby("allele"):
    tmpdf["colname"] = tmpdf["p_len"].astype(str) + "_" + tmpdf["pos"].astype(int).astype(str) + "_" + tmpdf["feature"]
    tmpseries = pd.Series(data = tmpdf["roc_score"].tolist(), index = tmpdf["colname"], name = allele)
    heatmap_df = heatmap_df.append(tmpseries)

# fill na
heatmap_df = heatmap_df.fillna(0)
plt.figure(figsize=(45,30), dpi = 300)
plt.pcolor(heatmap_df)
plt.yticks(ticks = np.arange(heatmap_df.shape[0])+0.5, labels = heatmap_df.index)
plt.xticks(ticks = np.arange(heatmap_df.shape[1])+0.5, labels = heatmap_df.columns, rotation=90)
plt.colorbar()
plt.show()

clustermap(heatmap_df, figsize=(45, 30), cmap='mako')

heatmap_9 = heatmap_df[[x for x in heatmap_df.columns if x.startswith("9")]]
heatmap_10 = heatmap_df[[x for x in heatmap_df.columns if x.startswith("10")]]

clustermap(heatmap_9, figsize=(30, 20), cmap='mako')
clustermap(heatmap_10, figsize=(30, 20), cmap='mako')


roc_df.to_csv("/Users/duaghk/data/strumpi/new_data/allele_pep_feature_wise_ROC.csv", index=False)


roc_df_filtered = roc_df[roc_df["roc_score"] >= 0.15]
d = roc_df_filtered.groupby(["p_len", "pos", "feature"]).count()

roc_df_filtered["feature"].value_counts()
roc_df_filtered["colname"] = roc_df_filtered["p_len"].astype(str) + "_" + roc_df_filtered["pos"].astype(int).astype(str) + "_" + roc_df_filtered["feature"]
d = roc_df_filtered["colname"].value_counts()
len(roc_df_filtered["allele"].unique())

def draw_hist(df, tgtcol, title):
    bins = cal_bins(df[tgtcol], num=int(len(df[tgtcol])/5))
    plt.hist(tmpdf[tgtcol][tmpdf["Quality"].isin([0])], bins=bins, density=True,
             color="red", alpha=0.5, label="Non-binder")
    plt.hist(tmpdf[tgtcol][tmpdf["Quality"].isin([1])], bins=bins, density=True,
             color="blue", alpha=0.5, label="Binder")
    plt.legend()
    plt.title(title)
    plt.show()

tmpdf = chain_wise_df[[True if "A3001" in x else False for x in chain_wise_df["Pdb"]] & chain_wise_df["pos_len"].isin([9]) & chain_wise_df["Group2"].isin([3])]
draw_hist(tmpdf, "entropy sidechain", "A3001, length 9, position 3\nentropy sidechain")

tmpdf = chain_wise_df[[True if "A3001" in x else False for x in chain_wise_df["Pdb"]] & chain_wise_df["pos_len"].isin([9]) & chain_wise_df["Group2"].isin([2])]
draw_hist(tmpdf, "electrostatic kon", "A3001, length 9, position 2\nelectrostatic kon")

tmpdf = chain_wise_df[[True if "A3001" in x else False for x in chain_wise_df["Pdb"]] & chain_wise_df["pos_len"].isin([9]) & chain_wise_df["Group2"].isin([2])]
draw_hist(tmpdf, "Solvation Polar", "A3001, length 9, position -2\nSolvation Polar")

tmpdf = chain_wise_df[[True if "A3001" in x else False for x in chain_wise_df["Pdb"]] & chain_wise_df["pos_len"].isin([9]) & chain_wise_df["Group2"].isin([1])]
draw_hist(tmpdf, "Van der Waals", "A3001, length 9, position 1\nVan der Waals")





