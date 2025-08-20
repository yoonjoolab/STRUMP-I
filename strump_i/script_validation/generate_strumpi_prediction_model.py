#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 18:08:11 2021

@author: duaghk

modeling script generation.

"""

# %% import library
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot as plt
import seaborn as sns
from string import ascii_uppercase as abt
from sklearn.manifold import TSNE
from functions import strumpi_scaler as ss
from sklearn.ensemble import RandomForestClassifier as rfc
from matplotlib import pyplot as plt

# define function
def chain_pos_dict_gen(tgt_n: int):
    chain = abt[15:15+tgt_n]
    pos = list([x+1 for x in range(tgt_n)])
    returndict = {chain[i]: pos[i] for i in range(len(pos))}
    return returndict


# %% data load
maindir = "/Users/duaghk/data/strumpi/new_data/runned_output"
savedir = os.path.dirname(maindir)
# allele_list = os.listdir(maindir)

# data load
peptide_df = pd.read_csv(os.path.join(savedir, "peptide_data_with_stability_211006.csv"),
                         header=0)
peptide_df["Allele"] = [x.split("_")[0] for x in peptide_df["Pdb"]]
peptide_df["Peptide"] = [x.split("_")[1] for x in peptide_df["Pdb"]]
peptide_df["p_len"] = [len(x) for x in peptide_df["Peptide"]]

# column reordering.
infocol = ['Pdb', 'Allele', 'Peptide', 'p_len', 'Group1', 'Group2', 'Quality']
elsecol = [x for x in peptide_df.columns if x not in infocol]
peptide_df = peptide_df[infocol + elsecol]

# filtering data.
over_pep = peptide_df[peptide_df["StabilityGroup1"] > 100]
peptide_df = peptide_df[~peptide_df["Pdb"].isin(over_pep["Pdb"].tolist())]

# set features
features = [
    'Interaction Energy',
    'Backbone Hbond', 'Sidechain Hbond', 'entropy sidechain', 'entropy mainchain',
    'Van der Waals', 'Solvation Hydrophobic', 'Van der Waals clashes',
    'Solvation Polar', 'Electrostatics', 'electrostatic kon',
    'StabilityGroup1', 'StabilityGroup2'
    ]

# %% RF model generation.
# Quality only train test split.
# data split.

datalist = {}
for i in range(30):
    print(i)
    train = peptide_df.groupby("Quality").sample(frac=0.85)
    test = peptide_df[~peptide_df["Pdb"].isin(train["Pdb"].tolist())]
    if len(set(train["Pdb"]).intersection(set(test["Pdb"]))):
        continue
    else:
        scaler = ss.StrumpiScaler(len(infocol))
        train = scaler.fit_scaler(train)
        test = scaler.scaling(test)
        model = rfc()
        model.fit(train[features], train["Quality"])

        # predict
        test_pred = model.predict_proba(test[features])
        test_pred = pd.DataFrame(test_pred)
        test_pred = test_pred[1].to_list()
        test_returndf = pd.DataFrame(data={"Pred": test_pred,
                                           "Pdb": test["Pdb"].tolist(),
                                           "Quality": test["Quality"].tolist()})
        roc = roc_auc_score(test_returndf["Quality"], test_returndf["Pred"])
        tmpdict = {"test": test, "train": train, "scaler": scaler,
                   "model": model, "roc": roc}
        datalist[i] = tmpdict

# import pickle
# with open(os.path.join(maindir, "iter_30_all_pep_rfc_model.pickle.gz"), 'wb') as f:
#     pickle.dump(datalist, f)

# with open(os.path.join(maindir, "iter_30_all_pep_rfc_model.pickle.gz"), 'rb') as f:
#     datalist = pickle.load(f)

roclist = [v["roc"] for _,v in datalist.items()]
sum(roclist)/len(roclist)

plt.figure(figsize=(3, 4))
sns.violinplot(roclist)
plt.ylim(0.5,1)
plt.xlabel("Total")
plt.tick_params(bottom=False, labelbottom=False)
plt.title("All peptide model ")
plt.show()

# %% allele, quality, p_len stratified total model generation.
maindir = "/Users/duaghk/data/strumpi/new_data/runned_output"
savedir = os.path.dirname(maindir)
# allele_list = os.listdir(maindir)

# data load
peptide_df = pd.read_csv(os.path.join(savedir, "peptide_data_with_stability_211006.csv"),
                         header=0)
peptide_df["Allele"] = [x.split("_")[0] for x in peptide_df["Pdb"]]
peptide_df["Peptide"] = [x.split("_")[1] for x in peptide_df["Pdb"]]
peptide_df["p_len"] = [len(x) for x in peptide_df["Peptide"]]

# column reordering.
infocol = ['Pdb', 'Allele', 'Peptide', 'p_len', 'Group1', 'Group2', 'Quality']
elsecol = [x for x in peptide_df.columns if x not in infocol]
peptide_df = peptide_df[infocol + elsecol]

# filtering data.
# weird data.
over_pep = peptide_df[peptide_df["StabilityGroup1"] > 100]
peptide_df = peptide_df[~peptide_df["Pdb"].isin(over_pep["Pdb"].tolist())]
# peptide length only 9, 10
peptide_df = peptide_df[peptide_df["p_len"].isin([9, 10])]

# set features
features = [
    'Interaction Energy',
    'Backbone Hbond', 'Sidechain Hbond', 'entropy sidechain', 'entropy mainchain',
    'Van der Waals', 'Solvation Hydrophobic', 'Van der Waals clashes',
    'Solvation Polar', 'Electrostatics', 'electrostatic kon',
    'StabilityGroup1', 'StabilityGroup2'
    ]

iterdict = {}
for i in range(30):
    train = peptide_df.groupby(["Allele", 'p_len', 'Quality']).sample(frac=0.85)
    # check_fraction = train.groupby(["Allele", 'p_len', 'Quality']).count()
    test = peptide_df[~peptide_df['Pdb'].isin(train["Pdb"].tolist())]
    # sum(set(train["Pdb"]).intersection(test["Pdb"]))

    # model generation
    scaler = ss.StrumpiScaler(len(infocol))
    train = scaler.fit_scaler(train)
    test = scaler.scaling(test)
    model = rfc()
    model.fit(train[features], train["Quality"])
    
    pred = model.predict_proba(test[features])
    pred = pd.DataFrame(pred)
    pred = pred[1].to_list()
    pred = pd.DataFrame(data={"Pred": pred,
                              "Pdb": test["Pdb"].tolist(),
                              "Quality": test["Quality"].tolist()})
    # roc_auc_score(pred["Quality"], pred["Pred"])
    preddict = {}
    for ap, df in test.groupby(["Allele", 'p_len']):
        allele, p_len = ap
        if len(df["Quality"].unique()) == 1:
            continue
        pred = model.predict_proba(df[features])
        pred = pd.DataFrame(pred)
        pred = pred[1].to_list()
        pred = pd.DataFrame(data={"Pred": pred,
                                  "Pdb": df["Pdb"].tolist(),
                                  "Quality": df["Quality"].tolist()})
        roc = roc_auc_score(pred["Quality"], pred["Pred"])
        preddict[f"{allele}_{p_len}mer"] = {"roc": roc,
                                            "pos_ratio": sum(df["Quality"].isin([1]))/len(df),
                                            "neg_ratio": sum(df["Quality"].isin([0]))/len(df)}
    
    allele_plen_wise_roc_df = pd.DataFrame.from_dict(preddict, orient='index')
    iterdict[i] = {"train": train, "test": test, "model": model,
                   "pred": allele_plen_wise_roc_df}

iterdict[1]["pred"]

rocdf = pd.DataFrame()
for i, v in iterdict.items():
    preddf = v["pred"]
    preddf.columns = [f"{x}_{i}" for x in preddf.columns]
    tgtcol = f"roc_{i}"
    if len(rocdf) == 0:
        rocdf = preddf[tgtcol]
    else:
        rocdf = pd.merge(rocdf, preddf[tgtcol],
                         left_index=True, right_index=True)

rocdf.to_csv()

rocdf2 = rocdf.transpose()
rocdf2.boxplot()
len(df)



# roclist = []
# for i, valdict in datalist.items():
#     roclist.append(valdict["roc"])
# max(roclist)
# roclist

# # %% TESLA data test.
# maindir = "/Users/duaghk/data/strumpi/tesla_data/tesla_output"
# savedir = os.path.dirname(maindir)
# peptide_df = pd.read_csv(os.path.join(savedir, "TESLA_peptide_data_with_stability.csv"),
#                          header=0)

# tesla_data = pd.read_excel("/Users/duaghk/data/strumpi/tesla_data/TESLA_consortium/Supplementary/1-s2.0-S0092867420311569-mmc4.xlsx")
# tesla_data = tesla_data[tesla_data["PEP_LEN"] < 12]
# tesla_data["MHC"] = tesla_data["MHC"].apply(lambda x: x.replace(":", "").replace("*", ""))
# tesla_data["ap"] = tesla_data["MHC"] + "_" + tesla_data["ALT_EPI_SEQ"]
# tesla_data = tesla_data.drop_duplicates(subset=["PMHC"])
# tesla_data["Quality"] = tesla_data["MEASURED_BINDING_AFFINITY"].apply(lambda x: 1 if x <= 500 else 0)
# # set model.
# model = datalist[0]["model"]
# scaler = datalist[0]['scaler']
# tesla_scaled = scaler.scaling(peptide_df)

# tesla_pred = model.predict_proba(tesla_scaled[features])
# tesla_pred = pd.DataFrame(tesla_pred)
# tesla_pred = tesla_pred[1].to_list()
# tesla_pred = pd.DataFrame(data={"Pred": tesla_pred,
#                                 "Pdb": tesla_scaled["Pdb"].tolist()})

# tesla_pred = pd.merge(tesla_pred, tesla_data[["ap", "Quality"]],
#                       left_on="Pdb", right_on="ap")

# roc_auc_score(tesla_pred["Quality"], tesla_pred["Pred"])



