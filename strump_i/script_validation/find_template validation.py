#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 10:42:17 2021

@author: duaghk

Benchmark template finder.

"""

import os
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
from math import floor
from random import sample
from matplotlib import pyplot as plt
from functions import find_template


# %% define functions
def make_logger(outdir, name=None, consoleset=True):
    '''
        logger making.
    '''
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    loggerformat = "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(loggerformat)
    if consoleset:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(formatter)
        logger.addHandler(console)
    else:
        pass
    loggerfile = os.path.join(outdir, name)
    file_handler = logging.FileHandler(filename=f"{loggerfile}.log")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger


def parse_tf_result(tf_results_fn):
    tf_results = pd.DataFrame()
    with open(tf_results_fn) as f:
        for line in f.readlines():
            tokens = line.strip().split(" - ")
            ap = tokens[1]
            results = tokens[4]
            allele, peptide = ap.split("_")
            results = {x.strip().split(": ")[0]: x.strip().split(": ")[1]
                       for x in results.split("Matched")[1:]}
            tmpseries = pd.Series(results)
            tmpseries["allele"] = allele
            tmpseries["peptide"] = peptide
            tf_results = tf_results.append(tmpseries, ignore_index=True)
    return(tf_results)


def get_hla_seq(hla_df, allele):
    returndf = hla_df[[True if allele in x else False for x in hla_df.index]]
    # first seq used.
    returnseq = returndf.iloc[0]
    return returnseq


def check_dot(mhc_seq, pdb_seq):
    '''
        Check dot in MHC seq and PDB allele seq.
    '''
    for a, b in zip(mhc_seq, pdb_seq):
        yield True if (a == "." and b != ".") or (a != "." and b == ".") else False


def calculate_blosum(seq1, seq2, blosum, gap_s, gap_e, gap=True):
    for A, B in zip(list(seq1), list(seq2)):
        if A == "." and B == ".":
            continue
        diag = (A == '.') or (B == '.') or (A == " ") or (B == " ")
        yield (gap_e if gap else gap_s) if diag else blosum[A][B]
        gap = diag


def calculate_pam(pep1, pep2, pam):
    '''
        Peptide's BLOSUM score calculator
    '''
    for a, b in zip(list(pep1), list(pep2)):
        yield pam[a][b]


# %% filepath setting.
imgt_fn = "./reference/IMGT_structure_summary.csv"
hla_fn = "./reference/HLA_template.pickle"
pdb_fn = "./reference/PDB_template.pickle"
blosum_fn = "./reference/alignment_matrix/BLOSUM62"
pam_fn = "./reference/alignment_matrix/PAM30"
pep_fn = "./script_validation/peptide_list.txt"

# data calling
templateinfo = pd.read_csv(imgt_fn, index_col=0)
hla_df = pd.read_pickle(hla_fn, compression="gzip")
pdb_df = pd.read_pickle(pdb_fn, compression="gzip")
blosum = pd.read_csv(blosum_fn, header=0, skiprows=6,
                     index_col=0, sep="\\s+")
pam = pd.read_csv(pam_fn, header=0, skiprows=9,
                  index_col=0, sep="\\s+")

# set target allele, peptide.
tgt_alleles = ["A2402", 'A0201', "A0362", "A1101", "A7402", "A2301", "A3004", "A6835",
               'B5101', 'B1501', "B4402", "B5801", "B0762", "B1539", "B3530", "B8202",
               "C1802", "C0812", "C0303", "C0702"]  # 20 allele
with open(pep_fn) as f:
    tgt_peptides = [x.strip() for x in f.readlines()]
# peptide selection
tgt_peptides = sample(tgt_peptides, 100)
logger = make_logger("./script_validation", name="test_tf", consoleset=False)

# run class
tf_results = pd.DataFrame()
for allele in tqdm(tgt_alleles):
    for peptide in tgt_peptides:
        tf = find_template.TemplateFinder(allele, peptide, hla_fn, pdb_fn,
                                          imgt_fn, blosum_fn, pam_fn, logger)
        pdb, blosum_score, pam_score = tf()
        tmpseries = pd.Series([allele, peptide, pdb, blosum_score, pam_score],
                              index=["allele", "peptide", "pdb", "blosum_score", "pam_score"])
        tf_results = tf_results.append(tmpseries, ignore_index=True)

tf_results_bak = tf_results.copy()
# clear value.
del allele, peptide, tf, pdb, blosum_score, pam_score, tmpseries

# draw plot
fig, axes = plt.subplots(nrows=4, ncols=5, figsize=(16, 15))
allelelist = tf_results["allele"].unique()
for i in range(len(allelelist)):
    ncol = i % 5
    nrow = floor(i/5)
    allele = allelelist[i]
    tmp_results = tf_results[tf_results["allele"].isin([allele])]
    allele_2 = f"{allele[0]}*{allele[1:3]}:{allele[3:5]}"
    hla_seq = get_hla_seq(hla_df, allele_2)
    hla_seq2 = "".join(hla_seq.tolist()).replace(".", "")
    matched_id = [i for i, row in pdb_df.iterrows()
                  if sum(check_dot(hla_seq, row)) == 0]
    # dot matching template filtering.
    tmp_pdb = pdb_df.copy()
    tmp_pdb["dot_check"] = [sum(check_dot(hla_seq, row)) for i, row in tmp_pdb.iterrows()]
    tmp_pdb = tmp_pdb[tmp_pdb["dot_check"].isin([0])]
    matched_pdb = tmp_pdb.index.tolist()
    # sequence blosum check.
    tmp_info = templateinfo.loc[matched_pdb]
    # peptide length matcher get.
    tmp_info = tmp_info[tmp_info["p_len"].isin([9])]
    tmp_info["blosum_score"] = [sum(calculate_blosum(hla_seq2, x,
                                                     blosum, -11, 1))
                                for x in tmp_info["mhc_seq"]]
    bins = np.linspace(tmp_info["blosum_score"].min(),
                       tmp_info["blosum_score"].max(),
                       int(len(tmp_info)/5))
    axes[nrow][ncol].hist(tmp_info["blosum_score"], bins=bins, density=True,
             color="gray", alpha=0.5, label="PDB score")
    axes[nrow][ncol].hist(tmp_results["blosum_score"], bins=bins, density=True,
             color="blue", alpha=0.5, label="Matched PDB score")
    axes[nrow][ncol].legend()
    axes[nrow][ncol].set_title(allele)
plt.show()

def get_allele_pdb_list(tf_results):
    returndict = {}
    for allele in tf_results["allele"].unique():
        tmpdf = tf_results[tf_results['allele'].isin([allele])]
        pdblist = tmpdf["pdb"].unique().tolist()
        returndict[allele] = pdblist
    return returndict

def get_max_idx(tmpseries):
    maxval = tmpseries.max()
    idxlist = [x for x in tmpseries.index if tmpseries[x] == maxval]
    # join
    returnval = ",".join(idxlist)
    return returnval

def make_pdb_df(tf_results, templateinfo, pam):
    returndf = pd.DataFrame()
    returndict = {}
    allele_pdb_dict = get_allele_pdb_list(tf_results)
    for allele, templatelist in allele_pdb_dict.items():
        # allele results.
        tmpdf = tf_results[tf_results['allele'].isin([allele])]
        for pdb in templatelist:
            p = templateinfo["p_seq"][pdb]
            tmpdf[pdb] = [sum(calculate_pam(x, p, pam)) for x in tmpdf["peptide"]]
        # print(tmpdf.iloc[0])
        # print(templatelist)
        # print(tmpdf[templatelist])
        # d = tmpdf.iloc[0][templatelist]
        # print()
        # print(d)
        # print(type(d))
        # print()
        # print(d.idxmax())
        tmpdf["max_id"] = [get_max_idx(row[templatelist].astype('int64')) for i, row in tmpdf.iterrows()]
        tmpdf["matched_id_is_max_id"] = [True if row["pdb"] in row["max_id"] else False for i, row in tmpdf.iterrows()]
        returndict[allele] = tmpdf
        tmpseries = pd.Series([allele, len(templatelist), len(tmpdf), sum(tmpdf["matched_id_is_max_id"])],
                              index=["allele", "n_template", "n_peptide", "n_matched"])
        returndf = returndf.append(tmpseries, ignore_index=True)
    return returndf, returndict

tf_pep_summary, tf_dict = make_pdb_df(tf_results, templateinfo, pam)

# %% test code
tmp_info["mhc_seq"]["4u1h_A"]

allele = "B0789"
peptide = tgt_peptides[0]
tf = find_template.TemplateFinder(allele, peptide, hla_fn, pdb_fn,
                                  imgt_fn, blosum_fn, pam_fn, logger)
pdb, blosum, pam = tf()



allele_df = hla_df.loc[[True if x.endswith(":01:01") else False for x in hla_df.index]]
sample(allele_df.index.tolist(), 20)

tf_results["pam_score"].idxmax()



test_series = pd.Series([1,1,2,3,4,1,2,4,1],
                        index=["A","B", "C", "D", 'E', 'F', 'G', "H", "I"])
test_series.idxmin()











