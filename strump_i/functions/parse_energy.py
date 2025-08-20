#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 08:33:58 2021

@author: duaghk

STRUMP-I Energy data parser

"""

# import library
import os
import pandas as pd
from string import ascii_uppercase as abt


# define functions.

class ParseEnergy:
    '''
        STRUMP-I Energy output parser.
    '''

    def __init__(self, input_dir, output_dir):
        self.indir = input_dir
        self.outdir = output_dir

    def check_skiplines(self, fn):
        count = 0
        with open(fn) as f:
            while True:
                line = f.readline()
                if line.startswith("Pdb"):
                    break
                count += 1
        return count

    def parse_data(self, fn):
        skipnum = self.check_skiplines(fn)
        df = pd.read_csv(fn, sep="\t", skiprows=skipnum)
        df["Pdb"] = [x.split("/")[-2] for x in df["Pdb"]]
        df = df[df["Group1"].isin(["A"])]
        return df

    def chain_pos_dict_gen(self, tgt_n: int):
        chain = abt[15:15+tgt_n]
        pos = list([x+1 for x in range(tgt_n)])
        returndict = {chain[i]: pos[i] for i in range(len(pos))}
        return returndict

    def parse_all(self, start_n, chain_n):
        returndf = [self.parse_data(os.path.join(root, [x for x in files if x.startswith(start_n) and chain_n in x][0]))
                    for root, dir, files in os.walk(self.indir) if len([x for x in files if x.startswith(start_n) and chain_n in x]) > 0]
        returndf = pd.concat(returndf)
        returndf = returndf.reset_index(drop=True)
        return returndf

    def merge_and_reordering(self, e_df, sum_e_df):
        sum_cols = ["Pdb", "Group1", "Group2", 'StabilityGroup1', 'StabilityGroup2']
        return_df = pd.merge(e_df, sum_e_df[sum_cols], on=["Pdb", "Group1", "Group2"])
        return_df["allele"] = [x.split("_")[0] for x in return_df["Pdb"]]
        return_df["peptide"] = [x.split("_")[1] for x in return_df["Pdb"]]
        infocol = ["Pdb", 'allele', 'peptide', "Group1", "Group2"]
        elsecol = return_df.columns.drop(infocol).tolist()
        return_df = return_df[infocol + elsecol]
        return return_df

    def __call__(self):
        pep_e = self.parse_all("Interaction", "Repair")
        pep_sum_e = self.parse_all("Summary", "Repair")
        pos_e = self.parse_all("Interaction", "modified")
        pos_sum_e = self.parse_all("Summary", "modified")

        pep_e = self.merge_and_reordering(pep_e, pep_sum_e)
        pos_e = self.merge_and_reordering(pos_e, pos_sum_e)

        pep_e_fn = os.path.join(self.outdir, f"peptide_energy_data.csv")
        pos_e_fn = os.path.join(self.outdir, f"chain_wise_energy_data.csv")
        pep_e.to_csv(pep_e_fn, header=True, index=False)
        pos_e.to_csv(pos_e_fn, header=True, index=False)
        return pep_e_fn


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", "-I", required=True)
    parser.add_argument("--output_dir", "-O", required=True)
    args = parser.parse_args()
    e_parser = ParseEnergy(args.input_dir, args.output_dir)
    e_parser()
    pass
