#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 08:45:19 2021

@author: duaghk

STRUMP-I Quality data parser.

"""

# import library
import os
import numpy as np
import pandas as pd


# define functions.

class ParseQuality:
    '''
        STRUMP-I Quality data parser
    '''

    def __init__(self, input_dir, output_dir):
        '''
            Initiation.
            need value : inputdir, outputdir
        '''
        self.indir = input_dir
        self.outdir = output_dir

    def parse_final_rms(self, lines):
        '''
            Line parser for optimized RMS value.
        '''
        returndict = {}
        tgt_line = [x for x in lines if "Final RMS Gradient" in x]
        if len(tgt_line):
            final_rms = tgt_line[0].split()[-1]
            returndict["minimized_RMS"] = final_rms
        else:
            returndict["minimized_RMS"] = ""
        return returndict

    def parse_repair_energy(self, lines):
        '''
            Line parser for repaired energy value and iteration times.
        '''
        returndict = {}
        # get first energy
        tgt_lines = [x for x in lines if "Improvement" in x]
        iter_lines = [x for x in lines if "Attempt" in x]
        if len(tgt_lines):
            tokens = [x.split("> ")[1].split(" (")[0] for x in tgt_lines]
            primary_e = tokens[0]
            final_e = min(map(float, tokens))
            iter_ns = [int(x.split("Attempt ")[1].split(")")[0]) for x in iter_lines]
            iter_n = max(iter_ns)
            returndict["primary_e"] = primary_e
            returndict["final_e"] = final_e
            returndict["iter_n"] = iter_n
        else:
            returndict["primary_e"] = ""
            returndict["final_e"] = ""
            returndict["iter_n"] = ""
        return returndict

    def parse_times(self, lines):
        '''
            Line parser for each step's elapsed times
        '''
        returndict = {}
        tgt_modules = ["build_template", "optimize_template", "strumpi"]
        tgt_lines = [x for x in lines if "Time elapsed:" in x]
        for module in tgt_modules:
            tgt_line = [x for x in tgt_lines if module in x]
            module_name = f"{module}_time"
            if len(tgt_line):
                times = tgt_line[0].split("elapsed: ")[-1]
                times = times.split()[0]
                returndict[module_name] = times
            else:
                returndict[module_name] = ""
        return returndict

    def parse_log(self, ap_dir):
        '''
            one allele-peptide pair log parser
        '''
        prefix = os.path.basename(ap_dir)
        tmpdict = {"Pdb": prefix}
        tgt_fn = [x for x in os.listdir(ap_dir) if x.endswith(".log")][0]
        with open(os.path.join(ap_dir, tgt_fn)) as f_read:
            lines = f_read.readlines()
        optimize_dict = self.parse_final_rms(lines)
        repair_dict = self.parse_repair_energy(lines)
        time_dict = self.parse_times(lines)
        tmpdict.update(optimize_dict)
        tmpdict.update(repair_dict)
        tmpdict.update(time_dict)
        returnseries = pd.Series(tmpdict, name=prefix)
        return returnseries

    def parse_all_log(self):
        '''
            Parse all allele-peptide pairs.
        '''

        qual_df = {os.path.basename(root): self.parse_log(root) for root, dir, files in os.walk(self.indir) if len([x for x in files if x.endswith(".log")]) > 0}
        qual_df = pd.DataFrame.from_dict(qual_df, orient="index")
        qual_df["Pdb"] = qual_df.index
        # remained parsing.
        qual_df = qual_df.replace("", np.NaN)
        qual_df["allele"] = qual_df["Pdb"].apply(lambda x: x.split("_")[0])
        qual_df["peptide"] = qual_df["Pdb"].apply(lambda x: x.split("_")[1])
        qual_df["primary_e"] = qual_df["primary_e"].astype(float)
        qual_df["final_e"] = qual_df["final_e"].astype(float)
        qual_df["repair_improved"] = qual_df["primary_e"] - qual_df["final_e"]
        # column reordering
        column_order = [
            'Pdb', 'allele', 'peptide', 'minimized_RMS',
            'primary_e', 'final_e', 'repair_improved', 'iter_n',
            'build_template_time', 'optimize_template_time', 'strumpi_time'
            ]
        qual_df = qual_df[column_order]
        return qual_df

    def __call__(self):
        quality_df = self.parse_all_log()
        save_fn = os.path.join(self.outdir, f"quality_data.csv")
        quality_df.to_csv(save_fn, header=True, index=False)
        pass

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", "-I", required=True)
    parser.add_argument("--output_dir", "-O", required=True)
    args = parser.parse_args()
    q_parser = ParseQuality(args.input_dir, args.output_dir)
    q_parser()