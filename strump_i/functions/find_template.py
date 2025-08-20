#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 11:16:08 2021

@author: duaghk

STRUMP-I template finding modules.

"""

import pandas as pd
from random import sample


class TemplateFinder:
    '''
        Most matched template finder.
        return dataframe
        peptide:template
    '''

    def __init__(self, allele, peptide,
                 mhc_seq_mat_fn, pdb_seq_mat_fn, templateinfo_fn,
                 blosum_fn, pam_fn, logger):
        # peptide list input.
        self.allele = f"{allele[0]}*{allele[1:3]}:{allele[3:5]}"
        self.mhc = None
        self.peptide = peptide
        self.mhc_seq_mat = pd.read_pickle(mhc_seq_mat_fn, compression="gzip")
        self.pdb_seq_mat = pd.read_pickle(pdb_seq_mat_fn, compression="gzip")
        self.templateinfo = pd.read_csv(templateinfo_fn, index_col=0)
        self.blosum = pd.read_csv(blosum_fn, header=0, skiprows=6,
                                  index_col=0, sep="\\s+")
        self.pam = pd.read_csv(pam_fn, header=0, skiprows=9,
                               index_col=0, sep="\\s+")
        self.gap_s = -11
        self.gap_e = -1
        self.logger = logger

    def find_mhc_seq(self):
        '''
            Find mhc seq in mhc_seq_mat
        '''
        self.mhc_seq_mat = self.mhc_seq_mat.loc[[x for x in self.mhc_seq_mat.index
                                                 if self.allele in x]]
        self.mhc = self.mhc_seq_mat.iloc[0]

    def check_dot(self, mhc_seq, pdb_seq):
        '''
            Check dot in MHC seq and PDB allele seq.
        '''
        for a, b in zip(mhc_seq, pdb_seq):
            yield True if (a == "." and b != ".") or (a != "." and b == ".") else False

    def calculate_pam(self, pep1, pep2):
        '''
            Peptide's BLOSUM score calculator
        '''
        for a, b in zip(list(pep1), list(pep2)):
            yield self.pam[a][b]

    def calculate_blosum(self, seq1, seq2, gap_s, gap_e, gap=True):
        for A, B in zip(list(seq1), list(seq2)):
            if A == "." and B == ".":
                continue
            diag = (A == '.') or (B == '.') or (A == " ") or (B == " ")
            yield (gap_e if gap else gap_s) if diag else self.blosum[A][B]
            gap = diag

    def filter_pdb(self, pdb_list):
        '''
            If seq position not matched, need to filter.
        '''
        # mhc seq filtering first.
        matched_id = [i for i, row in self.pdb_seq_mat.iterrows()
                      if sum(self.check_dot(self.mhc, row)) == 0]
        # filter check pep len equal.
        matched_id = [x for x in matched_id if x in self.templateinfo.index]
        matched_id = [x for x in matched_id if x not in pdb_list]
        self.templateinfo = self.templateinfo.loc[matched_id]
        self.templateinfo = self.templateinfo[self.templateinfo['p_len'].isin([len(self.peptide)])]
        # intersect templateinfo.index to pdb_seq_mat
        self.pdb_seq_mat = self.pdb_seq_mat[self.pdb_seq_mat.index.isin(self.templateinfo.index)]
        # template blosum score.
        self.pdb_score = {i: sum(self.calculate_blosum(self.mhc, row, self.gap_s, self.gap_e))
                          for i, row in self.pdb_seq_mat.iterrows()}
        # calculate all equal sequence's score.
        self.best_blosum_score = sum(self.calculate_blosum(self.mhc, self.mhc,
                                                           self.gap_s, self.gap_e))
        self.max_pdb_score = max(self.pdb_score.values())
        self.matched_pdb = [k for k, v in self.pdb_score.items()
                            if v == self.max_pdb_score]
        self.templateinfo = self.templateinfo[self.templateinfo.index.isin(self.matched_pdb)]

        # calculate peptide pam score
        self.templateinfo['peptide_pam'] = [sum(self.calculate_pam(self.peptide, x))
                                            for x in self.templateinfo['p_seq']]
        self.max_pam_score = self.templateinfo['peptide_pam'].max()
        self.best_pam_score = sum(self.calculate_pam(self.peptide, self.peptide))
        self.templateinfo = self.templateinfo[
            self.templateinfo['peptide_pam'].isin([self.max_pam_score])
            ]
        self.matched_pdb = self.templateinfo.index.tolist()
        if len(self.matched_pdb) > 1:
            self.matched_pdb = sample(self.matched_pdb, 1)
        self.matched_pdb = self.matched_pdb[0]

    def __call__(self, pdb_list: list):
        self.find_mhc_seq()
        self.filter_pdb(pdb_list)
        returndict = {
            "pdb": self.matched_pdb,
            "HLA_score": self.max_pdb_score,
            "HLA_score_ratio": round(self.max_pdb_score/self.best_blosum_score, 3),
            "Peptide_score": self.max_pam_score,
            "Peptide_score_ratio": round(self.max_pam_score/self.best_pam_score, 3)
            }
        tf_results_log = (
            f"Matched pdb id: {returndict['pdb']}, "
            f"Matched HLA BLOSUM score: {returndict['HLA_score']}, "
            f"Matched HLA BLOSUM ratio: {returndict['HLA_score_ratio']}, "
            f"Matched Peptide PAM score: {returndict['Peptide_score']}, "
            f"Matched Peptide PAM ratio: {returndict['Peptide_score_ratio']}"
            )
        self.logger.info(tf_results_log)
        return self.matched_pdb


if __name__ == "__main__":
    import os
    import logging
    import argparse
    parser = argparse.ArgumentParser()
    maindir = os.path.dirname(os.path.abspath(__file__))
    maindir = os.path.dirname(maindir)
    refdir = os.path.join(maindir, "reference")
    parser = argparse.ArgumentParser()
    parser.add_argument("--allele")
    parser.add_argument("--peptide")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mhc_seq_mat_fn", default=os.path.join(refdir,
                                                                 "HLA_template.pickle"))
    parser.add_argument("--pdb_seq_mat_fn", default=os.path.join(refdir,
                                                                 "PDB_template.pickle"))
    parser.add_argument("--templateinfo_fn", default=os.path.join(refdir,
                                                                  "IMGT_structure_summary.csv"))
    parser.add_argument("--blosum_fn", default=os.path.join(refdir, "alignment_matrix",
                                                            "BLOSUM62"))
    parser.add_argument("--pam_fn", default=os.path.join(refdir, "alignment_matrix",
                                                         "PAM30"))
    args = parser.parse_args()
    
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
    
    logger = make_logger(args.outdir, name="test_find_template")
    tf = TemplateFinder(args.allele, args.peptide, 
                        args.mhc_seq_mat_fn, args.pdb_seq_mat_fn, args.templateinfo_fn,
                        args.blosum_fn, args.pam_fn, logger)
    tf()



