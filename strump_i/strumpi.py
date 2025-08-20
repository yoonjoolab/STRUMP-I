#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 23:10:21 2021

@author: duaghk

run eachmodules.

"""

# import library
import os
import logging
import pandas as pd
from time import time
from pathlib import Path
from strump_i.functions import find_template, build_template, optimize_template, calculate_energy


class STRUMPI:
    def __init__(self, pdbdir, outdir, tinker_dir, foldx,
                 mhc_seq_mat_fn, pdb_seq_mat_fn, templateinfo_fn,
                 blosum_fn, pam_fn, keyfile_fn, rms):
        self.allele = ""
        self.peptide = ""
        self.pdbdir = pdbdir
        self.outdir = outdir
        self.pdbxyz = os.path.join(tinker_dir, "pdbxyz")
        self.optimizer = os.path.join(tinker_dir, "minimize")
        self.xyzpdb = os.path.join(tinker_dir, "xyzpdb")
        self.foldx = foldx
        self.mhc_seq_mat_fn = mhc_seq_mat_fn
        self.pdb_seq_mat_fn = pdb_seq_mat_fn
        self.templateinfo_fn = templateinfo_fn
        self.blosum_fn = blosum_fn
        self.pam_fn = pam_fn
        self.keyfile_fn = keyfile_fn
        self.rms = rms

    def make_logger(self, ap_dir, name=None, consoleset=True):
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
        loggerfile = os.path.join(ap_dir, name)
        file_handler = logging.FileHandler(filename=f"{loggerfile}.log")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        return logger

    def parse_matching_template(self, pdb):
        hla_df = pd.read_pickle(self.mhc_seq_mat_fn, compression="gzip")
        pdb_df = pd.read_pickle(self.pdb_seq_mat_fn, compression="gzip")
        allele_name = f"{self.allele[0]}*{self.allele[1:3]}:{self.allele[3:5]}"
        # find allele.
        hla_name = [x for x in hla_df.index if allele_name in x][0]
        pdb_name = [x for x in pdb_df.index if pdb in x][0]
        returndf = pd.DataFrame({"hla": hla_df.loc[hla_name],
                                 "pdb": pdb_df.loc[pdb_name]})
        return returndf

    def generate_template_mutate(self, template_df: pd.DataFrame):
        # need to remove . aa
        # make condition.
        dot_cond = (template_df['hla'] != ".") & (template_df['pdb'] != ".")
        template_df = template_df[dot_cond]
        for i in range(len(template_df)):
            if template_df['hla'][i] != template_df['pdb'][i]:
                yield f"{template_df['pdb'][i]}A{i+1}{template_df['hla'][i]}"

    def get_template_peptide(self, template_fn):
        # set aa_dict
        aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        # get template_peptide list.
        peplist = [x for x in open(template_fn).readlines() if x[21] == "P"]
        pepdict = {}
        for line in peplist:
            aa = aa_dict[line[17:20]]
            pos = line[22:26].strip()
            if pos not in pepdict.keys():
                pepdict[int(pos)] = aa
        pepseries = pd.Series(pepdict)
        template_pep = "".join(pepseries)
        return template_pep

    def generate_peptide_mutate(self, seq1: str, seq2: str):
        for i in range(len(seq1)):
            yield f"{seq1[i]}P{i+1}{seq2[i]}"

    def generate_indiv_fn(self, pdb):
        # make individual_list.
        template_df = self.parse_matching_template(pdb)
        template_mut = list(self.generate_template_mutate(template_df))

        # get pdb file name.
        pdbfile = [x for x in os.listdir(self.pdbdir) if pdb in x][0]
        template_pep = self.get_template_peptide(os.path.join(self.pdbdir, pdbfile))
        pep_mut = list(self.generate_peptide_mutate(template_pep, self.peptide))

        # make individual files.
        mut_all = template_mut + pep_mut
        return pdbfile, mut_all

    def __call__(self, allele, peptide):
        self.allele = allele
        self.peptide = peptide
        start = time()
        prefix = f"{self.allele}_{self.peptide}"
        # make allele dir
        allele_dir = os.path.join(self.outdir, self.allele)
        ap_dir = os.path.join(allele_dir, prefix)
        Path(ap_dir).mkdir(parents=True, exist_ok=True)
        # self.outdir = ap_dir
        self.logger = self.make_logger(ap_dir, name=prefix)
        self.logger.info("STRUMP-I start. "
                         f"target allele: {self.allele} "
                         f"target peptide: {self.peptide}"
                         )

        # code start.
        # template finding.
        primary_e = 999
        count = 0
        pdb_list = []
        while primary_e > 0:
            # file remove.
            file_list = list(Path(ap_dir).iterdir())
            remove_file_list = [x for x in file_list if not str(x).endswith(".log")]
            [x.unlink() for x in remove_file_list]
            count += 1
            if count > 5:
                break
            TF = find_template.TemplateFinder(
                self.allele, self.peptide,
                self.mhc_seq_mat_fn, self.pdb_seq_mat_fn, self.templateinfo_fn,
                self.blosum_fn, self.pam_fn, self.logger)
            pdb = TF(pdb_list)
            pdb_list.append(pdb)
            pdb = pdb.split("_")[0]

            # make individual files and copy pdb to outdir
            pdb_fn, mut_all = self.generate_indiv_fn(pdb)
            # pdb file copy
            pdb_cp_cmd = f"cp {os.path.join(self.pdbdir, pdb_fn)} {ap_dir}"
            os.system(pdb_cp_cmd)
            indiv_fn = os.path.join(ap_dir, "individual_list.txt")
            with open(indiv_fn, 'w') as f:
                f.write(",".join(mut_all) + ";")

            # build template
            BT = build_template.BuildTemplate(pdb_fn, indiv_fn, ap_dir,
                                            self.foldx, self.logger)
            builded = BT()

            # optimize template
            OT = optimize_template.OptimizeTemplate(
                builded, self.pdbxyz, self.optimizer, self.xyzpdb, self.foldx,
                self.keyfile_fn, self.rms, self.logger)
            optimized, primary_e = OT()


        # calculate energy
        CE = calculate_energy.CalculateEnergy(
            optimized, self.foldx, self.logger)
        CE()
        end = time()
        self.logger.info("STRUMP-I completed. "
                         f"Time elapsed: {round((end-start)/60, 3)} min")
        self.logger.handlers.clear()


if __name__ == "__main__":
    import argparse
    maindir = os.path.dirname(os.path.abspath(__file__))
    refdir = os.path.join(maindir, "reference")
    parser = argparse.ArgumentParser()
    parser.add_argument("--allele")
    parser.add_argument("--peptide")
    parser.add_argument("--pdbdir", default=os.path.join(refdir, "pdb"))
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--tinker_dir", default=os.path.dirname("/tinker"))
    parser.add_argument("--foldx", default=os.path.dirname("/foldx"))
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
    parser.add_argument("--keyfile_fn", default=os.path.join(refdir, "force.key"))
    parser.add_argument("--rms", default=1)
    args = parser.parse_args()
    strumpi = STRUMPI(
        args.pdbdir, args.outdir, args.tinker_dir, args.foldx,
        args.mhc_seq_mat_fn, args.pdb_seq_mat_fn, args.templateinfo_fn,
        args.blosum_fn, args.pam_fn, args.keyfile_fn, args.rms)
    strumpi(args.allele, args.peptide)
