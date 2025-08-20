#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 15:47:51 2021

@author: duaghk

Build template.

find_template return:
        returnlist = {
            "pdb": self.matched_pdb,
            "HLA_score": self.max_pdb_score,
            "HLA_score_ratio": round(self.max_pdb_score/self.best_blosum_score, 3),
            "Peptide_score": self.max_pam_score,
            "Peptide_score_ratio": round(self.max_pam_score/self.best_pam_score, 3)
            }

"""

# import library
import os
import subprocess as sp
from time import time
import psutil


# define class
class BuildTemplate:
    '''
        Build template.
        STRUMP-I need to modify template allele, peptide chain sequence.
        Using this module, best-matched pdb were changed via foldx buildModel.
    '''

    def __init__(self, pdb: str, indiv_fn, outdir: str, foldx: str, logger):
        self.pdb = pdb
        self.indiv_fn = indiv_fn
        self.outdir = outdir
        self.foldx = foldx
        self.logger = logger

    def run_buildmodel(self):
        '''
            Foldx BuildModel running.
        '''
        cmd = (
            f"{self.foldx} --command=BuildModel "
            f"--pdb-dir={self.outdir} --pdb={self.pdb} "
            f"--mutant-file={self.indiv_fn} "
            f"--output-dir={self.outdir}"
            )
        self.logger.info(f"Build Template command: {cmd}")
        # run process.
        timeout_check = 0
        with sp.Popen(cmd, shell=True, stdout=sp.PIPE, universal_newlines=True) as proc:
            try:
                outs = proc.communicate(timeout=3600)
            except sp.TimeoutExpired:
                proc.kill()
                outs = proc.communicate()
                timeout_check = 1
        self.logger.info("foldx BuildModel log:")
        self.logger.info("\n" + outs[0] + "\n")
        if timeout_check:
            p_proc = psutil.Process(os.getppid())
            p_proc.terminate()
            raise SystemExit("Process terminated due to Timeout")
        # moving template
        template_name = os.path.splitext(self.pdb)[0]
        buildmodel_out = os.path.join(self.outdir, "builded_template.pdb")
        os.rename(os.path.join(self.outdir, f"{template_name}_1.pdb"), buildmodel_out)
        self.logger.info(f"foldx BuildModel out: {buildmodel_out}")
        return buildmodel_out

    def __call__(self):
        self.logger.info("Template building start...")
        start = time()
        return_fn = self.run_buildmodel()
        end = time()
        self.logger.info(f"Template building completed. "
                         f"Time elapsed: {round((end-start)/60, 3)} min")

        return return_fn


if __name__ == "__main__":
    import argparse
    import logging
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb")
    parser.add_argument("--indiv_fn")
    parser.add_argument("--outdir")
    parser.add_argument("--foldx")
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
    logger = make_logger(os.path.dirname(args.outdir), name="test_build_energy")
    build_tempt = BuildTemplate(args.pdb, args.indiv_fn, args.outdir, args.foldx, logger)
    build_tempt()
