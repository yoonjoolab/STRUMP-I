#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 15:47:51 2021

@author: duaghk

energy calculation

"""

# import library
import os
import subprocess as sp
from time import time
from string import ascii_uppercase as abt


# define functions
class CalculateEnergy:
    def __init__(self, input_fn: str, foldx: str, logger):
        '''

        Args:
            input_fn (str): repair completed pdb file path. need to full path.
            foldx (str): foldx tool path.
            logger (logger): logging inherited from main script.

        Returns:
            None.

        '''
        self.outdir = os.path.dirname(input_fn)
        self.input_fn = os.path.basename(input_fn)
        self.foldx = foldx
        self.logger = logger

    def chain_pos_dict_gen(self, tgt_n: int):
        '''

        Args:
            tgt_n (int): Peptide length setting. target peptide length: 8~11.
                        so, do not over length than 12.

        Returns:
            returndict (dict): chain alphabet: integer.
                        peptide position to chain alphabet.

        '''
        chain = abt[15:15+tgt_n]
        pos = [x+1 for x in range(tgt_n)]
        returndict = {str(pos[i]): chain[i] for i in range(len(pos))}
        return returndict

    def change_chain(self, tgt_lines: list):
        '''

        Args:
            tgt_lines (list): target lines need to change chain ID.
                            chain ID == P (position 21) are in.

        Returns:
            returnline (list): chain ID changed by position number.
                            e.g. pos 2: P -> Q
                                 pos 3: P -> R
                                 etc.

        functions change peptide chain ID.

        '''
        returnline = []
        chain_num_dict = self.chain_pos_dict_gen(11)
        for line in tgt_lines:
            pos_num = line[22:26].strip()
            # change chain.
            line = list(line)
            line[21] = chain_num_dict[pos_num]
            line = "".join(line)
            returnline.append(line)
        return returnline

    def modify_pdblines(self, pdb_lines: list):
        '''

        Args:
            pdb_lines (list): all pdb lines.

        Returns:
            returnlines (list): return peptide chain ID changed all pdb lines.

        '''
        orig_lines = []
        tgt_lines = []
        for line in pdb_lines:
            if len(line) < 50:
                orig_lines.append(line)
            else:
                if line[21] == "P":
                    tgt_lines.append(line)
                else:
                    orig_lines.append(line)
        tgt_lines = self.change_chain(tgt_lines)
        # aggregate
        returnlines = orig_lines + tgt_lines
        return returnlines

    def rewrite_pdb(self):
        prefix = self.input_fn.split("_best")[0]
        pdb_fn = os.path.join(self.outdir, self.input_fn)
        # parse lines
        with open(pdb_fn) as f_read:
            orig_lines = f_read.readlines()
        modified_lines = self.modify_pdblines(orig_lines)
        modified_fn = f"{prefix}_peptide_chain_modified.pdb"
        with open(os.path.join(self.outdir, modified_fn), 'w') as f_write:
            f_write.write("".join(modified_lines))
        return modified_fn

    def calculate_energy(self, pdb):
        cmd = (
            f"{self.foldx} --command=AnalyseComplex "
            f"--pdb-dir={self.outdir} --pdb={pdb} "
            f"--output-dir={self.outdir} "
            f"--screen 0"
            )
        self.logger.info(f"Calculate energy command: {cmd}")
        with sp.Popen(cmd, shell=True) as proc:
            proc.wait()

    def __call__(self):
        self.logger.info("Calculate energy start...")
        start = time()
        # chain not modified energy calculation.
        self.logger.info("Repaired pdb energy calculation.")
        self.calculate_energy(self.input_fn)
        # chain modifying
        self.logger.info("Peptide chain modification.")
        chain_modified = self.rewrite_pdb()
        # chain modified energy calculation.
        self.logger.info("Chain-modified pdb energy calculation.")
        self.calculate_energy(chain_modified)
        end = time()
        self.logger.info("Calculate energy completed. "
                         f"Time elapsed: {round((end-start)/60, 3)} min")


if __name__ == "__main__":
    import argparse
    import logging
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fn")
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
    logger = make_logger(os.path.dirname(args.input_fn), name="calculate_energy_test")
    cal_e = CalculateEnergy(args.input_fn, args.foldx, logger)
    cal_e()
