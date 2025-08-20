#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 17:32:03 2021

@author: duaghk

test repairing script.

"""

import os
import sys
import logging
import subprocess as sp
import multiprocessing as mp
from itertools import repeat


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


def run_repair(foldx, atom_changed_fn, logger):
    '''
        Repairing process.
        iterated functions
    '''
    prev_e = sys.maxsize
    count = 1
    outdir = os.path.dirname(atom_changed_fn)
    optimized_fn = os.path.basename(atom_changed_fn)
    target = os.path.splitext(optimized_fn)[0]
    repaired_fn = f"{target}_Repair.pdb"
    repaired_fxout = os.path.join(outdir, f"{target}_Repair.fxout")
    while True:
        cmd = (
            f"{foldx} --command=RepairPDB "
            f"--pdb-dir={outdir} --pdb={optimized_fn} "
            f"--screen 0 --output-dir={outdir} "
            )
        logger.info(f"Repairing command: {cmd}")
        with sp.Popen(cmd, shell=True) as proc:
            proc.wait()
        # proc = sp.Popen(cmd, shell=True)
        # proc.wait()
        with open(repaired_fxout) as f_read:
            line = f_read.readlines()[-1]
        cur_e = float(line.split()[1])
        cur_e = round(cur_e, 4)
        logger.info(
            f"Previous energy value: ({prev_e}) "
            f"and the current energy: ({cur_e})"
            )
        if cur_e > prev_e - 0.001:
            logger.info(
                f"{prev_e} < {cur_e} (Improvement {round(prev_e-cur_e, 4)}). "
                "Repairing is done."
                )
            break
        count += 1
        logger.info(
            f"{prev_e} > {cur_e} (Improvement {round(prev_e-cur_e, 4)})."
            f"Foldx Repair is performed again (Attempt {count})"
            )
        prev_e = cur_e
        mvcmd = (
            f"mv {os.path.join(outdir, repaired_fn)} "
            f"{os.path.join(outdir, optimized_fn)}"
            )
        logger.info(f"Repaired file moving command: {mvcmd}")
        os.system(mvcmd)
    return os.path.join(outdir, repaired_fn)


def run(foldx, allele_dir, ap):
    allele, peptide = ap.split("_")
    ap_dir = os.path.join(allele_dir, ap)
    logger = make_logger(ap_dir, name=ap)
    # find input fn
    tgtlist = os.listdir(ap_dir)
    # get minim.pdb
    tgt_fn = [x for x in tgtlist if x.endswith("minim.pdb")][0]
    tgt_fn = os.path.join(ap_dir, tgt_fn)
    run_repair(foldx, tgt_fn, logger)
    logger.handlers.clear()
    pass


if __name__ == "__main__":
    maindir = "/data1/duaghk/strumpi/minimized/cut_0_01"
    foldx = "/home/duaghk/bin/foldx"
    tgtalleles = ['A2402', "B5101"]
    for allele in tgtalleles:
        allele_dir = os.path.join(maindir, allele)
        ap_list = os.listdir(allele_dir)
        with mp.Pool(50) as pool:
            pool.starmap(run, zip(repeat(foldx), repeat(allele_dir), ap_list))
