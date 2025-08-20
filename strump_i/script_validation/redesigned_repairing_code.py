#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 10:30:46 2021

@author: duaghk
"""

import os
import sys
import pandas as pd
import subprocess as sp
from random import sample


def run_repair(foldx, atom_changed_fn):
    '''
        Repairing process.
        iterated functions
    '''
    best_e = sys.maxsize
    count = 1
    stopping_count = 0
    outdir = os.path.dirname(atom_changed_fn)
    optimized_fn = os.path.basename(atom_changed_fn)  # need to backup
    target = os.path.splitext(optimized_fn)[0]
    repaired_fn = f"{target}_Repair.pdb"
    repaired_fxout = os.path.join(outdir, f"{target}_Repair.fxout")
    best_fn = f"{target}_best_Repair.pdb"
    while True:
        print()
        print('**** step 1 *****')
        cmd = (
            f"{foldx} --command=RepairPDB "
            f"--pdb-dir={outdir} --pdb={optimized_fn} "
            f"--screen 0 --output-dir={outdir} "
            )
        # logger.info(f"Repairing command: {cmd}")
        with sp.Popen(cmd, shell=True) as proc:
            proc.wait()
        print('**** step 2 *****')
        with open(repaired_fxout) as f_read:
            line = f_read.readlines()[-1]
        cur_e = float(line.split()[1])
        cur_e = round(cur_e, 3)
        # logger.info(
        #     f"Previous energy value: ({prev_e}) "
        #     f"and the current energy: ({cur_e})"
        #     )
        print()
        print(f"{best_e} -> {cur_e} (Improvement {round(best_e-cur_e, 3)}). ")
        print()
        print('**** step 3 *****')
        if cur_e < best_e:
            print('**** step 4 *****')
            best_e = cur_e
            cpcmd = f"cp {os.path.join(outdir, repaired_fn)} {os.path.join(outdir, best_fn)}"
            os.system(cpcmd)
            stopping_count = 0
        else:
            print('**** step 5 *****')
            stopping_count += 1
            if stopping_count >= 5:
                print('**** step 6 *****')
                print("Repairing is done.")
                break
        print('**** step 7 *****')
        mvcmd = (
            f"mv {os.path.join(outdir, repaired_fn)} "
            f"{os.path.join(outdir, optimized_fn)}"
            )
        os.system(mvcmd)
        count += 1
        print()
        print(
            f"Foldx Repair is performed again (Attempt {count}) /"
            f" (Stopping count {stopping_count})"
            )        
        print()

    return os.path.join(outdir, best_fn), count, best_e

foldx = "/Users/duaghk/program/foldx"
atom_changed_fn = "/Users/duaghk/data/strumpi/script_test/B5101_RIRSERPAF/builded_template_minim.pdb"

run_repair(foldx, atom_changed_fn)

# random selected 10 file copy testdir.
src_dir = "/Users/duaghk/data/strumpi/new_strumpi_output/output/"
tgt_dir = "/Users/duaghk/data/strumpi/script_test/"
allele_list = ["A2402", "B5101"]
for allele in allele_list:
    allele_dir = os.path.join(src_dir, allele)
    ap_list = os.listdir(allele_dir)
    sampled_list = sample(ap_list, 10)
    for samples in sampled_list:
        cpcmd = f"cp -r {os.path.join(allele_dir, samples)} {os.path.join(tgt_dir, allele)}"
        os.system(cpcmd)

# run parallel

# itered data calling.
tgt_dir = "/Users/duaghk/data/strumpi/script_test/"
allele_list = ["A2402", "B5101"]
results_df = pd.DataFrame()
for allele in allele_list:
    allele_dir = os.path.join(src_dir, allele)
    ap_list = os.listdir(allele_dir)
    for ap in ap_list:
        




