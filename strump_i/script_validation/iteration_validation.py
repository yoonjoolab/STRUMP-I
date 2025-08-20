#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 13:57:45 2021

@author: duaghk
"""

import os
import sys
import subprocess as sp
import multiprocessing as mp
from itertools import repeat

def run_repair(foldx, outdir):
    '''
        Repairing process.
        iterated functions
    '''
    atom_changed_fn = os.path.join(outdir, "builded_template_minim.pdb")
    # outdir = os.path.dirname(atom_changed_fn)
    best_e = sys.maxsize
    count = 1
    stopping_count = 0
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
    with open(os.path.join(outdir, "iter_change_results.txt"), 'w') as f_write:
        f_write.write(f"{count},{best_e}")
    # return os.path.join(outdir, best_fn), count, best_e


if __name__ == "__main__":
    tgt_dir = "/Users/duaghk/data/strumpi/script_test/"
    foldx = "/Users/duaghk/program/foldx"
    allele_list = ["A2402", "B5101"]
    for allele in allele_list:
        allele_dir = os.path.join(tgt_dir, allele)
        ap_list = os.listdir(allele_dir)
        ap_dir = [os.path.join(allele_dir, x) for x in ap_list]
        with mp.Pool(5) as pool:
            pool.starmap(run_repair, zip(repeat(foldx), ap_dir))    