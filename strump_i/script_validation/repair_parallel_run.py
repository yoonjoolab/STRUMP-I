#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 08:49:52 2021

@author: duaghk

Repair compare

run_repair : previous
run_repair2 : new

"""

import os
import sys
import subprocess as sp
import multiprocessing as mp
from itertools import repeat

def run_repair(foldx, atom_changed_fn):
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
        print(f"Repairing command: {cmd}")
        with sp.Popen(cmd, shell=True) as proc:
            proc.wait()
        # proc = sp.Popen(cmd, shell=True)
        # proc.wait()
        with open(repaired_fxout) as f_read:
            line = f_read.readlines()[-1]
        cur_e = float(line.split()[1])
        cur_e = round(cur_e, 4)
        if count == 1:
            start_e = cur_e
        print(
            f"Previous energy value: ({prev_e}) "
            f"and the current energy: ({cur_e})"
            )
        if cur_e > prev_e - 0.001:
            print(
                f"{prev_e} -> {cur_e} (Improvement {round(prev_e-cur_e, 4)}). "
                "Repairing is done."
                )
            break
        count += 1
        print(
            f"{prev_e} -> {cur_e} (Improvement {round(prev_e-cur_e, 4)}). "
            f"Foldx Repair is performed again (Attempt {count})"
            )
        prev_e = cur_e
        mvcmd = (
            f"mv {os.path.join(outdir, repaired_fn)} "
            f"{os.path.join(outdir, optimized_fn)}"
            )
        print(f"Repaired file moving command: {mvcmd}")
        os.system(mvcmd)
    # save time.
    maindir = os.path.dirname(outdir)
    with open(os.path.join(maindir, "repair_results.csv"), 'a') as f_write:
        f_write.write(f"{os.path.join(outdir, repaired_fn)},{start_e},{cur_e},{count}")


def run_repair2(foldx, atom_changed_fn):
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
        cmd = (
            f"{foldx} --command=RepairPDB "
            f"--pdb-dir={outdir} --pdb={optimized_fn} "
            f"--screen 0 --output-dir={outdir} "
            )
        with sp.Popen(cmd, shell=True) as proc:
            proc.wait()
        with open(repaired_fxout) as f_read:
            line = f_read.readlines()[-1]
        cur_e = float(line.split()[1])
        cur_e = round(cur_e, 3)
        if count == 1:
            start_e = cur_e
        print(
            f"Previous energy value: ({best_e}) / "
            f"and the current energy: ({cur_e})"
            )
        print(
            f"{best_e} -> {cur_e} (Improvement {round(best_e-cur_e, 3)}). "
            )

        if cur_e < best_e:
            best_e = cur_e
            cpcmd = f"cp {os.path.join(outdir, repaired_fn)} {os.path.join(outdir, best_fn)}"
            os.system(cpcmd)
            stopping_count = 0
        else:
            stopping_count += 1
            if stopping_count >= 5:
                print("Repairing is done.")
                break
        mvcmd = (
            f"mv {os.path.join(outdir, repaired_fn)} "
            f"{os.path.join(outdir, optimized_fn)}"
            )
        print(f"Repaired file moving command: {mvcmd}")
        os.system(mvcmd)
        count += 1
        print(
            f"Foldx Repair is performed again (Attempt {count}) /"
            f" (Stopping count {stopping_count})"
            )
    # save time.
    maindir = os.path.dirname(outdir)
    with open(os.path.join(maindir, "repair_results.csv"), 'a') as f_write:
        f_write.write(f"{os.path.join(outdir, best_fn)},{start_e},{best_e},{count}")
    # return os.path.join(outdir, best_fn), start_e, best_e, count


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtype", choices=["previous", "new"])
    parser.add_argument("--inputdir")  # allele_dir
    parser.add_argument("--threads", default=10, type=int)

    args = parser.parse_args()
    foldx = "~/bin/foldx"
    ap_list = os.listdir(args.inputdir)
    tgt_list = [os.path.join(args.inputdir, x, "builded_template_minim.pdb") for x in ap_list]
    if args.runtype == "previous":
        with mp.Pool(args.threads) as pool:
            pool.starmap(run_repair, zip(repeat(foldx), tgt_list))
    else:
        with mp.Pool(args.threads) as pool:
            pool.starmap(run_repair2, zip(repeat(foldx), tgt_list))
