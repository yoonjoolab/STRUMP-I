#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 15:47:51 2021

@author: duaghk

template optimization.

"""

# import library
import os
import sys
from time import time
import subprocess as sp
import psutil


class OptimizeTemplate:
    '''
        Run Tinker Optimize functions.
    '''

    def __init__(self, input_fn: str, pdbxyz: str, optimizer: str,
                 xyzpdb: str, foldx: str, keyfile: str, rms: float, logger):
        '''

        Args:
            inputfn (str): Template build completed file path.
            pdbxyz (str): tinker pdbxyz tool path.
            minimize (str): tinker minimize tool path.
            xyzpdb (str): tinker xyzpdb tool path
            keyfile (str): force.key file path. (tinker config file)
            logger : logger class for logging. logger inherited from main class.

        Returns:
            None.

        '''
        self.input_fn = input_fn
        self.pdbxyz = pdbxyz
        self.optimizer = optimizer
        self.xyzpdb = xyzpdb
        self.foldx = foldx
        self.keyfile = keyfile
        self.rms = rms
        self.logger = logger

    def run_pdbxyz(self):
        '''
            Run pdbxyz using class input.
        '''
        cmd = f"{self.pdbxyz} {self.input_fn} ALL ALL -k {self.keyfile}"
        self.logger.info(f"pdbxyz command: {cmd}")
        # process running.
        std_output = sp.check_output(cmd, shell=True, universal_newlines=True)
        self.logger.info("pdbxyz log:")
        self.logger.info(std_output + "\n")
        output_fn = os.path.splitext(self.input_fn)[0] + ".xyz"
        return output_fn

    def run_optimizer(self, pdbxyz_out: str):
        '''

        Args:
            pdbxyz_out (str): pdbxyz out file path. extension: xyz

        Raises:
            SystemExit: If optimizer run not well, e.g. BadIntpln error,
            program raise error and terminated.

        Returns:
            optimized (TYPE): optimized completed file path

        '''
        if self.optimizer.endswith("minimize"):
            cmd = f"{self.optimizer} {pdbxyz_out} -k {self.keyfile} {self.rms}"
        elif self.optimizer.endswith("newton"):
            cmd = f"{self.optimizer} {pdbxyz_out} -k {self.keyfile} A A {self.rms}"
        self.logger.info(f"optimize command: {cmd}")
        # process running.
        timeout_check = 0
        with sp.Popen(cmd, shell=True, stdout=sp.PIPE, universal_newlines=True) as proc:
            try:
                outs = proc.communicate(timeout=7200)
            except sp.TimeoutExpired:
                proc.kill()
                timeout_check = 1
                outs = proc.communicate()
        # output logging.
        self.logger.info("optimizer log:")
        self.logger.info(outs[0])
        # output processing.
        stdout = outs[0].split("\n")
        final_rms = stdout[-3].split(" ")[-1]
        final_rms = float(final_rms)
        if timeout_check:
            p_proc = psutil.Process(os.getppid())
            p_proc.terminate()
            raise SystemExit("Process termiate due to minimize timeout error")

        if final_rms > float(self.rms):
            error_reason = stdout[-6].split("--  ")[-1]
            self.logger.info(f"Optimization terminated due to: {error_reason}")
            raise SystemExit()
        optimized = f"{pdbxyz_out}_2"
        return optimized

    def run_xyzpdb(self, optimized_fn: str):
        '''

        Args:
            optimized_fn (str): optimized completed file path.

        Returns:
            return_fn (str): xyzpdb completed file path. extension: pdb

        '''
        cmd = f"{self.xyzpdb} {optimized_fn} -k {self.keyfile} ALL ALL"
        self.logger.info(f"xyzpdb command: {cmd}")
        # process running.
        std_output = sp.check_output(cmd, shell=True, universal_newlines=True)
        self.logger.info("pdbxyz log:")
        self.logger.info(std_output + "\n")
        output_fn = os.path.splitext(optimized_fn)[0] + ".pdb_2"
        return_fn = os.path.splitext(optimized_fn)[0] + "_minim.pdb"
        mv_cmd = f"mv {output_fn} {return_fn}"
        self.logger.info(f"xyzpdb out moving. command: {mv_cmd}")
        os.system(mv_cmd)
        return return_fn

    def change_atom(self, xyzpdb_fn):
        '''
            Atom name change function.
            Tinker tools use cystein as CYX instead CYS.
            Foldx not used CYX, so need to change CYX to CYS.
        '''
        with open(xyzpdb_fn) as f_read:
            lines = f_read.read()
        lines = lines.replace("CYX", "CYS")
        with open(xyzpdb_fn, 'w') as f_write:
            f_write.write(lines)
        return xyzpdb_fn

    def run_repair(self, atom_changed_fn):
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
                f"{self.foldx} --command=RepairPDB "
                f"--pdb-dir={outdir} --pdb={optimized_fn} "
                f"--screen 0 --output-dir={outdir} "
                )
            self.logger.info(f"Repairing command: {cmd}")
            with sp.Popen(cmd, shell=True) as proc:
                proc.wait()
            with open(repaired_fxout) as f_read:
                line = f_read.readlines()[-1]
            cur_e = float(line.split()[1])
            cur_e = round(cur_e, 3)
            self.logger.info(
                f"Previous energy value: ({best_e}) / "
                f"and the current energy: ({cur_e})"
                )
            self.logger.info(
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
            self.logger.info(f"Repaired file moving command: {mvcmd}")
            os.system(mvcmd)
            if count == 1:
                primary_e = cur_e
            if primary_e >= 0:
                break
            count += 1
            self.logger.info(
                f"Foldx Repair is performed again (Attempt {count}) /"
                f" (Stopping count {stopping_count})"
                )
        return os.path.join(outdir, best_fn), primary_e

    def __call__(self):
        self.logger.info("Template optimization start...")
        start = time()
        pdbxyz_out = self.run_pdbxyz()
        optimized_out = self.run_optimizer(pdbxyz_out)
        xyzpdb_out = self.run_xyzpdb(optimized_out)
        atom_changed_out = self.change_atom(xyzpdb_out)
        repaired_out, primary_e = self.run_repair(atom_changed_out)
        end = time()
        self.logger.info("Template optimization completed. "
                         f"Time elapsed: {round((end-start)/60, 3)} min")
        return repaired_out, primary_e


if __name__ == "__main__":
    import argparse
    import logging
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fn")
    parser.add_argument("--pdbxyz")
    parser.add_argument("--optimizer")
    parser.add_argument("--xyzpdb")
    parser.add_argument("--foldx")
    parser.add_argument("--keyfile")
    parser.add_argument("--rms", default=1)
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
    logger = make_logger(os.path.dirname(args.input_fn), name="optimizer_test")
    opti_tempt = OptimizeTemplate(args.input_fn, args.pdbxyz, args.optimizer,
                                  args.xyzpdb, args.foldx, args.keyfile, args.rms, logger)
    opti_tempt()
