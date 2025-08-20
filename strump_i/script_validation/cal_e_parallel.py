#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 18:45:20 2021

@author: duaghk
"""

# parallel processing.

import os
import logging
import multiprocessing as mp
from . import calculate_energy3
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

if __name__ == "__main__":
    maindir = "/data1/duaghk/test_new_strumpi/output/new/A0201"
    foldx = "~/bin/foldx"
    logger = make_logger(maindir, name="calculate_energy_test")
    # make input_fn list.
    ap_list = os.listdir(maindir)
    ap_list = [x for x in ap_list if "A0201" in x]
    fn_list = [os.path.join(maindir, x, "builded_template_minim_best_Repair.pdb") for x in ap_list]
    with mp.Pool(30) as pool:
        cal_e_func = calculate_energy3.CalculateEnergy(foldx)
        pool.map(cal_e_func, fn_list)


