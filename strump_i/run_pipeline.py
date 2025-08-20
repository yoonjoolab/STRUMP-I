#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 18:05:25 2021

@author: duaghk

STRUMP-I multi peptides.

"""

# import library
import os
from pathlib import Path
import pandas as pd
import argparse
import multiprocessing as mp
from itertools import repeat

import strumpi
# from functions import parse_energy, parse_quality, parse_sasa, predict_binder, calculate_sasa
from strump_i.functions import parse_energy, parse_quality


# define functions
class STRUMP_I:
    def __init__(self, ref_dir: Path) -> None:
        self.ref_dir = ref_dir
        pass

    def run_strumpi_modeler(self, args):
        '''
        run_strumpi_binder 

        STRUMP-I PDB structure make and predict whether binder or not.

        Args:
            args (parser): multi args parser.
                allele, pepfile, pdbdir, outdir, tinker_dir, foldx
        '''
        # peptide list parsing
        run_strumpi = strumpi.STRUMPI(
            args.pdb_dir, args.output_dir, args.tinker_dir, args.foldx,
            args.mhc_seq_mat_fn, args.pdb_seq_mat_fn, args.templateinfo_fn,
            args.blosum_fn, args.pam_fn, args.keyfile_fn, args.rms
            )
        with open(args.pepfile) as f_read:
            peplist = [x.strip() for x in f_read.readlines()]
        with mp.Pool(args.threads) as pool:
            pool.starmap(run_strumpi, zip(repeat(args.allele), peplist))
    
    def run_strumpi_binding_prdictor(self, args):
    	pass
#         q_parser = parse_quality.ParseQuality(args.input_dir, Path(args.output_dir))
#         quality_fn = q_parser()
#         e_parser = parse_energy.ParseEnergy(args.input_dir, Path(args.output_dir))
#         energy_fn = e_parser()
#         binding_predictor = predict_binder.StrumpiBindingClassifier(args.model_dir)
#         binding_prediction_results = binding_predictor(args.output_dir, args.output_dir)
#         print(f"Final result path: {binding_prediction_results}")
    
    def run_strumpi_sasa_calculator(self, args):
    	pass
#         sasa_dir = Path(args.input_dir)
#         pdb_list = sasa_dir.glob("*/*/builded_template_minim_best_Repair.pdb")
#         sasa_calculator = calculate_sasa.CalSASA()
#         with mp.Pool(args.threads) as pool:
#             pool.starmap(sasa_calculator, zip(pdb_list, repeat("P")))
#         sasa_parser = parse_sasa.ParseSASA()
#         sasa_parser(args.input_dir)

    def run_strumpi_immunogen(self, args):
        pass


if __name__ == "__main__":
    main_dir = Path(__file__).absolute().parent
    ref_dir = Path(main_dir, "reference")
    runner = STRUMP_I(ref_dir)
    parser = argparse.ArgumentParser(prog='STRUMP-I')
    subparser = parser.add_subparsers()
    subparser_modeler = subparser.add_parser(name="Modeler")
    subparser_modeler.add_argument("--allele", required=True)
    subparser_modeler.add_argument("--pepfile", required=True)
    subparser_modeler.add_argument("--pdb_dir", default=ref_dir.joinpath("pdb"))
    subparser_modeler.add_argument("--output_dir", required=True)
    subparser_modeler.add_argument("--tinker_dir", default=os.path.dirname("/app/tinker"))
    subparser_modeler.add_argument("--foldx", default=os.path.dirname("/app/foldx/foldx"))
    subparser_modeler.add_argument("--mhc_seq_mat_fn", default=ref_dir.joinpath("HLA_template.pickle"))
    subparser_modeler.add_argument("--pdb_seq_mat_fn", default=ref_dir.joinpath("PDB_template.pickle"))
    subparser_modeler.add_argument("--templateinfo_fn", default=ref_dir.joinpath("IMGT_structure_summary.csv"))
    subparser_modeler.add_argument("--blosum_fn", default=ref_dir.joinpath("alignment_matrix", "BLOSUM62"))
    subparser_modeler.add_argument("--pam_fn", default=ref_dir.joinpath("alignment_matrix", "PAM30"))
    subparser_modeler.add_argument("--keyfile_fn", default=ref_dir.joinpath("force.key"))
    subparser_modeler.add_argument("--rms", default=0.01)
    subparser_modeler.add_argument("--threads", default=4, type=int)
    subparser_modeler.set_defaults(func=runner.run_strumpi_modeler)

    subparser_binding_predictor = subparser.add_parser("BindingPredictor")
    subparser_binding_predictor.add_argument("--input_dir")
    subparser_binding_predictor.add_argument("--output_dir")
    subparser_binding_predictor.add_argument("--model_dir", default=ref_dir.joinpath("model", "binding"))
    subparser_binding_predictor.set_defaults(func=runner.run_strumpi_binding_prdictor) 

    subparser_sasa_calculator = subparser.add_parser("SASACalculator")
    subparser_sasa_calculator.add_argument("--input_dir")
    subparser_sasa_calculator.add_argument("--threads", type=int)
    subparser_sasa_calculator.set_defaults(func=runner.run_strumpi_sasa_calculator) 

    args = parser.parse_args()
    print(args)
    args.func(args)