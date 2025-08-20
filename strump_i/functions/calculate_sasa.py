#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
    For SASA calculation
'''

# import os
import freesasa
import pandas as pd
from pathlib import Path


class CalSASA:
    '''
        Peptide's SASA calculation module.
        SASA : Surface Area Solvent Accessability.
    '''

    def __init__(self):
        '''
            pdb_fn : full path of pdb files.
            p_chain : peptide chain labeling in PDB.
        '''

    @staticmethod
    def get_sasa_results(pdb_fn, p_chain) -> freesasa.Result:
        structure = freesasa.Structure(str(pdb_fn))
        results = freesasa.calc(structure)
        p_chain_results = results.residueAreas()[p_chain]
        return p_chain_results

    @staticmethod
    def get_single_aa_sasa(pdb: str, allele: str, peptide: str, i: str, aa_area: freesasa.ResidueArea) -> pd.Series:
        return_series = pd.Series(
            data=[
                pdb, allele, peptide, len(peptide), int(i),
                aa_area.total, aa_area.polar, aa_area.apolar,
                aa_area.relativeTotal, aa_area.relativePolar, aa_area.relativeApolar
            ],
            index=[
                'Pdb', 'allele', 'peptide', 'p_len', 'position',
                'abs_total_sasa', 'abs_polar_sasa', 'abs_apolar_sasa',
                'rel_total_sasa', 'rel_polar_sasa', 'rel_apolar_sasa'
            ]
        )
        return return_series

    def make_sasa_df(self, pdb: str, allele: str, peptide: str, p_chain_results: freesasa.Result) -> pd.DataFrame:
        sasa_list = [self.get_single_aa_sasa(pdb, allele, peptide, k, v) for k, v in p_chain_results.items()]
        sasa_df = pd.DataFrame(sasa_list)
        sasa_df = sasa_df.sort_values(by='position')
        return sasa_df

    def __call__(self, pdb_fn: str, p_chain: str):
        pdb_fn = Path(pdb_fn)
        pdb_dir = pdb_fn.parent
        pdb = pdb_dir.name
        allele, peptide = pdb.split("_")
        p_chain = p_chain
        p_chains = self.get_sasa_results(pdb_fn, p_chain)
        sasa_df = self.make_sasa_df(pdb, allele, peptide, p_chains)
        sasa_df.to_csv(Path(pdb_dir, "SASA_results.csv"),
                       header=True, index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_fn", required=True)
    parser.add_argument("--p_chain", default="P")
    args = parser.parse_args()
    sasa = CalSASA()
    sasa(args.pdb_fn, args.p_chain)
