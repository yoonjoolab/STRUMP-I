#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 08:33:58 2021

@author: duaghk

STRUMP-I Energy data parser

"""

# import library
import pandas as pd
from tqdm import tqdm
from pathlib import Path

class ParseSASA:
    def __init__(self) -> None:
        pass
    
    @staticmethod
    def parse_single_row(row: pd.Series, drop_indexes: list = ['Pdb','allele','peptide','p_len','position']) -> pd.Series: 
        row.index = [f"{x}_{row['position']}" if x.endswith("sasa") else x for x in row.index]
        return_row = row.drop(index=drop_indexes)
        return_row.index = [f"{x}" for x in return_row.index]
        return return_row

    def parse_sasa(self, sasa_fn: str):
        sasa = pd.read_csv(sasa_fn)
        info_cols = ['allele','peptide','p_len']
        info_series = sasa.iloc[0]
        info_series = info_series[info_cols]
        sasa_list = [self.parse_single_row(row) for _, row in sasa.iterrows()]
        sasa_list = [info_series] + sasa_list
        return_series = pd.concat(sasa_list)
        return return_series

    def parse_sasa_dir(self, sasa_dir: Path) -> pd.DataFrame:
        print("Checking ", sasa_dir.absolute())
        sasa_fn_list = sasa_dir.absolute().glob("**/SASA_results.csv")
        sasa_df = [self.parse_sasa(x) for x in tqdm(sasa_fn_list)]
        sasa_df = pd.DataFrame(sasa_df)
        print("Number of files: ", len(sasa_df))
        return sasa_df
    
    def __call__(self, strump_i_dir) -> Path:
        strump_i_dir = Path(strump_i_dir)
        output_fn = strump_i_dir.parent.joinpath("parsed_SASA.csv")
        parsed_sasa = self.parse_sasa_dir(strump_i_dir)
        parsed_sasa.to_csv(output_fn, header=True, index=False)
        return output_fn

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir")
    args = parser.parse_args()
    runner = ParseSASA()
    runner(args.input_dir)



