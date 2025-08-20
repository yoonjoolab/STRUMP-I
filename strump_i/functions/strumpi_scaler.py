#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 14:59:07 2021

@author: duaghk

Arranging code.
"""
import gzip
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
# pd.set_option("max_column", 999)



# define functions.

class StrumpiScaler:
    def __init__(self, startcol, scale_fn = None):
        
        if scale_fn:
            with gzip.open(scale_fn, 'rb') as f:
                self.scaler = pickle.load(f)
        else:
            pass
        self.startcol = startcol
        pass

    @staticmethod
    def scale_column(col_series):
        # condition-wise scaling.
        # outlier check. 
        # min-max check first.
        if col_series.max() < 1.1 and col_series.min() > -1.1:
            scale_factor = {col_series.name:1}
        elif len(col_series.value_counts()) == 1:
            scale_factor = {col_series.name:1}
        else:
            q_val = [round(x,2) for x in list(np.arange(0,1,0.01))] + [1]
            data_q = col_series.abs().quantile(q_val)
            # gap check.
            if data_q[1]/data_q[0.99] < 3:
                pass
            elif data_q[0.99] == 0:
                pass
            else:
                col_series = col_series[col_series.abs() <= data_q[0.99]]
    
            # now, max normalization
            maxval = col_series.abs().max()
            col_series = col_series/maxval
            scale_factor = {col_series.name:maxval}
    
        return col_series, scale_factor
    
    def fit_scaler(self, df):
        # need to set index to PDB
        collist = df.columns
        returndf = df[collist[0:self.startcol]]
        scale_dict = {}
        for col in tqdm(collist[self.startcol:]):
            tmpseries = df[col]
            tmpseries, tmp_scale = self.scale_column(tmpseries)
            returndf = pd.merge(returndf, tmpseries, how = "outer", 
                                left_index = True, right_index = True)
            scale_dict.update(tmp_scale)
        
        self.scale_dict = scale_dict
        return returndf
        
    def scaling(self, scaling_df):
        # need to set index to PDB
        collist = scaling_df.columns
        returndf = scaling_df[collist].copy()
        # returndf = scaling_df[collist[0:self.startcol]]
        # returndf = scaling_df[collist[self.startcol:]]
        for col in collist[self.startcol:]:
            returndf[col] = returndf[col]/self.scale_dict[col]
        
        return returndf







