# -*- coding: utf-8 -*-

import os.path

import numpy as np
import pandas as pd

import dGinsertion

def dg_score_compare(old_output_fn, **kwargs):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    old_output_fn = os.path.join(dir_path, old_output_fn)

    calc = pd.read_csv(old_output_fn, sep='\t', names=['seq', 'olddG'])

    calc = pd.concat([calc,
        calc['seq'].str.replace('.', '').
                    apply(dGinsertion.scan_for_best_TM_dGraw, **kwargs).
                    apply(pd.Series, index=['newdG', 'start', 'length'])],
                                    axis=1)

    # skipping cases where old script produced an error
    calc['olddG'] = pd.to_numeric(calc['olddG'], errors='coerce')
    calc = calc.dropna(subset=['olddG'])

    # rtol gives issues when atol is satisfied but the abs(number) is small
    assert np.allclose(calc['olddG'].values,
                       calc['newdG'].values, atol=.001, rtol=1)

    return


def test_calcdG_local_compare():
    # with_length on by default
    dg_score_compare("Daley_gfp.segs.dGpred_local", with_length=True)
    return

def test_calcdG_web_compare():
    dg_score_compare("Daley_gfp.segs.dGpred_web")
    return

def test_calcdG_web_lencorr_compare():
    dg_score_compare("Daley_gfp.segs.dGpred_web_lencorr", with_length=True)
    return

def test_calcdG_web_subseq_compare():
    dg_score_compare("Daley_gfp.segs.dGpred_web_subseq", allow_sub=True)
    return

def test_calcdG_web_lencorr_subseq_compare():
    dg_score_compare("Daley_gfp.segs.dGpred_web_lencorr_subseq",
                     allow_sub=True, with_length=True)
    return

