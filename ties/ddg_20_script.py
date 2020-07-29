#!/usr/bin/env python3
import os
import time
import pickle as pkl
from pathlib import Path

from ties.ddg20 import extract_energies, merge_datasets, analyse

analysis_dir = 'analysis'
data_dir = 'data'

# make directories if they do not exist
if not os.path.isdir(analysis_dir):
    os.mkdir(analysis_dir)
if not os.path.isdir(data_dir):
    os.mkdir(data_dir)

calc_aga_err = True

protein = 'thrombin'
cases = ['l1_l9', 'l3_l5', 'l5_l6', 'l1_l8', 'l4_l10']

for case in cases:
    ligs = []
    for i in range(1, 4+1):
        lig = extract_energies(Path('/home/dresio/ucl/validation/replica20/') / protein / f'set{i}' / case / 'lig')
        ligs.append(lig)
    lig_all = merge_datasets(ligs)

    laele_int, lavdw_int, ldvdw_int, ldele_int, lig_data = analyse(lig_all, 'lig', calc_aga_err=calc_aga_err, verbose=True)
    lig_delta = laele_int + lavdw_int - ldvdw_int - ldele_int

    complexes = []
    for i in range(1, 4+1):
        complex = extract_energies(Path('/home/dresio/ucl/validation/replica20/') / protein / f'set{i}' / case / 'complex')
        complexes.append(complex)
    complex_all = merge_datasets(complexes)

    caele_int, cavdw_int, cdvdw_int, cdele_int, complex_data = analyse(complex_all, 'complex', calc_aga_err=calc_aga_err, verbose=True)
    complex_delta = caele_int + cavdw_int - cdvdw_int - cdele_int

    # Give the overall results
    print(f'{protein} {case} 20')
    print(f"Delta Delta: {complex_delta - lig_delta:.4f}")
    print (f"Agastya Error {complex_data['sigma_2017'] + lig_data['sigma_2017']:.4f}")