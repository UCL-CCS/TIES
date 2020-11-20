#!/usr/bin/env python3
"""
Visualise the bootstrapped replicas. We have 20 replicas altogether.
Show as a functin of a number of bootstrapped replicas how the values change
"""
import os
from pathlib import Path
from collections import OrderedDict
import glob
import time
import itertools
import random
import sys
import json
import re

import numpy as np
from numpy import genfromtxt
import scipy.stats as st
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

def combine(transformations):
    """
    Combine the different replicas (4 sets of 5)
    """

    # sort by transformation
    transformations.sort(key=lambda k: str(k).split('/')[-1].rsplit('_', maxsplit=1)[0])

    # group by transformation
    trans_dgs = {}
    for transformation, replica_sets in itertools.groupby(transformations,
                                                          key=lambda t: str(t).split('/')[-1].rsplit('_', maxsplit=1)[0]):
        print(f'Next: {transformation}')
        transformation = transformation.split('_', maxsplit=1)[1]

        # create keys
        dg_replicas_samples = {}  # repsNum: [sample1, sample2, ]

        for repset in replica_sets:
            # load pickled lig and complex data
            with open(repset) as F:
                lig_all = json.load(F)
                # verify that the length is the same for each "replica number"
                sample_num = len(lig_all['1'])
                for next_rep_num, samples in lig_all.items():
                    int_next_rep = int(next_rep_num)
                    if sample_num != len(samples):
                        print(f'In file {repset}, replica {next_rep_num} has {len(samples)} samples, '
                              f'but replica 1 has {sample_num}')

                    # merge the data
                    dg_replicas_samples.setdefault(int_next_rep, [])
                    dg_replicas_samples[int_next_rep].extend(samples)

        print(f'in each set there is {len(dg_replicas_samples[1])} samples')
        trans_dgs[transformation] = dg_replicas_samples

    # for this protein, across the ligand transformation,
    # plot the information in improvement as a function of an increase in the number of replicas
    # ie for each case, calculate how much "range" is removed when increasing the replica
    sd_improvements = OrderedDict()
    sd_improvement_range = []
    for system, data in trans_dgs.items():
        sds = []
        for rep_no, dgs in data.items():
            interval95 = st.t.interval(0.95, len(dgs) - 1, loc=np.mean(dgs), scale=st.sem(dgs))
            sd = np.std(dgs)
            # print(f'{system} with {rep_no}, sd {sd}')
            # print(f'{system} with {rep_no}, interval range {np.abs(interval95[0]-interval95[1])}')
            sds.append(sd)
        sd_improvements[system] = np.array(sds)

    return trans_dgs, sd_improvements


def plot_dgs(trans_dgs, sd_improvements, work_dir, filename, yrange=0.5, plots=[5, 5]):
    plt.figure(figsize=(15, 10))  # , dpi=80) facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 9})

    counter = 1
    for tran, data in trans_dgs.items():
        # plot the boxplots
        plt.subplot(plots[0], plots[1], counter)
        plt.xlim([0.5, 20 + 0.5])
        plt.title(tran.rsplit('/')[-1])

        # give files their own specific name
        tseed = str(time.time()).split('.')[1]

        for number_of_reps, dGs in data.items():
            plt.boxplot(dGs, positions=[number_of_reps, ], showfliers=False)
        plt.ylabel('$\\rm Bootstrapped~\Delta G $')
        plt.xlabel('$\\rm Replicas~\# ~/~ 20  $')
        plt.xticks(range(0, 20 + 1, 2), range(0, 20 + 1, 2))

        reps_largest = data[max(data.keys())]
        reps_largest_mean = np.mean(reps_largest)
        # plt.ylim(reps_largest_mean - yrange, reps_largest_mean + yrange)

        counter += 1

    # plt.legend()
    plt.tight_layout()
    # plt.suptitle('TYK2 ligand')
    plt.savefig(work_dir / f'{filename}_y{yrange*2:0.2f}.png', dpi=300)
    # plt.show()

    # for each system plot how much it improves
    plt.figure(figsize=(15, 10))
    plt.rcParams.update({'font.size': 9})

    counter = 1
    for system, sds in sd_improvements.items():
        plt.subplot(plots[0], plots[1], counter)
        # plt.xlim([0.5, 20 + 0.5])
        plt.title(system.rsplit('/')[-1])
        plt.plot(sds)
        plt.tight_layout()
        plt.xlim([-0.5, 20.5])
        plt.xlabel('Number of replicas')
        plt.ylabel('$\\rm  \Delta G ~ \sigma $')
        plt.savefig(work_dir / f'sds_decrease_{filename}.png', dpi=300)
        counter += 1

def plot_sds(sdsL, sdsC, work_dir, filename, plots=[5, 5]):
    # for each system plot how much it improves
    plt.figure(figsize=(15, 10))
    plt.rcParams.update({'font.size': 9})

    counter = 1
    for system, sds1 in sdsL.items():
        plt.subplot(plots[0], plots[1], counter)
        # plt.xlim([0.5, 20 + 0.5])
        plt.ylim([0, 3])
        plt.title(system.rsplit('/')[-1])
        plt.plot(sds1, label='In Water')
        sds2 = sdsC[system]
        plt.plot(sds2, label='In Protein')
        plt.legend()
        plt.tight_layout()
        plt.xlim([-0.5, 20.5])
        plt.xlabel('Number of replicas')
        plt.ylabel('$\\rm  \Delta G ~ \sigma $')
        plt.savefig(work_dir / f'sds_decrease_{filename}.png', dpi=300)
        counter += 1


def get_exp_data():
    # load the data
    # columns 0 are names, and 13 are experiments
    data = genfromtxt("/home/dresio/pubs/ties20/figures/energies.csv",
                      delimiter=',', dtype=None, encoding='UTF-8', usecols=[0,13])

    # separate the proteins
    exp_data = {}
    tyk = data[4:15].T
    exp_data['tyk2'] = {t: float(ddg) for t, ddg in zip(tyk[0], tyk[1])}
    mcl = data[16:32].T
    exp_data['mcl1'] = {t: float(ddg) for t, ddg in zip(mcl[0], mcl[1])}
    thrombin = data[34:44].T
    exp_data['thrombin'] = {t: float(ddg) for t, ddg in zip(thrombin[0], thrombin[1])}
    ptp1b = data[46:55].T
    exp_data['ptp1b'] = {t: float(ddg) for t, ddg in zip(ptp1b[0], ptp1b[1])}
    cdk = data[57:63].T
    exp_data['cdk2'] = {t: float(ddg) for t, ddg in zip(cdk[0], cdk[1])}

    return exp_data

def extract_every5_merge_cases(trans_dgs, exp_protein):
    """
    Take the distribution of dG across the different transformations,
    and subtract the exp value.
    """
    fives = []
    tens = []
    fifteens = []
    twenties = []
    for sys, data in trans_dgs.items():
        if 'l6_l14' in sys:
            continue
        # look up the experimental value for this case
        case_parts = sys.split('/')[-1].split('_')
        case = case_parts[1].upper() + ' ' + case_parts[2].upper()
        exp = exp_protein[case]
        fives.extend(np.array(data[5]) - exp)
        tens.extend(np.array(data[10]) - exp)
        fifteens.extend(np.array(data[15]) - exp)
        twenties.extend(np.array(data[max(data.keys())]) - exp)
        # extract
    return fives, tens, fifteens, twenties


exp_data = get_exp_data()
ligs = OrderedDict()
complexs = OrderedDict()

root_work = Path('/home/dresio/ucl/ties/ties20/replica20')

# tyk2
# ligand
transformations = list(root_work.glob('tyk2/analysis/lig_l*_l*_*.json'))
dgs, sdsL = combine(transformations)
# plot_dgs(dgs, sds, root_work, filename='tyk2_lig', plots=[2,3])
ligs['tyk2'] = sdsL
# complex
transformations = list(root_work.glob('tyk2/analysis/complex_l*_l*_*.json'))
dgs, sdsC = combine(transformations)
# plot_dgs(dgs, sds, root_work, filename='tyk2_complex',yrange=1, plots=[2,3])
# plot_sds(sdsL, sdsC, root_work, filename='tyk2', plots=[2,3])
complexs['tyk2'] = sdsC


# mcl1
# ligand
transformations = list(root_work.glob('mcl1/analysis/lig_l*_l*_*.json'))
dgs, sdsL = combine(transformations)
# plot_dgs(dgs, sds, root_work, filename='mcl1_lig', yrange=1, plots=[4, 2])
ligs['mcl1'] = sdsL
# complex
transformations = list(root_work.glob('mcl1/analysis/complex_l*_l*_*.json'))
dgs, sdsC = combine(transformations)
# plot_dgs(dgs, sds, root_work, filename='mcl1_complex', yrange=2.5, plots=[4,2])
# plot_sds(sdsL, sdsC, root_work, filename='mcl1', plots=[4,2])
complexs['mcl1'] = sdsC

# thrombin
# lig
transformations = list(root_work.glob('thrombin/analysis/lig_l*_l*_*.json'))
dgs, sds = combine(transformations) # , root_work, filename='thrombin_lig', yrange=0.5, plots=[2, 3])
ligs['thrombin'] = sds
# complex
transformations = list(root_work.glob('thrombin/analysis/complex_l*_l*_*.json'))
dgs, sds = combine(transformations) # , root_work, filename='thrombin_complex', yrange=2.5, plots=[2,3])
complexs['thrombin'] = sds

# ptp1b
# lig
transformations = list(root_work.glob('ptp1b/analysis/lig_l*_l*_*.json'))
dgs, sds = combine(transformations)
ligs['ptp1b'] = sds
# plot_dgs(trans_dgs, ptp1b_lig_sds, root_work, filename='ptp1b_lig', yrange=3, plots=[2, 3])
# complex
transformations = list(root_work.glob('ptp1b/analysis/complex_l*_l*_*.json'))
dgs, sds = combine(transformations)
# plot_dgs(trans_dgs, ptp1b_complex_sds, root_work, filename='ptp1b_complex', yrange=5, plots=[2, 3])
complexs['ptp1b'] = sds

# cdk2
# lig
transformations = list(root_work.glob('cdk2/analysis/lig_l*_l*_*.json'))
dgs, sds = combine(transformations) #, root_work, filename='cdk2_lig', yrange=0.5, plots=[1, 3])
ligs['cdk2'] = sds
# complex
transformations = list(root_work.glob('cdk2/analysis/complex_l*_l*_*.json'))
dgs, sds = combine(transformations) # , root_work, filename='cdk2_complex', yrange=1, plots=[1,3])
complexs['cdk2'] = sds

symbol_mapping = {
    'tyk2': 'x',
    'mcl1': '.',
    'thrombin': '*',
    'ptp1b': '+',
    'cdk2': 'o',
}
colour_mapping = {
    'tyk2': '#E13D77',
    'mcl1': '#C9E13D',
    'thrombin': '#3DE1A7',
    'ptp1b': '#553DE1',
    'cdk2': 'black',
}

# this was used before
color = {
    'tyk2' : 'red',
    'mcl1' : 'blue',
    'thrombin' : 'green',
    'ptp1b' : 'orange',
    'cdk2' : 'purple'
}

# Reduce the data type to take the means of the bootstrapping.
# So for each protein, each transformation, LA and LB take the data
# which is the dG bootstrapped as a function of the number of replicas,
# and do RMSD between LA and LB to find the overall difference,
# then plot it as part of the heatmap

# count transformations
transformation_count = sum(len(v) for k, v in ligs.items())
# create a heatmap
heatmap = np.zeros([transformation_count, transformation_count])

# each transformation has a unique plotting column/row, which will be the same
header_mapping = {}
counter = 0
for prot, transformations in ligs.items():
    for tran1 in transformations.keys():
        header_mapping[prot + tran1] = counter
        counter += 1

plt.figure(figsize=(5, 5))
plt.rcParams.update({'font.size': 12})

com_sds_one = []
lig_sds_one = []
for prot, transformations in ligs.items():
    # each trans to each other
    complex_trans = complexs[prot]
    label_kwarg = {'label': prot.upper()}  # use this only once
    for tran1, lig_sds in transformations.items():
        if prot == 'ptp1b' and tran1 == 'l6_l14':
            continue
        complex_sds = complex_trans[tran1]
        com_sds_one.append(complex_sds[0])
        lig_sds_one.append(lig_sds[0])
        print(complex_sds[0], lig_sds[0])
        plt.scatter(lig_sds[0], complex_sds[0], color=color[prot], **label_kwarg)
        label_kwarg['label'] = ''
corr_sds = pearsonr(com_sds_one, lig_sds_one)
print('corr', corr_sds)
plt.ylabel(r'$\rm SD(\Delta G_{alch}^{bound}) (kcal/mol)  $')
plt.xlabel(r'$\rm SD(\Delta G_{alch}^{aq}) (kcal/mol) $')
plt.text(0.01, 2.9, f'$\\rm \\rho$ = {corr_sds[0]:.2f}')
plt.legend(loc='lower right')
plt.savefig(root_work / 'dg_sds_corr_lig_complex.png', dpi=300)
# plt.show()
# ---------------------------------------------------------------
# combine all sds, so show the points and the distribution
# ligands first
# plt.figure(figsize=(15, 10))
# plt.rcParams.update({'font.size': 9})
# plt.subplot(1,2,1)
# plt.title('Ligands')
# plt.ylabel('$\\rm  \Delta G ~ \sigma $')
# plt.xlabel('Replica number')
#
# for i in range(20):
#     # extract from each case just the one index
#     for lig in ligands:
#         for system, sds in lig.items():
#             if 'l6_l14' in system:
#                 continue
#             name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
#             sym = symbol_mapping[name]
#             color = colour_mapping[name]
#             plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color)
#
# # complex
# plt.subplot(1, 2, 2)
# plt.title('Complexes')
# plt.ylabel('$\\rm  \Delta G ~ \sigma $')
# plt.xlabel('Replica number')
#
# for i in range(20):
#     # extract from each case just the one index
#     for complex in complexes:
#         for system, sds in complex.items():
#             if 'l6_l14' in system:
#                 continue
#             name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
#             sym = symbol_mapping[name]
#             color = colour_mapping[name]
#             plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color, label=name)
#
# # remove duplicate labels
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys())
#
# plt.tight_layout()
# plt.savefig(root_work / 'combined_sds.png')

# ---------------------------------------------------------------
# now we do the same thing but use the experimental data instead. We take the distribution, and subtract the exp value
# as a function of 5, 10, 15, 20 replicas,
# so we'll see if the values come closer to the experimental value, rather than just decreasing error

# combine all sds, so show the points and the distribution
# ligands first
# plt.figure(figsize=(15, 10))
# plt.rcParams.update({'font.size': 9})
# plt.subplot(1,2,1)
# plt.title('Ligands')
# plt.ylabel('$\\rm  \Delta G ~ \sigma $')
# plt.xlabel('Replica number')
#
# for i in range(20):
#     # extract from each case just the one index
#     for lig in ligands:
#         for system, sds in lig.items():
#             if 'l6_l14' in system:
#                 continue
#             name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
#             sym = symbol_mapping[name]
#             color = colour_mapping[name]
#             plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color)

# complex
# plt.subplot(1, 2, 2)
# plt.title('Complexes')
# plt.ylabel('$\\rm  \Delta G ~ \sigma $')
# plt.xlabel('Replica number')
#
# for i in range(20):
#     # extract from each case just the one index
#     for complex in complexes:
#         for system, sds in complex.items():
#             if 'l6_l14' in system:
#                 continue
#             name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
#             sym = symbol_mapping[name]
#             color = colour_mapping[name]
#             plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color, label=name)
#
# # remove duplicate labels
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys())
#
# plt.tight_layout()
# plt.savefig(root_work / 'combined_sds.png')