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

def combine(transformations):
    """
    Combine the different replicas (4 sets of 5)
    """

    # sort by transformation
    transformations.sort(key=lambda k: str(k).rsplit('/', maxsplit=1)[1])

    # group by transformation
    trans_dgs = {}
    for transformation, replica_sets in itertools.groupby(transformations,
                                                          key=lambda t: str(t).rsplit('_', maxsplit=1)[0]):
        print(f'Next: {transformation}')

        # create keys
        dg_replicas_samples = {}  # repsNum: [sample1, sample2, ]

        for repset in replica_sets:
            # load pickled lig and complex data
            with open(repset) as F:
                lig_all = json.load(F)
                # verify that the length is the same for each "replica number"
                sample_num = len(list(lig_all.values())[0])
                for next_rep_num, samples in lig_all.items():
                    int_next_rep = int(next_rep_num)
                    if sample_num != len(samples):
                        print(f'In file {repset}, replica {next_rep_num} has {len(samples)} samples, '
                              f'but replica 1 has {sample_num}')

                    # merge the data
                    dg_replicas_samples.setdefault(int_next_rep, [])
                    dg_replicas_samples[int_next_rep].extend(samples)

        print(f'in each set there is {len(list(dg_replicas_samples.values())[0])} samples')
        trans_dgs[transformation] = dg_replicas_samples

    # for this protein, across the ligand transformation,
    # plot the information in improvement as a function of an increase in the number of replicas
    # ie for each case, calculate how much "range" is removed when increasing the replica
    sd_improvements = {}
    sd_improvement_range = []
    for system, data in trans_dgs.items():
        sds = []
        for rep_no, dgs in data.items():
            interval95 = st.t.interval(0.95, len(dgs) - 1, loc=np.mean(dgs), scale=st.sem(dgs))
            sd = np.std(dgs)
            # print(f'{system} with {rep_no}, sd {sd}')
            # print(f'{system} with {rep_no}, interval range {np.abs(interval95[0]-interval95[1])}')
            sds.append(sd)
        sd_improvements[system] = sds

    return trans_dgs, sd_improvements


def plot_ddgs(exp, trans_dgs, sd_improvements, work_dir, filename, yrange=0.5, plots=[5, 5], experimental=None):
    plt.figure(figsize=(15, 10))  # , dpi=80) facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 9})

    counter = 1
    for tran, data in trans_dgs.items():
        # extract the name of the transformation
        _, left, right = tran.split('/')[-1].split('_')
        parsed_trans = left.upper() + ' ' + right.upper()

        # plot the boxplots
        plt.subplot(plots[0], plots[1], counter)
        plt.xlim([0.5, 20 + 0.5])
        plt.title(parsed_trans)

        for number_of_reps, ddGs in data.items():
            plt.boxplot(ddGs, positions=[number_of_reps, ], showfliers=False)
        plt.ylabel('$\\rm Bootstrapped~\Delta \Delta G $')
        plt.xlabel('$\\rm Replicas~\# ~/~ 20  $')
        plt.xticks(range(0, 20 + 1, 2), range(0, 20 + 1, 2))

        # plot the exp value
        exp_value = exp[parsed_trans]
        plt.axhline(y=exp_value)

        # ensure that the experimental line can be seen
        # ylim_min, ylim_max = plt.axes().get_ylim()
        # plt.ylim([min(exp_value, ylim_min), max(exp_value, ylim_max)])

        counter += 1

    plt.tight_layout()
    plt.savefig(work_dir / f'ddg_{filename}_y{yrange*2:0.2f}.png', dpi=300)
    # plt.show()

    # ----------------------------------
    # for each system plot show how much the system decreases the uncertainty
    plt.figure(figsize=(15, 10))
    plt.rcParams.update({'font.size': 9})

    counter = 1
    for system, sds in sd_improvements.items():
        plt.subplot(plots[0], plots[1], counter)
        # plt.xlim([0.5, 20 + 0.5])
        plt.title(system.rsplit('/')[-1])
        plt.plot([5, 10, 15, 20], sds)
        plt.tight_layout()
        plt.xlim([-0.5, 20.5])
        plt.xlabel('Number of replicas')
        plt.ylabel('$\\rm  \Delta\Delta G ~ \sigma $')
        plt.savefig(work_dir / f'ddg_sds_decrease_{filename}.png', dpi=300)
        counter += 1

def merge_protein_cases(exp, data):
    """
    For the system or protein, take all the transformations t and
    F = merge each | ddGs(t) - exp | for each transformation t

    So subtract the experiment from each transformation, then take the absolute values of that distributions,
    and merge the distributions together
    """
    merged_abs_dists = {5: [], 10: [], 15: [], 20: []}
    for tran, data in trans_ddGs.items():
        # extract the name of the transformation
        _, left, right = tran.split('/')[-1].split('_')
        parsed_trans = left.upper() + ' ' + right.upper()

        exp_value = exp[parsed_trans]

        for number_of_reps, ddGs in data.items():
            d = np.abs(np.array(ddGs) - exp_value)
            merged_abs_dists[number_of_reps].extend(d)

    return merged_abs_dists


def get_exp_data():
    # load the data
    # columns 0 are names, and 13 are experiments
    data = genfromtxt("/home/dresio/pubs/ties_replication_validation/figures/energies.csv",
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


exp_data = get_exp_data()

systems = []
complexes = []
ddg_minus_exp = {}

root_work = Path('/home/dresio/ucl/validation/replica20')

# tyk2
# ligand
transformations = list(root_work.glob('tyk2/analysis/ddg_l*_l*_*.json'))
trans_ddGs, sds = combine(transformations)
plot_ddgs(exp_data['tyk2'], trans_ddGs, sds, root_work, filename='tyk2', plots=[2, 3])
ddg_minus_exp['tyk2'] = merge_protein_cases(exp_data['tyk2'], trans_ddGs)

# mcl1
# ligand
transformations = list(root_work.glob('mcl1/analysis/ddg_l*_l*_*.json'))
trans_ddGs, sds = combine(transformations)
plot_ddgs(exp_data['mcl1'], trans_ddGs, sds, root_work, filename='mcl1', yrange=1, plots=[4, 2])
ddg_minus_exp['mcl1'] = merge_protein_cases(exp_data['mcl1'], trans_ddGs)

# thrombin
# lig
transformations = list(root_work.glob('thrombin/analysis/ddg_l*_l*_*.json'))
trans_ddGs, sds = combine(transformations)
plot_ddgs(exp_data['thrombin'], trans_ddGs, sds, root_work, filename='thrombin', yrange=0.5, plots=[2, 3])
ddg_minus_exp['thrombin'] = merge_protein_cases(exp_data['thrombin'], trans_ddGs)

# ptp1b
# lig
transformations = list(root_work.glob('ptp1b/analysis/ddg_l*_l*_*.json'))
trans_ddGs, sds = combine(transformations)
plot_ddgs(exp_data['ptp1b'], trans_ddGs, sds, root_work, filename='ptp1b', yrange=3, plots=[2, 3])
ddg_minus_exp['ptp1b'] = merge_protein_cases(exp_data['ptp1b'], trans_ddGs)


# cdk2
# lig
transformations = list(root_work.glob('cdk2/analysis/ddg_l*_l*_*.json'))
trans_ddGs, sds = combine(transformations)
plot_ddgs(exp_data['cdk2'], trans_ddGs, sds, root_work, filename='cdk2', yrange=0.5, plots=[1, 3])
ddg_minus_exp['cdk2'] = merge_protein_cases(exp_data['cdk2'], trans_ddGs)


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


# ------------------------
# plot for each protein separately the progress
plt.figure(figsize=(15, 10))
plt.rcParams.update({'font.size': 9})

protein_counter = 1
for protein, dists in ddg_minus_exp.items():
    plt.subplot(2, 3, protein_counter)
    plt.title(protein.upper())
    plt.ylabel('$\\rm |\Delta\Delta G - exp| $')
    plt.xlabel('# Replicas used in bootstrapping / 20')
    plt.ylim([-0.1, 3])

    case_counter = 1
    for rep_no, ddgs in dists.items():
        plt.boxplot(ddgs, positions=[case_counter, ], widths=0.13, showfliers=False)
        case_counter += 1

    plt.xticks([1, 2, 3, 4], [5, 10, 15, 20])


    protein_counter += 1

plt.tight_layout()
plt.savefig(root_work / 'ddgs_byprot.png')
# plt.show()

# -----------------------------------
# lot the ddg-exp joined distributions as a bar plot
# plt.figure(figsize=(15, 10))
# plt.rcParams.update({'font.size': 9})
# plt.ylabel('$\\rm Bootstrapped~ |\Delta\Delta G - exp| $')
# plt.xlabel('Replica number')
#
# # rearrange the data, such that each internal list represents different proteins
#
# protein_factor = -0.3
# for protein, dists in ddg_minus_exp.items():
#     case_counter = 1
#     for rep_no, ddgs in dists.items():
#         plt.boxplot(ddgs, positions=[case_counter + protein_factor, ], widths=0.13, showfliers=False)
#         case_counter += 1
#
#
#     protein_factor += 0.15
#
# plt.xticks([1, 2, 3, 4], [5, 10, 15, 20])
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
# plt.title('ddG')
# plt.ylabel('$\\rm  \Delta G ~ \sigma $')
# plt.xlabel('Replica number')
#
# for i in range(20):
#     # extract from each case just the one index
#     for lig in systems:
#         for system, sds in lig.items():
#             if 'l6_l14' in system:
#                 continue
#             name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
#             sym = symbol_mapping[name]
#             color = colour_mapping[name]
#             plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color)
#
#
# plt.tight_layout()
# plt.savefig(root_work / 'combined_sds.png')