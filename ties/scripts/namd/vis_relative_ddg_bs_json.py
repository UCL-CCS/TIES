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
    trans_ddgs_change = {}
    for transformation, replica_sets in itertools.groupby(transformations,
                                                          key=lambda t: str(t).rsplit('_', maxsplit=1)[0]):
        print(f'Next: {transformation}')

        # create keys
        ddg_replicas_samples = {}  # repsNum: [sample1, sample2, ]

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
                    ddg_replicas_samples.setdefault(int_next_rep, [])
                    ddg_replicas_samples[int_next_rep].extend(samples)

        print(f'in each set there is {len(list(ddg_replicas_samples.values())[0])} samples')
        trans_ddgs_change[transformation] = ddg_replicas_samples

    return trans_ddgs_change


def plot_ddgs_relative(exp, trans_dgs, sd_improvements, work_dir, filename, yrange=0.5, plots=[5, 5], experimental=None):
    """
    Plot for each case what the average change in ddG is.
    """

    plt.figure(figsize=(15, 10))  # , dpi=80) facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 9})

    counter = 1
    # for each transformation
    for tran, data in trans_dgs.items():
        # extract the name of the transformation
        _, _, left, right = tran.split('/')[-1].split('_')
        parsed_trans = left.upper() + ' ' + right.upper()

        # plot the boxplots
        plt.subplot(plots[0], plots[1], counter)
        plt.xlim([0.5, 20 + 0.5])
        plt.title(parsed_trans)

        # for each replica
        for number_of_reps, dddGs in data.items():
            # apply absolute value to understand the "change" to ddG regardless of direction
            plt.boxplot(np.abs(dddGs), positions=[number_of_reps, ], showfliers=False)
        plt.ylabel('$\\rm Bootstrapped~\Delta \Delta \Delta G $')
        plt.xlabel('$\\rm Replicas~\# ~/~ 20  $')
        plt.xticks(range(0, 20 + 1, 2), range(0, 20 + 1, 2))

        # plot the exp value
        # exp_value = exp[parsed_trans]
        # plt.axhline(y=exp_value)

        # ensure that the experimental line can be seen
        # ylim_min, ylim_max = plt.axes().get_ylim()
        # plt.ylim([min(exp_value, ylim_min), max(exp_value, ylim_max)])

        counter += 1

    plt.tight_layout()
    plt.savefig(work_dir / f'relative_ddg_{filename}_y{yrange*2:0.2f}.png', dpi=300)
    # plt.show()


def merge_protein_cases(exp, data):
    """
    For the system or protein, take all the transformations t and
    F = merge each | ddGs(t) - exp | for each transformation t

    So subtract the experiment from each transformation, then take the absolute values of that distributions,
    and merge the distributions together
    """
    merged_abs_dists = {5: [], 10: [], 15: [], 20: []}
    for tran, data in trans_dddGs.items():
        # extract the name of the transformation
        _, _, left, right = tran.split('/')[-1].split('_')
        parsed_trans = left.upper() + ' ' + right.upper()

        exp_value = exp[parsed_trans]

        for number_of_reps, ddGs in data.items():
            d = np.abs(np.array(ddGs) - exp_value)
            merged_abs_dists[number_of_reps].extend(d)

    return merged_abs_dists


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


exp_data = get_exp_data()

systems = []
complexes = []
dddgs = {}

root_work = Path('/home/dresio/ucl/ties/ties20/replica20')

# tyk2
# ligand
transformations = list(root_work.glob('tyk2/analysis/direct_ddg_l*_l*_*.json'))
trans_dddGs = combine(transformations)
# plot_ddgs_relative(exp_data['tyk2'], trans_ddGs, sds, root_work, filename='tyk2', plots=[2, 3])
dddgs['tyk2'] = trans_dddGs #  merge_protein_cases(exp_data['tyk2'], trans_ddGs)

# mcl1
# ligand
transformations = list(root_work.glob('mcl1/analysis/direct_ddg_l*_l*_*.json'))
trans_dddGs = combine(transformations)
# plot_ddgs_relative(exp_data['mcl1'], trans_dddGs, sds, root_work, filename='mcl1', yrange=1, plots=[4, 2])
dddgs['mcl1'] = trans_dddGs # merge_protein_cases(exp_data['mcl1'], trans_ddGs)

# thrombin
# lig
transformations = list(root_work.glob('thrombin/analysis/direct_ddg_l*_l*_*.json'))
trans_dddGs = combine(transformations)
# plot_ddgs_relative(exp_data['thrombin'], trans_dddGs, sds, root_work, filename='thrombin', yrange=0.5, plots=[2, 3])
dddgs['thrombin'] = trans_dddGs # merge_protein_cases(exp_data['thrombin'], trans_dddGs)

# ptp1b
# lig
transformations = list(root_work.glob('ptp1b/analysis/direct_ddg_l*_l*_*.json'))
trans_dddGs = combine(transformations)
# plot_ddgs_relative(exp_data['ptp1b'], trans_dddGs, sds, root_work, filename='ptp1b', yrange=3, plots=[2, 3])
dddgs['ptp1b'] = trans_dddGs# merge_protein_cases(exp_data['ptp1b'], trans_dddGs)


# cdk2
# lig
transformations = list(root_work.glob('cdk2/analysis/direct_ddg_l*_l*_*.json'))
trans_dddGs = combine(transformations)
# plot_ddgs_relative(exp_data['cdk2'], trans_dddGs, sds, root_work, filename='cdk2', yrange=0.5, plots=[1, 3])
dddgs['cdk2'] = trans_dddGs # merge_protein_cases(exp_data['cdk2'], trans_dddGs)


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
fig = plt.figure(figsize=(10, 7))
global_ax = fig.add_subplot(111)
plt.rcParams.update({'font.size': 12, 'legend.fontsize': 9})

protein_counter = 1
avg357 = []
for protein, dists in dddgs.items():
    plt.subplot(3, 2, protein_counter)
    plt.title(protein.upper())
    # if protein_counter in {3}:
    #     plt.ylabel(r'$\rm \langle |f(x) - f(x+1)| \rangle_{average} (kcal/mol) $')
    # if protein_counter in (4, 5):
    plt.xlabel('Replica #')

    # add the experimental range?
    plt.ylim([0, 1.6])
    plt.xlim([0.5, 19.5])

    # each case should be plotted separately
    for tran_name, all_reps in dists.items():
        _, _, ligA, ligB = tran_name.split('/')[-1].split('_')
        if ligA.upper() == 'L6' and ligB.upper() == 'L14':
            continue
        # use only average dddgs for each replica
        data = sorted(all_reps.items())
        replica_numbers = [x[0] for x in data]
        dddg_avgs = [np.mean(np.abs(x[1])) for x in data]
        avg357.append([dddg_avgs[2], dddg_avgs[4], dddg_avgs[6]])
        plt.plot(replica_numbers, dddg_avgs, label=f'{ligA.upper()} {ligB.upper()}')

    plt.xticks(range(1, 20, 2), fontsize=10)
    plt.yticks(fontsize=10)

    plt.legend()
    protein_counter += 1


# add one more plot showing the distribution of how much improvement is seen with 3, 5 and 7 replicas,
plt.subplot(3, 2, protein_counter)
just3 = [d[0] for d in avg357]
just5 = [d[1] for d in avg357]
just7 = [d[2] for d in avg357]
alpha = 0.7
bins = np.linspace(0, 0.75, 30)
plt.hist(just3, bins=bins, label='3 + 1 Replicas', alpha=alpha)
plt.hist(just5, bins=bins, label='5 + 1 Replicas', alpha=alpha)
plt.hist(just7, bins=bins, label='7 + 1 Replicas', alpha=alpha)
plt.xticks(np.linspace(0, 0.7, 8), [f'{x:.1f}' for x in np.linspace(0, 0.7, 8)])
plt.xticks(fontsize=10)
plt.tick_params(axis='y', which='both', left=False, labelleft=False)
plt.xlabel(r'average change in $\rm \Delta \Delta G$ (kcal/mol)')
plt.ylabel('P(x)')
plt.legend()

plt.tight_layout()
plt.subplots_adjust(left=0.08)
fig.text(0.01, 0.3, r'$\rm \langle | \Delta\Delta G(r) - \Delta\Delta G(r+1)| \rangle_{mean}  $ (kcal/mol)', rotation='vertical')

plt.savefig(root_work / 'relative_ddgs_byprot.png')
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