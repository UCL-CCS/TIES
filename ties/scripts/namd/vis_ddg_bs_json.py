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
import scipy.stats as st
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

def combine(transformations):
    """
    Combine the different replicas (4 sets of 5)
    """

    # sort by transformation
    transformations.sort(key=lambda k: str(k).split('/')[1])

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
            plt.boxplot(dGs, positions=[number_of_reps, ], showfliers=True)
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


ligands = []
complexes = []

root_work = Path('/home/dresio/ucl/validation/replica20')

# tyk2
# ligand
# transformations = list(root_work.glob('tyk2/analysis/lig_l*_l*_*.json'))
# _, tyk_lig_sds = combine_and_plot(transformations, root_work, filename='tyk2_lig', plots=[2,3])
# ligands.append(tyk_lig_sds)
# # complex
# transformations = list(root_work.glob('tyk2/analysis/complex_l*_l*_*.json'))
# _, tyk_complex_sds = combine_and_plot(transformations, root_work, filename='tyk2_complex', yrange=1, plots=[2,3])
# complexes.append(tyk_complex_sds)


# mcl1
# ligand
# transformations = list(root_work.glob('mcl1/analysis/lig_l*_l*_*.json'))
# _, mcl_lig_sds = combine_and_plot(transformations, root_work, filename='mcl1_lig', yrange=1, plots=[4, 2])
# ligands.append(mcl_lig_sds)
# # complex
# transformations = list(root_work.glob('mcl1/analysis/complex_l*_l*_*.json'))
# _, mcl_complex_sds = combine_and_plot(transformations, root_work, filename='mcl1_complex', yrange=2.5, plots=[4,2])
# complexes.append(mcl_complex_sds)
#
#
# thrombin
# lig
# transformations = list(root_work.glob('thrombin/analysis/lig_l*_l*_*.json'))
# _, thrombin_lig_sds = combine_and_plot(transformations, root_work, filename='thrombin_lig', yrange=0.5, plots=[2, 3])
# ligands.append(thrombin_lig_sds)
# # complex
# transformations = list(root_work.glob('thrombin/analysis/complex_l*_l*_*.json'))
# _, thrombin_complex_sds = combine_and_plot(transformations, root_work, filename='thrombin_complex', yrange=2.5, plots=[2,3])
# complexes.append(thrombin_complex_sds)
#
#
# ptp1b
# lig
transformations = list(root_work.glob('ptp1b/analysis/lig_l*_l*_*.json'))
trans_dgs, ptp1b_lig_sds = combine(transformations)
plot_dgs(trans_dgs, ptp1b_lig_sds, root_work, filename='ptp1b_lig', yrange=3, plots=[2, 3])
ligands.append(ptp1b_lig_sds)
# complex
transformations = list(root_work.glob('ptp1b/analysis/complex_l*_l*_*.json'))
trans_dgs, ptp1b_complex_sds = combine(transformations)
plot_dgs(trans_dgs, ptp1b_complex_sds, root_work, filename='ptp1b_complex', yrange=5, plots=[2, 3])
complexes.append(ptp1b_complex_sds)
#
#
# cdk2
# lig
# transformations = list(root_work.glob('cdk2/analysis/lig_l*_l*_*.json'))
# _, cdk_lig_sds = combine_and_plot(transformations, root_work, filename='cdk2_lig', yrange=0.5, plots=[1, 3])
# ligands.append(cdk_lig_sds)
#
# # complex
# transformations = list(root_work.glob('cdk2/analysis/complex_l*_l*_*.json'))
# _, cdk_complex_sds = combine_and_plot(transformations, root_work, filename='cdk2_complex', yrange=1, plots=[1,3])
# complexes.append(cdk_complex_sds)


# combine all sds, so show the points and the distribution
# ligands first
plt.figure(figsize=(15, 10))
plt.rcParams.update({'font.size': 9})
plt.subplot(1,2,1)
plt.title('Ligands')
plt.ylabel('$\\rm  \Delta G ~ \sigma $')
plt.xlabel('Replica number')

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

for i in range(20):
    # extract from each case just the one index
    for lig in ligands:
        for system, sds in lig.items():
            name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
            sym = symbol_mapping[name]
            color = colour_mapping[name]
            plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color)

# complex
plt.subplot(1, 2, 2)
plt.title('Complexes')
plt.ylabel('$\\rm  \Delta G ~ \sigma $')
plt.xlabel('Replica number')

for i in range(20):
    # extract from each case just the one index
    for complex in complexes:
        for system, sds in complex.items():
            name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
            sym = symbol_mapping[name]
            color = colour_mapping[name]
            plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color, label=name)

# remove duplicate labels
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

plt.tight_layout()
plt.savefig(root_work / 'combined_sds.png')