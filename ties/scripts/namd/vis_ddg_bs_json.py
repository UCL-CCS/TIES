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

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt


def plot(transformations, filename, yrange=0.5):
    # sort by transformation
    transformations.sort(key = lambda k: str(k).split('/')[1])

    # group by transformation
    trans_dgs = {}
    for transformation, replica_sets in itertools.groupby(transformations, key=lambda t:str(t).rsplit('_', maxsplit=1)[0]):
        print(f'Next: {transformation}')

        # create keys
        dg_replicas_samples = {i:[] for i in range(1, 20 + 1)} # repsNum: [sample1, sample2, ]

        for repset in replica_sets:
            # load pickled lig and complex data
            with open(repset) as F:
                lig_all = json.load(F)
                # verify that the length is the same for each "replica number"
                sample_num = len(lig_all['1'])
                for next_rep_num, samples in lig_all.items():
                    if sample_num != len(samples):
                        print(f'In file {repset}, replica {next_rep_num} has {len(samples)} samples, '
                              f'but replica 1 has {sample_num}')

                    # merge the data
                    dg_replicas_samples[int(next_rep_num)].extend(samples)

        print(f'Altogether in each set there is {len(dg_replicas_samples[1])}')
        trans_dgs[transformation] = dg_replicas_samples


    plt.figure(figsize=(15, 10))  # , dpi=80) facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 9})

    counter = 1
    for tran, data in trans_dgs.items():
        # plot the boxplots
        plt.subplot(2, 3, counter)
        plt.xlim([0.5, 20 + 0.5])
        plt.title(tran.split('/')[1])

        # give files their own specific name
        tseed = str(time.time()).split('.')[1]

        for number_of_reps, dGs in data.items():
            plt.boxplot(dGs, positions=[number_of_reps, ], showfliers=False)
        plt.ylabel('$\\rm \Delta G $')
        plt.xlabel('$\\rm Replicas~bootstrapped~\#  $')
        plt.xticks(range(0, 20 + 1, 2), range(0, 20 + 1, 2))

        reps1 = data[1]
        reps20 = data[20]
        r1min = np.min(reps1)
        r1max = np.max(reps1)
        # if np.abs(r1min - r1max) > 1:
        #     plt.ylim([r1min, r1max])
        # else:
        #     # set the range to span 1 kcal/mol with the centre of the 20th
        resp20_mean = np.mean(reps20)
        plt.ylim(resp20_mean - yrange, resp20_mean + yrange)

        counter += 1

    # plt.legend()
    plt.tight_layout()
    # plt.suptitle('TYK2 ligand')
    plt.savefig(f'{filename}_y{yrange*2:0.2f}.png', dpi=300)
    # plt.show()

# ligand
transformations = list(Path('.').glob('analysis/lig_l*_l*_*.json'))
plot(transformations, filename='tyk2_lig')

# complex
transformations = list(Path('.').glob('analysis/complex_l*_l*_*.json'))
plot(transformations, filename='tyk2_complex', yrange=1)