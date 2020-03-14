"""
# get from each lambda/replica the energies, sum the average vdw and avg ele
# todo - upgrade to pathlib
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import sem
import pandas as pd
from itertools import accumulate
from pymbar import timeseries
from collections import OrderedDict


def choder_get_eqpart(datapoints):
    [t0, g, Neff_max] = timeseries.detectEquilibration(datapoints)
    return datapoints[t0:]


def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[int(result.size / 2):]


def populate_stats(dataset):
    meta = {}
    for lambda_val, reps in dataset.items():
        if f'merged_mean' not in meta:
            meta[f'merged_mean'] = {}

        # flatten the list
        merged_replicas = []
        for rep in reps:
            # merged_replicas.extend(rep)
            # break
            # use Chodera's EQ code to only extract the prod values after equilibration
            merged_replicas.extend(choder_get_eqpart(rep))
        # item[f'merged_{itype}_sem'] = sem(merged_replicas)
        meta[f'merged_mean'][lambda_val] = np.mean(merged_replicas)
        # item[f'merged_{itype}_std'] = np.std(merged_replicas)
        # print('sem is overall', item['merged_avdw_sem'])
    return meta


def get_energies(location):
    # take only the lambda directories
    lambda_dirs = filter(lambda d: d.is_dir(), Path(location).glob(r'lambda_[0-1].[0-9]*'))
    # sort lambda directories in the increasing order
    lambda_dirs = sorted(lambda_dirs, key=lambda d: float(d.name.split('_')[1]))

    # different datasets
    data = {'avdw': OrderedDict(), 'dvdw': OrderedDict(), 'aele': OrderedDict(), 'dele': OrderedDict()}

    for lambda_dir in lambda_dirs:
        for rep in lambda_dir.glob('rep[0-9]*'):
            if not rep.is_dir():
                continue

            prod_alch = rep / 'prod.alch'

            # partition1 is appearing
            # partition2 is disappearing

            # COLUMNS:
            # 5 is AVGELECT1
            # 7 is AVGVDW1
            # 11 is AVGELECT2
            # 13 is AVGVDW2
            # we take the last row because these are the running averages computed by NAMD
            energies = np.loadtxt(prod_alch, comments='#', usecols=[5, 7, 11, 13])[-1]
            energies_cum_roll = np.loadtxt(prod_alch, comments='#', usecols=[5, 7, 11, 13]).T

            # extract the raw datapoints
            # 4 is ELECT1
            # 6 is VDW1
            # 10 is ELECT2
            # 12 is VDW2
            energies_datapoints = np.loadtxt(prod_alch, comments='#', usecols=[4, 6, 10, 12])

            # load metadata
            # NEW TI WINDOW: LAMBDA 0.3
            with open(prod_alch) as myfile:
                partition1, partition2 = [myfile.readline() for x in range(4)][-2:]
                # lines that we are parsing:
                # PARTITION 1 SCALING: BOND 1 VDW 0.3 ELEC 0
                # PARTITION 2 SCALING: BOND 1 VDW 0.7 ELEC 0.454545
                assert partition1.startswith('#PARTITION 1')
                assert partition2.startswith('#PARTITION 2')
                # fixme - use regex to get these
                app_vdw_x = float(partition1.split('VDW')[1].split('ELEC')[0])
                assert app_vdw_x == float(str(lambda_dir).split('_')[1])
                dis_vdw_x = float(partition2.split('VDW')[1].split('ELEC')[0])
                assert np.isclose(app_vdw_x, 1 - dis_vdw_x)
                app_ele_x = float(partition1.split('ELEC')[1])
                dis_ele_x = float(partition2.split('ELEC')[1])

                if app_vdw_x not in data['avdw']:
                    data['avdw'][app_vdw_x] = []
                if dis_vdw_x not in data['dvdw']:
                    data['dvdw'][dis_vdw_x] = []
                if app_ele_x not in data['aele']:
                    data['aele'][app_ele_x] = []
                if dis_ele_x not in data['dele']:
                    data['dele'][dis_ele_x] = []

                # add to the right dataset
                data['avdw'][app_vdw_x].append(energies_datapoints[:, 1])
                data['dvdw'][dis_vdw_x].append(energies_datapoints[:, 3])
                data['aele'][app_ele_x].append(energies_datapoints[:, 0])
                data['dele'][dis_ele_x].append(energies_datapoints[:, 2])

    return data


def analyse(data, location):
    # do some processing
    meta = {}
    for inter_type, dataset in data.items():
        if type(inter_type) is float:
            continue
        # calculate useful value for each data set
        # calc the mean value for all the replicas merged together
        stats = populate_stats(dataset)
        meta[inter_type] = stats

        continue

        # get running average
        # running_avg = pd.Series(merged_replicas).rolling(window=len(merged_replicas)).mean().iloc[len(merged_replicas) - 1:].values
        pd_frame = pd.DataFrame({'B': merged_replicas})
        window = int(0.1 * len(merged_replicas))
        rolling_mean = pd_frame.rolling(window).sum() / window
        cumu_sum = accumulate(merged_replicas)
        cum_rolling_avg = (accu / i for i, accu in enumerate(cumu_sum, 1))
        # plt.plot(list(cum_rolling_avg), label='my calc')
        # plt.plot(item['cumavg_avdw'][0], label='namd')
        # plt.plot(rolling_mean, label='rolling mean 10%')
        # plt.legend()
        # plt.show()

    # plot the average of the entire datasets now
    plt.figure()
    avdw_means = np.array(list(meta['avdw']['merged_mean'].items())).T
    plt.plot(avdw_means[0], avdw_means[1], label='Appearing VdW means', linestyle='--')
    dvdw_means = np.array(list(meta['dvdw']['merged_mean'].items())).T
    plt.plot(dvdw_means[0], dvdw_means[1], label='Disappearing VdW', linestyle='--')
    aele_means = np.array(list(meta['aele']['merged_mean'].items())).T
    plt.plot(aele_means[0], aele_means[1], label='Appearing q')
    dele_means = np.array(list(meta['dele']['merged_mean'].items())).T
    plt.plot(dele_means[0], dele_means[1], label='Disappearing q')
    plt.title(location)
    plt.legend()
    # plt.show()
    plt.savefig('/home/dresio/' + location + '.png')

    # integrate
    print(location)
    # sort by lambda before integrating
    avdw_means.sort()
    print('avdw means sorted', avdw_means)
    avdw_int = np.trapz(avdw_means[1], x=avdw_means[0])
    print('int avdw', avdw_int)

    dvdw_means.sort()
    print('avdw means sorted', dvdw_means)
    dvdw_int = np.trapz(dvdw_means[1], x=dvdw_means[0])
    print('int dvdw_means', dvdw_int)

    aele_means.sort()
    print('avdw means sorted', aele_means)
    aele_int = np.trapz(aele_means[1], x=aele_means[0])
    print('int aele_means', aele_int)

    dele_means.sort()
    print('avdw means sorted', dele_means)
    dele_int = np.trapz(dele_means[1], x=dele_means[0])
    print('int dele_means', dele_int)

    # we subtract partition 1 (appearing) from partition 2 (disappearing)
    # this means that there is some energy change when the atoms disappear,
    # and there is again some energy change when the atoms appear,

    # we know that if the ligand2 is more energetically favourable
    # than the sign should be -, because we are talking about the transition
    # from l1 to l2 (left to right)
    # therefore, we should do l1-l2,
    # so let's go through some examples,
    # if l1 is -10, but l2 is -12, then -10--12=2, but -12--10=-2, so it is just sign,

    # return aele_scratch - dele_scratch + avdw_scratch - dvdw_scratch
    # return (dele_int - aele_int) + (dvdw_int - avdw_int)
    return dele_int + aele_int + dvdw_int + avdw_int

complex_all = get_energies('complex')
complex_delta = analyse(complex_all, 'complex')
lig_all = get_energies('lig')
lig_delta = analyse(lig_all, 'lig')
print("Delta Delta: ", complex_delta - lig_delta)
