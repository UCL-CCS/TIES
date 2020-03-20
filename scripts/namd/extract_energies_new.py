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


def extract_energies(location):
    """
    @location - referes to the main locations where the lambda directories reside
    """

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

            # extract the raw datapoints
            # 4 is ELECT1
            # 6 is VDW1
            # 10 is ELECT2
            # 12 is VDW2
            energies_datapoints = np.loadtxt(prod_alch, comments='#', usecols=[4, 6, 10, 12])

            # load metadata from the file
            with open(prod_alch) as myfile:
                partition1, partition2 = [myfile.readline() for x in range(4)][-2:]
                # lines that we are parsing:
                # PARTITION 1 SCALING: BOND 1 VDW 0.3 ELEC 0
                # PARTITION 2 SCALING: BOND 1 VDW 0.7 ELEC 0.454545
                assert partition1.startswith('#PARTITION 1')
                assert partition2.startswith('#PARTITION 2')
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


def choder_get_eqpart(datapoints):
    """
    Extracts the equilibriated part
    """
    [t0, g, Neff_max] = timeseries.detectEquilibration(datapoints)
    print ('Chodera: t0 is', t0)
    return datapoints[t0:]


def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[int(result.size / 2):]


def get_replicas_stats(dataset, choderas_cut=False):
    meta = {
        'merged_mean': {},
        'merged_sem': {},
        'merged_std': {},
    }
    for lambda_val, replicas in dataset.items():
        # extract all the data points
        merged_replicas = []
        for rep in replicas:
            if choderas_cut:
                merged_replicas.extend(choder_get_eqpart(rep))
            else:
                merged_replicas.extend(rep)

        # record the observables
        meta['merged_mean'][lambda_val] = np.mean(merged_replicas)
        meta['merged_sem'][lambda_val] = sem(merged_replicas)
        meta['merged_std'][lambda_val] = np.std(merged_replicas)
        # fixme - add other errors, add bootstrapping here
    return meta


def analyse(data, location, choderas_cut=False):
    """
    Process the timeseries from each replica
    """

    print('Choderas cut turned on:', choderas_cut)

    # apply to each dataset
    stats = {}
    for interaction_type, dataset in data.items():
        stats[interaction_type] = get_replicas_stats(dataset, choderas_cut=choderas_cut)

    # plot the average of the entire datasets now
    plt.figure()
    avdw_means = np.array(list(stats['avdw']['merged_mean'].items())).T
    plt.plot(avdw_means[0], avdw_means[1], label='Appearing VdW means', linestyle='-', alpha=0.7)
    dvdw_means = np.array(list(stats['dvdw']['merged_mean'].items())).T
    plt.plot(dvdw_means[0][::-1], dvdw_means[1], label='Disappearing VdW', linestyle='--', alpha=0.7)
    aele_means = np.array(list(stats['aele']['merged_mean'].items())).T
    plt.plot(aele_means[0], aele_means[1], label='Appearing q', linestyle='-', alpha=0.7)
    dele_means = np.array(list(stats['dele']['merged_mean'].items())).T
    plt.plot(dele_means[0][::-1], dele_means[1], label='Disappearing q', linestyle='--', alpha=0.7)
    plt.title(location)
    plt.legend()
    # plt.show()
    plt.savefig(location + '.png')

    # integrate over the means from each replica
    print(location)
    # appearing vdw should have lambda values from 0 to 1,
    assert all([x < y for x, y in zip(avdw_means[0], avdw_means[0][1:])])
    avdw_int = np.trapz(avdw_means[1], x=avdw_means[0])
    print('int avdw', avdw_int)

    # disappearing vdw should have lambda values from 1 to 0,
    assert all([x > y for x, y in zip(dvdw_means[0], dvdw_means[0][1:])])
    dvdw_int = np.trapz(dvdw_means[1], x=dvdw_means[0])
    print('int dvdw_means', dvdw_int)

    # appearing ele should have lambda values from 0 to 1,
    assert all([x < y for x, y in zip(aele_means[0], aele_means[0][1:])])
    aele_int = np.trapz(aele_means[1], x=aele_means[0])
    print('int aele_means', aele_int)

    # disappearing ele should have lambda values from 1 to 0,
    assert all([x > y for x, y in zip(dele_means[0], dele_means[0][1:])])
    dele_int = np.trapz(dele_means[1], x=dele_means[0])
    print('int dele_means', dele_int)

    # return the final Delta G. Note that the sign in each delta G depends on the atoms contribution.
    return dele_int + aele_int + dvdw_int + avdw_int

complex_all = extract_energies('complex')
complex_delta = analyse(complex_all, 'complex')
lig_all = extract_energies('lig')
lig_delta = analyse(lig_all, 'lig')

print("Delta ligand", lig_delta)
print("Delta complex", complex_delta)
print("Delta Delta: ", complex_delta - lig_delta)
