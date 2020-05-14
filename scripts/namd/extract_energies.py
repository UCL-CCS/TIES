"""
# get from each lambda/replica the energies, sum the average vdw and avg ele
# todo - upgrade to pathlib

Note that the previous way to calculate error also takes place here.
It is described here. The gist:
- average all energies in a single replica to a single value
** if each energy contrib has 2500 datapoints, then that's 10k datapoints / 2500 time.
-
Paper: https://pubs.acs.org/doi/10.1021/acs.jctc.6b00979
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import sem
# import pandas as pd
from itertools import accumulate
# from pymbar import timeseries
from collections import OrderedDict
from scipy import interpolate


eq_steps = 500 # steps, 500 steps translates to 1 ns

def extract_energies(location):
    """
    @location - referes to the main locations where the lambda directories reside
    """

    # take only the lambda directories
    lambda_dirs = filter(lambda d: d.is_dir(), Path(location).glob(r'lambda_[0-1].[0-9]*'))
    # sort lambda directories in the increasing order
    lambda_dirs = sorted(lambda_dirs, key=lambda d: float(d.name.split('_')[1]))

    # different datasets, add bonded information for the future. It is not used now.
    data = {
        'dvdw': OrderedDict(), 'dele': OrderedDict(), 'dbon': OrderedDict(),
        'avdw': OrderedDict(), 'aele': OrderedDict(), 'abon': OrderedDict(),

        # this is for backward compatiblity with the previous error quantification
        'total_average': {}
    }

    for lambda_dir in lambda_dirs:
        dir_lambda_val = float(str(lambda_dir).split('_')[1])
        data['total_average'][dir_lambda_val] = []

        fresh_lambda = True
        ignore_dele_lambda = False
        for rep in lambda_dir.glob('rep[0-9]*'):
            if not rep.is_dir():
                continue

            # prod_alch = rep / 'prod_merged.alch'
            prod_alch = rep / 'prod.alch'
            if not prod_alch.is_file():
                print("A missing file: ", prod_alch)
                continue

            # partition1 is appearing
            # partition2 is disappearing

            # extract the raw datapoints
            # 2 is BOND1
            # 4 is ELECT1
            # 6 is VDW1
            # 8 is BOND2
            # 10 is ELECT2
            # 12 is VDW2
            energies_datapoints = np.loadtxt(prod_alch, comments='#', usecols=[2, 4, 6, 8, 10, 12])

            # load metadata from the file
            with open(prod_alch) as myfile:
                partition1, partition2 = [myfile.readline() for x in range(4)][-2:]
                # lines that we are parsing:
                # PARTITION 1 SCALING: BOND 1 VDW 0.3 ELEC 0
                # PARTITION 2 SCALING: BOND 1 VDW 0.7 ELEC 0.454545
                assert partition1.startswith('#PARTITION 1')
                assert partition2.startswith('#PARTITION 2')
                app_vdw_lambda = float(partition1.split('VDW')[1].split('ELEC')[0])
                assert app_vdw_lambda == float(str(lambda_dir).split('_')[1])
                dis_vdw_lambda = float(partition2.split('VDW')[1].split('ELEC')[0])
                assert np.isclose(app_vdw_lambda, 1 - dis_vdw_lambda)
                app_ele_lambda = float(partition1.split('ELEC')[1])
                dis_ele_lambda = float(partition2.split('ELEC')[1])

                if app_vdw_lambda not in data['avdw']:
                    data['avdw'][app_vdw_lambda] = []
                if dis_vdw_lambda not in data['dvdw']:
                    data['dvdw'][dis_vdw_lambda] = []
                if app_ele_lambda not in data['aele']:
                    data['aele'][app_ele_lambda] = []
                if dis_ele_lambda not in data['dele']:
                    data['dele'][dis_ele_lambda] = []

                if fresh_lambda:
                    # part 1 / aele, if this lambda is 0, then discard the previous derivative, the previous lambda 0
                    # was not the real the first lambda 0
                    if app_ele_lambda == 0:
                        data['aele'][app_ele_lambda] = []


                # to make it consitent with the previous calculataions, take data points from all
                this_replica_total = np.mean(energies_datapoints[eq_steps:, 2]) - np.mean(energies_datapoints[eq_steps:, 5]) + \
                                     np.mean(energies_datapoints[eq_steps:, 1]) - np.mean(energies_datapoints[eq_steps:, 4])
                data['total_average'][dir_lambda_val].append(this_replica_total)
                # ---------------------------------------------

                # add to the right dataset
                # add all values of VDW since they are appearing/disappearing throughout the lambda window
                data['avdw'][app_vdw_lambda].append(energies_datapoints[eq_steps:, 2])
                data['dvdw'][dis_vdw_lambda].append(energies_datapoints[eq_steps:, 5])

                # the appearing electrostatics being to appear around midway through the simulation,
                # we keep the very first lambda 0, not the previous states
                data['aele'][app_ele_lambda].append(energies_datapoints[eq_steps:, 1])

                # add part2/dele only the first time.
                # this means that ocne dele disappears, than any changes to dV/dl later are due to something else
                if fresh_lambda and len(data['dele'][dis_ele_lambda]) != 0:
                    ignore_dele_lambda = True
                if not ignore_dele_lambda:
                    data['dele'][dis_ele_lambda].append(energies_datapoints[eq_steps:, 4])

            fresh_lambda = False

    return data


def choder_get_eqpart(datapoints):
    """
    Extracts the equilibriated part
    """
    from pymbar import timeseries
    [t0, g, Neff_max] = timeseries.detectEquilibration(datapoints)
    # print('Chodera: t0 is', t0)
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


def get_int(xs, ys, interp=True):
    assert np.all(np.diff(xs) > 0)

    if interp:
        xnew = np.linspace(0, 1, num=50)
        # cubic_interp = interpolate.interp1d(dvdw_means[0], dvdw_means[1], kind='cubic')
        tck = interpolate.splrep(xs, ys, s=0)
        ynew = interpolate.splev(xnew, tck, der=False)

        plt.plot(xs, ys, label='original')
        # plt.plot(xnew, ynew, label='spline der0')
        # plt.legend()
        # plt.show()
        return np.trapz(ynew, x=xnew)

    # the xs (or lambdas) should be in a growing order
    return np.trapz(ys, x=xs)

def bootstrap_replica_averages(data):
    """
    Reproducing the way the error was calculated before.
    To do this, we have to extract for each replica the "derivative mean".
    This means that we clump all the derivatives together.
    Paper: https://pubs.acs.org/doi/10.1021/acs.jctc.6b00979
    Here we bootstrap the data, and create 10k means for each of the lambda window.
    """

    # fixme - so this is incorrect, but also the only way it is relevant to the past?
    bootstrapped_sem = OrderedDict()
    for lambda_val, tot_mean in data['total_average'].items():
        # for each lambda value, sample the means
        means = []
        # create 10 thousand of means
        for i in range(10 * 1000):
            means.append(np.mean(np.random.choice(tot_mean, size=len(tot_mean), replace=True)))
        bootstrapped_sem[lambda_val] = np.std(means)
        # it should be the STD, but have you found it in his code?
        # bootstrapped_sem[lambda_val] = np.std(means)

    # multiply by each
    k = list(bootstrapped_sem.keys())
    sigma_sum = 0
    for lambda_from, lambda_to, bsem in zip(k, k[1:], bootstrapped_sem.values()):
        assert lambda_from < lambda_to
        sigma_sum += (lambda_to - lambda_from)**2 * bsem**2

    data['sigma_sum'] = sigma_sum


def analyse(data, location, choderas_cut=False):
    """
    Process the timeseries from each replica
    """
    if choderas_cut:
        print('Choderas cut turned on')

    bootstrap_replica_averages(data)

    # apply to each dataset
    stats = {}
    for interaction_type, dataset in data.items():
        if interaction_type in ['avdw', 'dvdw', 'aele', 'dele']:
            stats[interaction_type] = get_replicas_stats(dataset, choderas_cut=choderas_cut)


    # plot the average of the entire datasets now
    # sort all lambdas from 0 to 1
    plt.figure()
    avdw_before_sort = list(stats['avdw']['merged_mean'].items())
    avdw_means = np.array(sorted(avdw_before_sort, key=lambda x: x[0])).T
    plt.plot(avdw_means[0], avdw_means[1], label='Appearing VdW means', linestyle='-', alpha=0.7)

    dvdw_before_sort = list(stats['dvdw']['merged_mean'].items())
    dvdw_means = np.array(sorted(dvdw_before_sort, key=lambda x: x[0])).T
    plt.plot(dvdw_means[0][::-1], dvdw_means[1], label='Disappearing VdW', linestyle='--', alpha=0.7)

    aele_before_sort = list(stats['aele']['merged_mean'].items())
    aele_means = np.array(sorted(aele_before_sort, key=lambda x: x[0])).T
    plt.plot(aele_means[0], aele_means[1], label='Appearing q', linestyle='-', alpha=0.7)

    dele_before_sort = list(stats['dele']['merged_mean'].items())
    dele_means = np.array(sorted(dele_before_sort, key=lambda x: x[0])).T
    plt.plot(dele_means[0][::-1], dele_means[1], label='Disappearing q', linestyle='--', alpha=0.7)
    plt.title(location)
    plt.legend()
    # plt.show()
    plt.savefig(location + '.png')
    plt.cla()

    # integrate over the means from each replica
    avdw_int = get_int(avdw_means[0], avdw_means[1])
    dvdw_int = get_int(dvdw_means[0], dvdw_means[1])
    aele_int = get_int(aele_means[0], aele_means[1])
    dele_int = get_int(dele_means[0], dele_means[1])

    out = f"""-------------------------{location:^10s}----------------------
                  Elec         vdW       Subtotal
---------------------------------------------------------
Part 1(app)     {aele_int:7.4f}  |  {avdw_int:7.4f}  | {aele_int + avdw_int:7.4f}     
Part 2(disapp)  {dele_int:7.4f}   |  {dvdw_int:7.4f}  | {dele_int + dvdw_int:7.4f}
---------------------------------------------------------
Subtotal        {aele_int - dele_int:7.4f}  |  {avdw_int - dvdw_int:7.4f}  | {aele_int + avdw_int - dvdw_int - dele_int:7.4f}
---------------------------------------------------------"""
    print(out)

    # return the final Delta G. Note that the sign in each delta G depends on the atoms contribution.
    return aele_int, avdw_int, dvdw_int, dele_int, data

choderas_cut = False
lig_all = extract_energies('lig')
laele_int, lavdw_int, ldvdw_int, ldele_int, lig_data = analyse(lig_all, 'lig', choderas_cut=choderas_cut)
lig_delta = laele_int + lavdw_int - ldvdw_int - ldele_int
complex_all = extract_energies('complex')
caele_int, cavdw_int, cdvdw_int, cdele_int, complex_data = analyse(complex_all, 'complex', choderas_cut=choderas_cut)
complex_delta = caele_int + cavdw_int - cdvdw_int - cdele_int

# print("Delta ligand", lig_delta)
# print("Delta complex", complex_delta)
print("Delta Delta: ", complex_delta - lig_delta)
print ("Agastya Error", complex_data['sigma_sum'] + lig_data['sigma_sum'])
