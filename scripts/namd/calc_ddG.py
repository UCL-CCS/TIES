"""
Calculate the ddG

In addition, this script contains error calculation.
In one case, it carries out the error calculation as in the
paper: https://pubs.acs.org/doi/10.1021/acs.jctc.6b00979
Furthermore, additional energy is calculated by sampling
the dv/dl and creating a distribution of ddG to estimate its error
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
import glob
import time


def merge_prod_files(files, output_merged_filename):
    # find the original prod.alch file first
    orig_prod = list(filter(lambda f:f.endswith('prod.alch'), files))[0]
    other_prods =  list(filter(lambda f:not f.endswith('prod.alch'), files))

    # take the prod.alch as teh reference, from every other remove the comments #
    lines = open(orig_prod).readlines()
    for other_prod in other_prods:

        next_lines = open(other_prod).readlines()
        # remove the comments
        data = filter(lambda l: not l.startswith('#'), next_lines)
        lines.extend(data)
    # save the results
    open(output_merged_filename, 'w').writelines(lines)


def extract_energies(location, choderas_cut=False, eq_steps=1000):
    """
    @location - referes to the main locations where the lambda directories reside.
    @eq_steps = 500 to be discarded as EQ stage.
    @choderas_cut - use Chodera's script to determine which part is the EQ part
    """
    print('Extracting energies from', location)
    if choderas_cut:
        print('Choderas cut turned on: it will be used to decide how much of the initial rep should be discarded')
    else:
        print('Using EQ cutoff to discard this many first steps:', eq_steps)

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

            # check if there are multiple .alch files, this means the restart was used and needs to be accounted for,
            # you could merge the results in that case and use them instead
            prod_files = glob.glob(str(rep) + os.sep + "prod*.alch")
            if len(prod_files) == 0:
                print("A missing energy file: prod.alch in %s" % rep)
                continue
            elif len(prod_files) == 1:
                prod_alch = rep / 'prod.alch'
                assert prod_alch.is_file()
            elif len(prod_files) > 1:
                # merge the files into a single file
                merged_prod = rep / 'prod_merged.alch'
                merge_prod_files(prod_files, merged_prod)
                # use the merged .alch instead
                prod_alch = merged_prod

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

                # fixme - this should be moved to the stage where at which the data is read
                # fixme
                if choderas_cut:
                    # use Chodera equilibration cutoff to decie which part of the production time
                    # should be discarded
                    avdw_dvdl = choder_get_eqpart(energies_datapoints[eq_steps:, 2])
                    dvdw_dvdl = choder_get_eqpart(energies_datapoints[eq_steps:, 5])
                    aele_dvdl = choder_get_eqpart(energies_datapoints[eq_steps:, 1])
                    dele_dvdl = choder_get_eqpart(energies_datapoints[eq_steps:, 4])
                else:
                    avdw_dvdl = energies_datapoints[eq_steps:, 2]
                    dvdw_dvdl = energies_datapoints[eq_steps:, 5]
                    aele_dvdl = energies_datapoints[eq_steps:, 1]
                    dele_dvdl = energies_datapoints[eq_steps:, 4]

                if not all([len(avdw_dvdl), len(dvdw_dvdl), len(aele_dvdl), len(dele_dvdl)]):
                    print('Not enough points. Ignoring ', lambda_dir, ' Replica', rep)
                    continue

                # --------------------------------------------------------------------------
                # to make it consitent with the previous calculataions, take data points from all
                # we always use the VDW energies
                this_replica_total = np.mean(avdw_dvdl) - np.mean(dvdw_dvdl)
                # we want to to take the dele values up until 0.6 (inclusive) where for the first time it's 0
                # and aele has its first 0 at 0.4 value

                if dir_lambda_val < 0.4:
                    # take only the dele
                    this_replica_total += np.mean(dele_dvdl)
                elif dir_lambda_val >= 0.4 and dir_lambda_val <= 0.6:
                    # take the difference of both
                    this_replica_total += np.mean(aele_dvdl) - np.mean(dele_dvdl)
                else:
                    # take only the aele value
                    this_replica_total += np.mean(aele_dvdl)

                data['total_average'][dir_lambda_val].append(this_replica_total)
                # ---------------------------------------------

                # add to the right dataset
                # add all values of VDW since they are appearing/disappearing throughout the lambda window
                data['avdw'][app_vdw_lambda].append(avdw_dvdl)
                data['dvdw'][dis_vdw_lambda].append(dvdw_dvdl)

                # the appearing electrostatics being to appear around midway through the simulation,
                # we keep the very first lambda 0, not the previous states
                data['aele'][app_ele_lambda].append(aele_dvdl)

                # add part2/dele only the first time.
                # this means that ocne dele disappears, than any changes to dV/dl later are due to something else
                if fresh_lambda and len(data['dele'][dis_ele_lambda]) != 0:
                    ignore_dele_lambda = True
                if not ignore_dele_lambda:
                    data['dele'][dis_ele_lambda].append(dele_dvdl)

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


def get_replicas_stats(dataset, sample_reps):
    meta = {
        'merged_mean': {},
        'merged_sem': {},
        'merged_std': {},
    }
    for lambda_val, replicas in dataset.items():
        # extract all the data points
        merged_replicas = []
        for rep in replicas:
            merged_replicas.extend(rep)

        if sample_reps:
            merged_replicas = np.random.choice(merged_replicas, size=len(merged_replicas), replace=True)

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
    In other words, if we have 5 replicas per lambda, we have only 5 values per lambda.
    Paper: https://pubs.acs.org/doi/10.1021/acs.jctc.6b00979
    Here we bootstrap the data, and create 10k means for each of the lambda window.
    In this case, we sample the actual means.
    Then at the end we take the SD.
    So we end up with the bootstrapped SD of means of means.
    """
    bootstrapped_sem = OrderedDict()
    for lambda_val, tot_mean in data['total_average'].items():
        # for each lambda value, sample the means
        means = []
        # create 10 thousand of means (of means)
        for i in range(10 * 1000):
            means.append(np.mean(np.random.choice(tot_mean, size=len(tot_mean), replace=True)))
        bootstrapped_sem[lambda_val] = np.std(means)

    # integrate the error with trapz
    sigma_int = np.trapz(y=list(bootstrapped_sem.values()), x=list(bootstrapped_sem.keys()))
    data['sigma_int'] = sigma_int


def analyse(data, location, calc_aga_err=False, sample_reps=False, verbose=True, plot=True):
    """
    Process the timeseries from each replica
    """

    if calc_aga_err:
        bootstrap_replica_averages(data)

    # apply to each dataset
    stats = {}
    for interaction_type, dataset in data.items():
        if interaction_type in ['avdw', 'dvdw', 'aele', 'dele']:
            stats[interaction_type] = get_replicas_stats(dataset, sample_reps=sample_reps)

    # plot the average of the entire datasets now
    # sort all lambdas from 0 to 1

    avdw_before_sort = list(stats['avdw']['merged_mean'].items())
    avdw_means = np.array(sorted(avdw_before_sort, key=lambda x: x[0])).T

    dvdw_before_sort = list(stats['dvdw']['merged_mean'].items())
    dvdw_means = np.array(sorted(dvdw_before_sort, key=lambda x: x[0])).T

    aele_before_sort = list(stats['aele']['merged_mean'].items())
    aele_means = np.array(sorted(aele_before_sort, key=lambda x: x[0])).T

    dele_before_sort = list(stats['dele']['merged_mean'].items())
    dele_means = np.array(sorted(dele_before_sort, key=lambda x: x[0])).T

    if plot:
        plt.figure()
        plt.plot(avdw_means[0], avdw_means[1], label='Appearing VdW means', linestyle='-', alpha=0.7)
        plt.plot(dvdw_means[0][::-1], dvdw_means[1], label='Disappearing VdW', linestyle='--', alpha=0.7)
        plt.plot(aele_means[0], aele_means[1], label='Appearing q', linestyle='-', alpha=0.7)
        plt.plot(dele_means[0][::-1], dele_means[1], label='Disappearing q', linestyle='--', alpha=0.7)
        plt.title(location)
        plt.legend()
        # plt.show()
        plt.savefig(location + '_dvdl.png')
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
    if verbose:
        print(out)

    # return the final Delta G. Note that the sign in each delta G depends on the atoms contribution.
    return aele_int, avdw_int, dvdw_int, dele_int, data


def bootstrapped_ddG(ligand_data, complex_data, bootstrap_size_k=1):
    """
    Use bootstrapped data to estimate the DDG. This requries access to both, the ligand and the complex data.
    Each time we resample our dataset. This means that for ligand/complex, for each lambda window,
    we sample the different dV/dL. Then we recreate the means.
    """
    # do bootstrapping to find out error for each component
    # we want to know where the difference comes from
    # fixme - try bootstrapping to understand the error in aele, dele etc
    bootstrapped_ddGs = []
    for i in range(bootstrap_size_k * 1000):  # fixme
        laele_int, lavdw_int, ldvdw_int, ldele_int, lig_data = analyse(lig_all, 'lig', calc_aga_err=False,
                                                                       verbose=False, sample_reps=True,
                                                                       plot=False)
        lig_delta = laele_int + lavdw_int - ldvdw_int - ldele_int
        caele_int, cavdw_int, cdvdw_int, cdele_int, complex_data = analyse(complex_all, 'complex',
                                                                           calc_aga_err=False,
                                                                           verbose=False, sample_reps=True,
                                                                           plot=False)
        complex_delta = caele_int + cavdw_int - cdvdw_int - cdele_int

        bootstrapped_ddGs.append(complex_delta - lig_delta)

    # return the standard error
    return bootstrapped_ddGs

t_start = time.time()
choderas_cut = False
calc_aga_err = True
lig_all = extract_energies('lig', choderas_cut=choderas_cut)
laele_int, lavdw_int, ldvdw_int, ldele_int, lig_data = analyse(lig_all, 'lig', calc_aga_err=calc_aga_err, verbose=True)
lig_delta = laele_int + lavdw_int - ldvdw_int - ldele_int
complex_all = extract_energies('complex', choderas_cut=choderas_cut)
caele_int, cavdw_int, cdvdw_int, cdele_int, complex_data = analyse(complex_all, 'complex', calc_aga_err=calc_aga_err, verbose=True)
complex_delta = caele_int + cavdw_int - cdvdw_int - cdele_int

# Give the overall results
print("Delta Delta: ", complex_delta - lig_delta)
print ("Agastya Error", complex_data['sigma_int'] + lig_data['sigma_int'])

# now that we have the bootstrapped_ddGs, we take SD to find the standard error in the bootstrapped ddG
# se_bootstrapped_ddG = bootstrapped_ddG(lig_data, complex_data, 1)
# print('The bootstrapped mean of ddG is', np.mean(se_bootstrapped_ddG))
# print('The bootstrapped standard error of ddG is', np.std(se_bootstrapped_ddG))
# print(os.linesep + os.linesep + 'Altogether analysis time(s)', time.time() - t_start)

"""
Bootstrapping performance upgrade:
use the new numpy and select as my times as you want:

x = array([[1.        , 1.11111111, 1.22222222, 1.33333333, 1.44444444,
        1.55555556, 1.66666667, 1.77777778, 1.88888889, 2.        ],
       [2.        , 2.11111111, 2.22222222, 2.33333333, 2.44444444,
        2.55555556, 2.66666667, 2.77777778, 2.88888889, 3.        ]])
rng.choice(x, size=10, replace=True, axis=1)
rng = np.random.default_rng()

in other words, we can apply straight away 5 k and generate all the data in one go and avoid calling so many functions
then, ideally we would call the integral on all of them as well
"""