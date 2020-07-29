#!/usr/bin/env python3
"""
Calculate the ddG

In addition, this script contains error calculation.
In one case, it carries out the error calculation as in the
paper: https://pubs.acs.org/doi/10.1021/acs.jctc.6b00979
Furthermore, additional energy is calculated by sampling
the dv/dl and creating a distribution of ddG to estimate its error
"""
import os
from pathlib import Path
from collections import OrderedDict
import glob
import time
import pickle as pkl

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from scipy.stats import sem
from scipy import interpolate
# import pandas as pd
from itertools import accumulate
# from pymbar import timeseries


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


def extract_energies(work_path, choderas_cut=False, eq_steps=1000):
    """
    @location - referes to the main locations where the lambda directories reside.
    @eq_steps = 500 to be discarded as EQ stage.
    @choderas_cut - use Chodera's script to determine which part is the EQ part
    @same_lambda_added - temporary
    @data - use multiple times in order to get data from more locations
    """
    print('Extracting energies from', work_path)
    if choderas_cut:
        print('Choderas cut turned on: it will be used to decide how much of the initial rep should be discarded')
    else:
        print('Using EQ cutoff to discard this many first steps:', eq_steps)

    # take only the lambda directories
    lambda_dirs = filter(lambda d: d.is_dir(), work_path.glob(r'lambda_[0-1].[0-9]*'))
    # sort lambda directories in the increasing order
    lambda_dirs = sorted(lambda_dirs, key=lambda d: float(d.name.split('_')[1]))

    # different datasets, add bonded information for the future. It is not used now.
    data = {
        'dvdw': OrderedDict(), 'dele': OrderedDict(),
        'avdw': OrderedDict(), 'aele': OrderedDict(),

        # this is for backward compatiblity with the previous error quantification
        'total_average': {},
        # this is similar but instead of recording the entire means, it records all the data points
        'added_series': {}
    }

    for lambda_dir in lambda_dirs:
        dir_lambda_val = float(str(lambda_dir).rsplit('_', maxsplit=1)[1])
        if dir_lambda_val not in data['total_average']:
            data['total_average'][dir_lambda_val] = []
            data['added_series'][dir_lambda_val] = []

        fresh_lambda = True
        ignore_dele_lambda = False
        for rep in sorted(lambda_dir.glob('rep[0-9]*')):
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
            print(prod_alch)
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
                dis_vdw_lambda = float(partition2.split('VDW')[1].split('ELEC')[0])

                # this is the default check that does not apply if you change the way lambda is created
                # assert app_vdw_lambda == float(str(lambda_dir).split('_')[1]), lambda_dir
                # assert np.isclose(app_vdw_lambda, 1 - dis_vdw_lambda)
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
                # extract the values
                data['added_series'][dir_lambda_val].append(avdw_dvdl - dvdw_dvdl + (aele_dvdl - dele_dvdl)/0.55)
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

        # replicas were iterated, check if they each have at least 3 replicas out of 5 done
        # fixme / todo

    return data


def saveSeriesTotalPerLambda(data):
    # save the series for each lambda together
    for lambda_window, reps in data['added_series'].items():
        # merge the replicas into one numpy array
        pass


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
    """
    To do: obtains stats for each replica separately, as well as for all of them together merged,
    todo
    """
    meta = {
        # combined
        'merged_mean': {},
        'merged_sem': {},
        'merged_std': {},
        # for each replica
        # 'rep_means': [],
    }
    for lambda_val, replicas in dataset.items():
        # extract all the data points
        merged_replicas = []
        for rep in replicas:
            merged_replicas.extend(rep)

        # this should not be here, create a separate function that carries this out
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

def bootstrap_replica_averages_err2017(data):
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

    # multiply by each
    k = list(bootstrapped_sem.keys())
    sigma_sum = 0
    for lambda_from, lambda_to, boot_point in zip(k, k[1:], bootstrapped_sem.values()):
        assert lambda_from < lambda_to
        sigma_sum += (lambda_to - lambda_from) ** 2 * boot_point ** 2

    # note that for some reason, at the end sqrt was taken of the value
    data['sigma_2017'] = np.sqrt(sigma_sum)


def bootstrap_replica_averages_improved(data):
    """
    This is a variation that is more likely to show the actual values. We simply integrate over variances at each point.
    Note that the variance at each point is not an optimal strategy: we will have very different variance
    coming from Coulomb and VDW. So what happens when we naively add their points together?
    In this case it might not matter because we take a single point which is a "distribution".
    Paper: https://pubs.acs.org/doi/10.1021/acs.jctc.6b00979
    """
    bootstrapped_sem = OrderedDict()
    for lambda_val, tot_mean in data['total_average'].items():
        # for each lambda value, sample the means
        means = []
        # create 10 thousand of means (of means)
        for i in range(10 * 1000):
            means.append(np.mean(np.random.choice(tot_mean, size=len(tot_mean), replace=True)))
        bootstrapped_sem[lambda_val] = np.var(means)

    # integrate the error with trapz
    sigma_int = np.trapz(y=list(bootstrapped_sem.values()), x=list(bootstrapped_sem.keys()))
    data['sigma_int'] = sigma_int


def analyse(data, system_conf, calc_aga_err=False, sample_reps=False, verbose=True, plot=True):
    """
    Process the timeseries from each replica
    """

    if calc_aga_err:
        bootstrap_replica_averages_err2017(data)

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
        # This is the main summary graph that merges all the replicas
        # Plot the lines that are integrated, together with their overall standard deviation
        plt.figure(figsize=(10, 8) ) #, dpi=80) facecolor='w', edgecolor='k')
        plt.axhline(color='grey', linestyle='dotted')
        # todo draw standard deviation instead
        plt.plot(avdw_means[0], avdw_means[1], label='Appearing VdW means', linestyle='-', alpha=0.7, color='blue')
        plt.plot(dvdw_means[0][::-1], dvdw_means[1], label='Disappearing VdW', linestyle='--', alpha=0.7, color='blue')
        plt.plot(aele_means[0], aele_means[1], label='Appearing q', linestyle='-', alpha=0.7, color='red')
        plt.plot(dele_means[0][::-1], dele_means[1], label='Disappearing q', linestyle='--', alpha=0.7, color='red')
        plt.title(system_conf)
        plt.legend()
        # plt.show()
        plt.savefig(os.path.join(analysis_dir, system_conf + '_dvdl.png'))
        plt.cla()

        # In this part we create more detailed plots for energy contribution
        # ie appearing and disappearing vdw are separated
        # Then, we draw boxplots to show where the mean is coming from
        def plot_single_component(means, data, label, file_suffix):
            """
            Plots a single component (for example appearing VDW).
            Shows the error bars for each replica.
            """
            plt.figure(figsize=(8, 10))  # , dpi=80) facecolor='w', edgecolor='k')
            plt.rcParams.update({'font.size': 15})

            if system_conf == 'lig':
                plt.title('Ligand')
            elif system_conf == 'complex':
                plt.title('Complex')

            for lam, reps in data.items():
                # plot the dv/dl boxplots
                # stagger lambdas so that replicas do not overlap
                staggered_lambdas = np.linspace(lam-0.025, lam+0.025, len(reps))
                for new_lambda, rep in zip(staggered_lambdas, reps):
                    plt.boxplot(rep, positions=[new_lambda, ],
                                widths=staggered_lambdas[1] - staggered_lambdas[0],
                                sym='.', showfliers=False) # showfliers=False

            # plot the global mean values
            plt.plot(means[0], means[1], label=label, linestyle='-', alpha=1)

            plt.ylabel('$\\rm \\left\\langle \\frac{dU}{d\\lambda} \\right\\rangle $')
            plt.xlabel('$\\rm \\lambda$')
            plt.xticks(means[0], [f"{m:.2f}" for m in means[0]], rotation=45)
            plt.xlim([-0.05, 1.05])
            plt.legend()
            # plt.show()
            plt.tight_layout()
            plt.savefig(os.path.join(analysis_dir, system_conf + file_suffix))
            plt.cla()

        plot_single_component(avdw_means, data['avdw'], 'Appearing VdW means', '_avdw_dvdl.png')
        plot_single_component(dvdw_means, data['dvdw'], 'Disappearing VdW means', '_dvdw_dvdl.png')
        plot_single_component(aele_means, data['aele'], 'Appearing q means', '_aele_dvdl.png')
        plot_single_component(dele_means, data['dele'], 'Disappearing q', '_dele_dvdl.png')


    # integrate over the means from each replica
    avdw_int = get_int(avdw_means[0], avdw_means[1])
    dvdw_int = get_int(dvdw_means[0], dvdw_means[1])
    aele_int = get_int(aele_means[0], aele_means[1])
    dele_int = get_int(dele_means[0], dele_means[1])

    out = f"""-------------------------{system_conf:^10s}----------------------
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


def merge_datasets(data_list):
    # Merge all datasets with the first one

    merged = data_list[0]
    for data in data_list[1:]:
        for key in ['dvdw', 'dele', 'avdw', 'aele', 'total_average']:
            for lambda_val, reps in data[key].items():
                merged[key][lambda_val].extend(reps)

    return merged
