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

def get_energies(location):
    # extract from each replica the right
    vdw_lambdas = []
    dis_avgs_ele = []
    dis_ele_xs = []
    app_avgs_ele = []
    app_ele_xs = []
    averages_vdw = []

    # take only the lambda directories
    lambda_dirs = filter(lambda d: d.is_dir(), Path(location).glob(r'lambda_[0-1].[0-9]*'))
    # sort lambda directories in the increasing order
    lambda_dirs = sorted(lambda_dirs, key=lambda d: float(d.name.split('_')[1]))

    # data uses lambdas as keys
    data = {'avdw': OrderedDict(), 'dvdw': OrderedDict(), 'aele': OrderedDict(), 'dele': OrderedDict()}

    for lambda_dir in lambda_dirs:
        dis_ele = None
        app_ele = None

        # ? we have to average the results over the replicas

        # ensure that lambdas are in the right order
        next_lambda = float(lambda_dir.name.split('_')[1])
        assert 1 >= next_lambda >= 0
        if vdw_lambdas:
            assert next_lambda > vdw_lambdas[-1]
        vdw_lambdas.append(next_lambda)

        # each lambda is going to have {vdw, appearing ele, disappearing ele}
        # for which there are going to be datapoints
        # data[next_lambda] = {'avdw': [], 'dvdw': [], 'aele': [], 'dele': [],
        #                      'cumavg_avdw': []}
        data[next_lambda] = {'cumavg_avdw': []}

        dis_replicas_avgs_ele = []
        app_replicas_avgs_ele = []
        replicas_avgs_vdw = []

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
            data[next_lambda]['cumavg_avdw'].append(energies_cum_roll[1])

            # extract the raw datapoints
            # 4 is ELECT1
            # 6 is VDW1
            # 10 is ELECT2
            # 12 is VDW2
            energies_datapoints = np.loadtxt(prod_alch, comments='#', usecols=[4, 6, 10, 12])

            # add info about the

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
                assert app_vdw_x == next_lambda
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

            if app_ele_x == dis_ele_x == 0:
                # only take the VDW terms
                vdw_deriv = -energies[1] + energies[3]
                app_ele = dis_ele = None
            elif app_ele_x == 0:
                vdw_deriv = -energies[1] + energies[3]
                dis_ele = energies[2]
                app_ele = None
            elif dis_ele_x == 0:
                vdw_deriv = -energies[1] + energies[3]
                app_ele = -energies[0]
                dis_ele = None
            else:
                vdw_deriv = -energies[1] + energies[3]
                dis_ele = energies[2]
                app_ele = -energies[0]

            # change the electrostatic terms dependant on the lambda
            # add them row by row to get all the energies
            # namd_averaged_deriv = -energies[0] - energies[1] + energies[2] + energies[3]
            # get the average of the second half
            # half2_avg = np.average(energies_over_time[int(len(energies_over_time)/2):])
            if dis_ele is not None:
                dis_replicas_avgs_ele.append(dis_ele)

            if app_ele is not None:
                app_replicas_avgs_ele.append(app_ele)

            replicas_avgs_vdw.append(vdw_deriv)
            # plt.plot(energies_over_time)
            # plt.show()
            # break

        # use the average from the replicas
        averages_vdw.append(np.average(replicas_avgs_vdw))
        if dis_ele is not None:
            dis_avgs_ele.append(np.average(dis_replicas_avgs_ele))
            dis_ele_xs.append(dis_ele_x)

        if app_ele is not None:
            app_avgs_ele.append(np.average(app_replicas_avgs_ele))
            app_ele_xs.append(app_ele_x)

        # plt.plot(replicas_avgs)
        # plt.show()
        # break

    print('dis ele xs', dis_ele_xs)
    print('app ele xs', app_ele_xs)

    def populate_stats(dataset, itype):
        # chodera equilibration
        # reps = item[itype]
        # print(itype, 'has reps number ', len(reps))
        # for rep in reps:
        #     [t0, g, Neff_max] = timeseries.detectEquilibration(rep)
        #     # A_t_equlibrated = A_t[t0:]

        meta = {}
        for lambda_val, reps in dataset.items():
            if f'merged_mean' not in meta:
                meta[f'merged_mean'] = {}

            # interaction type
            merged_replicas = [datapoint for rep in reps for datapoint in rep]
            # item[f'merged_{itype}_sem'] = sem(merged_replicas)
            meta[f'merged_mean'][lambda_val] = np.mean(merged_replicas)
            # item[f'merged_{itype}_std'] = np.std(merged_replicas)
            # print('sem is overall', item['merged_avdw_sem'])
        return meta

    # do some processing
    meta = {}
    for inter_type, dataset in data.items():
        if type(inter_type) is float:
            continue
        # calculate useful value for each data set
        # calc the mean value for all the replicas merged together
        stats = populate_stats(dataset, inter_type)
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

        # plot altogether
        # plt.hist(merged_replicas)
        # plt.xlabel('Merged avdw')
        # plt.title(lambda_key)
        # plt.show()

    def autocorr(x):
        result = np.correlate(x, x, mode='full')
        return result[int(result.size / 2):]

    # plot the avdw vs lambda
    # test_xs = []
    # for lambda_key, item in data.items():
    #     reps = item['avdw']
    #     test_xs.append(lambda_key)

        # plot avdw with the different replicas
        # for rep in reps:
        #     plt.plot(autocorr(rep))
        # plt.xlabel('A vdw energies')
        # plt.show()

    # for each lambda, extract average avdw, dvdw, aele, dele and use that as the beginning in the integration
    # for inter_type, dataset in data.items():
    #     avdw_replicas = dataset['avdw']
    #     # combine all
    #     merged = [datapoint for rep in avdw_replicas for datapoint in rep]
    #     # get the average
    #     dataset['avdw_mean'] = np.mean(merged)
    #
    #     pass

    # plot the average of the entire datasets now

    plt.figure()
    avdw_means = np.array(list(meta['avdw']['merged_mean'].items())).T
    plt.plot(avdw_means[0], avdw_means[1], label='avdw_mean')
    dvdw_means = np.array(list(meta['dvdw']['merged_mean'].items())).T
    plt.plot(dvdw_means[0], dvdw_means[1], label='dvdw_mean')
    aele_means = np.array(list(meta['aele']['merged_mean'].items())).T
    plt.plot(aele_means[0], aele_means[1], label='aele_mean')
    dele_means = np.array(list(meta['dele']['merged_mean'].items())).T
    plt.plot(dele_means[0], dele_means[1], label='dele_mean')
    plt.title(location)
    plt.legend()
    # plt.show()
    plt.savefig('/home/dresio/' + location + '.png')
    # plt.cla()

    # integrate
    avdw_scratch = np.trapz(avdw_means[1], x=avdw_means[0])
    print('scratch avdw', avdw_scratch)
    dvdw_scratch = np.trapz(dvdw_means[1], x=dvdw_means[0])
    print('scratch dvdw_means', dvdw_scratch)
    aele_scratch = np.trapz(aele_means[1], x=aele_means[0])
    print('scratch aele_means', aele_scratch)
    dele_scratch = np.trapz(dele_means[1], x=dele_means[0])
    print('scratch dele_means', dele_scratch)

    # integrate
    int_vdw = np.trapz(averages_vdw, x=vdw_lambdas)
    # trapz_res = np.trapz(averages_vdw, x=[0.45])averages_vdw
    print("vdw integral", int_vdw)
    plt.plot(vdw_lambdas, averages_vdw, label='vdw')
    plt.plot(dis_ele_xs, dis_avgs_ele, label='dis ele')
    plt.plot(app_ele_xs, app_avgs_ele, label='app ele')

    int_dis_ele = np.trapz(dis_avgs_ele, x=dis_ele_xs)
    int_app_ele = np.trapz(app_avgs_ele, x=app_ele_xs)
    print("ele dis integral", int_dis_ele)
    print("ele app integral", int_app_ele)
    # plt.plot(lambdas, averages_ele, label='ele')
    print('Altogether:', int_dis_ele + int_app_ele + int_vdw)
    # plt.legend()
    # plt.show()

    # return int_dis_ele + int_app_ele + int_vdw
    return avdw_scratch + dvdw_scratch + aele_scratch + dele_scratch

complex_all = get_energies('complex')
lig_all = get_energies('lig')
print("Final: ", complex_all - lig_all)