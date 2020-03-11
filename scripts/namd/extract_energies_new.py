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


def get_energies(location):
    # extract from each replica the right
    vdw_lambdas = []
    dis_avgs_ele = []
    dis_ele_xs = []
    app_avgs_ele = []
    app_ele_xs = []
    averages_vdw = []

    # take only the lambda directories
    lambda_dirs = filter(lambda d: d.is_dir(), Path(location).glob(r'lambda_0.[0-9]*'))
    # sort lambda directories in the increasing order
    lambda_dirs = sorted(lambda_dirs, key=lambda d: float(d.name.split('_')[1]))

    # data contains lambdas
    data = {}

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
        data[next_lambda] = {'vdw': [], 'aele': [], 'dele': []}

        dis_replicas_avgs_ele = []
        app_replicas_avgs_ele = []
        replicas_avgs_vdw = []

        for rep in lambda_dir.glob('rep[0-9]*'):
            if not rep.is_dir():
                continue

            prod_alch = rep / 'prod.alch'

            # COLUMNS:
            # 5 is AVGELECT1
            # 7 is AVGVDW1
            # 11 is AVGELECT2
            # 13 is AVGVDW2
            # we take the last row because these are the running averages computed by NAMD
            energies = np.loadtxt(prod_alch, comments='#', usecols=[5, 7, 11, 13])[-1]

            # load metadata
            # NEW TI WINDOW: LAMBDA 0.3
            # PARTITION 1 SCALING: BOND 1 VDW 0.3 ELEC 0
            # PARTITION 2 SCALING: BOND 1 VDW 0.7 ELEC 0.454545
            with open(prod_alch) as myfile:
                partition1, partition2 = [myfile.readline() for x in range(4)][-2:]
                assert partition1.startswith('#PARTITION 1')
                assert partition2.startswith('#PARTITION 2')
                app_ele_x = float(partition1.split('ELEC')[1])
                dis_ele_x = float(partition2.split('ELEC')[1])

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
    plt.legend()
    plt.show()

    return int_dis_ele + int_app_ele + int_vdw

complex_all = get_energies('complex')
lig_all = get_energies('lig')
print("Final: ", complex_all - lig_all)