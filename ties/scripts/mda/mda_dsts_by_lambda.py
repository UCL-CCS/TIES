"""
A quick script to scan the 65 replicas for some property.
The distance between two atoms is computed for all replicas and
plotted as a distribution for two different simulations.
"""

import os
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

def get_sims(dirname):
    finished_sims = {}
    for lambda_dir in os.listdir(dirname):
        if not lambda_dir.startswith('lambda_'):
            continue

        lambda_val = float(lambda_dir.split('_')[1])
        finished_sims[lambda_val] = []

        # go into the replicas
        for rep in os.listdir(os.path.join(dirname, lambda_dir)):
            # for each replica, go in and schedule the simulation
            if not rep.startswith('rep'):
                continue

            rep_dir = os.path.join(dirname, lambda_dir, rep)
            # find the simulation data (including the .BAK files present after reruns)
            prod_files = glob(os.path.join(rep_dir, 'prod.dcd*'))
            # check if there are files with extensions different than .dcd or .BAK
            for filename in prod_files:
                if not filename.lower().endswith('dcd') and not filename.lower().endswith('bak'):
                    raise Exception(f'A new extension type!: {filename}')

            # rename the .BAK files so MDAnalysis will understand the extension
            for filename in filter(lambda f:f.lower().endswith('bak'), prod_files):
                # rename the .BAK file to prod_int.dcd
                # generate a new name
                for i in range(10):
                    new_filename = os.path.join(rep_dir, f'prod_{i:d}.dcd')
                    if not os.path.exists(new_filename):
                        break
                #
                os.rename(filename, new_filename)

            renamed_files = glob(os.path.join(rep_dir, 'prod*.dcd'))
            finished_sims[lambda_val].extend(renamed_files)

    return finished_sims

def get_dsts(top_path, sims):
    print('Top: ', top_path)
    dsts_for_lambda = []
    for lambda_val, replicas in sorted(sims.items()):
        # if lambda_val != 0.4:
        #     continue
        # load
        u = mda.Universe(os.path.join(top_path, "sys_solv.top") , *replicas)
        print(f'Lambda {lambda_val}, Frames {u.trajectory.n_frames}')

        # C7 and C59  / CL8 BR1
        cl8 = u.select_atoms('name C9')
        br1 = u.select_atoms('name Na+')
        # check their distance over time
        dsts = []
        for ts in u.trajectory:
            dst = distance_array(cl8.positions, br1.positions, box=ts.dimensions)[0]
            dsts.extend(dst)

        dsts_for_lambda.append(dsts)

    # flat = []
    # [flat.extend(one_d) for one_d in dsts_for_lambda]

    return dsts_for_lambda


good_lig = '/home/dresio/ucl/validation/resp_validation/ptp1b/l6_l14_atol0.131_hsp/lig'
good_sims = get_sims(good_lig)
good_dsts = get_dsts(good_lig, good_sims)

bad_lig = '/home/dresio/ucl/validation/resp_validation/ptp1b/l6_l14_atol0.101_hsp/lig'
bad_sims = get_sims(bad_lig)
bad_dsts = get_dsts(bad_lig, bad_sims)

lambdas = sorted(good_sims.keys())

# nocross_lig = '/home/dresio/ucl/validation/resp_validation/ptp1b/l6_l14_atol0.101_no_cross_angdi_hsp/lig'
# nocross_sims = get_sims(nocross_lig)
# nocross_dsts = get_dsts(nocross_lig, nocross_sims)

counter = 1
plt.figure(figsize=(15,10))
plt.rcParams.update({'font.size': 6})

for lam, gdst, bdst in zip(lambdas, good_dsts, bad_dsts):
    plt.subplot(3, 5, counter)

    bin_from = min(min(bdst), min(gdst))
    bin_to = max(max(bdst), max(gdst))
    bin_to = 7
    bins = np.linspace(bin_from, bin_to, 20)
    # bins=20
    plt.hist(gdst, label='0.131', bins=bins, alpha=0.8)
    plt.hist(bdst, label='0.101', bins=bins, alpha=0.8)
    # plt.hist(nocross_flat, label='no cross 0.101', bins=bins, alpha=0.8)
    plt.title(f'C9 to Na+ Lambda {lam}')
    plt.ylabel('Freq')
    plt.xlabel('A')

    counter += 1

plt.legend()
plt.show()

# flatten altogether and compare
