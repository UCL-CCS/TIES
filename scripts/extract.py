"""
Read the namd output and sum all the energies, t
"""
import os
import sys
import subprocess
import numpy as np

for lambda_dir in os.listdir('.'):
    if not lambda_dir.startswith('lambda_0.00'):
        continue

    # go into the replicas
    for rep in os.listdir(lambda_dir):
        # for each replica, go in and schedule the simulation
        if not rep.startswith('rep'):
            continue

        rep_dir = os.path.join(lambda_dir, rep)
        # lambda_0.00/rep3/eq.log:[Partition 0][Node 0] End of program
        # lambda_0.00/rep3/prod.log:[Partition 0][Node 0] End of program
        prod_alch = os.path.join(rep_dir, 'prod.alch')
        if not os.path.isfile(prod_alch):
            print('Not found .alch', prod_alch)
            continue

        alch = np.loadtxt('prod.alch', comments=['#'], usecols=list(range(1, 13 + 1)))



