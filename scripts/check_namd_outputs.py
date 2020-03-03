"""
# would through the lambda directories
# and schedule each simulation
# fixme - could check if the simulation is finished/running etc, in that case ignore
"""
import os
import sys
import subprocess

for lambda_dir in os.listdir('.'):
    if not lambda_dir.startswith('lambda_'):
        continue

    # go into the replicas
    for rep in os.listdir(lambda_dir):
        # for each replica, go in and schedule the simulation
        if not rep.startswith('rep'):
            continue

        rep_dir = os.path.join(lambda_dir, rep)
        # lambda_0.00/rep3/eq.log:[Partition 0][Node 0] End of program
        # lambda_0.00/rep3/prod.log:[Partition 0][Node 0] End of program
        if not 'End of program' in open(os.path.join(rep_dir, 'eq_step4.log')).read():
            print('Eq not finished: %s' % rep_dir)
        if not 'End of program' in open(os.path.join(rep_dir, 'prod.log')).read():
            print('Prod not finished: %s' % rep_dir)
