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

        submit_sh_loc = os.path.join(lambda_dir, rep, 'submit.sh')
        try:
            output = subprocess.check_output(['sbatch', submit_sh_loc])
        except Exception as e:
            print(e)
            sys.exit(0)

        # fixme - add a name to each job: lambda + pro + protein name

        # e.g. Submitted batch job 7599862
        if "Submitted batch job" not in str(output):
            print('There was a problem with submitting a job')
            print('Directory: ', submit_sh_loc)
            print('Full output', output)
