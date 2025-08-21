"""
Go to each directory and submit each replica separately
Walks through the typical lambda_X.XX/repX structures and submits the submit.sh script
"""

import os
import sys
import subprocess

for lambda_dir in os.listdir("."):
    if not lambda_dir.startswith("lambda_"):
        continue

    # go into the replicas
    for rep in os.listdir(lambda_dir):
        # for each replica, go in and schedule the simulation
        if not rep.startswith("rep"):
            continue

        rep_dir = os.path.join(lambda_dir, rep)
        try:
            # each script is submitted from the correct directory
            output = subprocess.check_output(
                [
                    "sbatch",
                    "--job-name=l%sr%s" % (lambda_dir.split("_")[1], rep[3:]),
                    "submit.sh",
                ],
                cwd=rep_dir,
            )
        except Exception as e:
            print(e)
            sys.exit(0)

        # e.g. Submitted batch job 7599862
        if "Submitted batch job" not in str(output):
            print("There was a problem with submitting a job")
            print("Directory: %s" % rep_dir)
            print("Full output: %s" % output)
        else:
            print(output)
