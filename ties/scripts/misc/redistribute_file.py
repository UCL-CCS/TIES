"""
Walks through the typical lambda_X.XX/repX structures and copies the file there
"""

import os
import sys
import shutil

redistribute_filename = sys.argv[1]
assert os.path.isfile(redistribute_filename)
print("Will redistribute the file: " + redistribute_filename)

for lambda_dir in os.listdir("."):
    if not lambda_dir.startswith("lambda_"):
        continue

    # go into the replicas
    for rep in os.listdir(lambda_dir):
        # for each replica, go in and schedule the simulation
        if not rep.startswith("rep"):
            continue

        dest_dir = os.path.join(lambda_dir, rep)
        print("Copying " + redistribute_filename + " to " + dest_dir)
        shutil.copy(redistribute_filename, dest_dir)
