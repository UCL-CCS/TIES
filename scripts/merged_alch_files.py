"""
Walks through the typical lambda_X.XX/repX structures and copies the file there
"""
import os
import sys
import shutil
import glob

for lambda_dir in os.listdir('.'):
    if not lambda_dir.startswith('lambda_'):
        continue

    # go into the replicas
    for rep in os.listdir(lambda_dir):
        # for each replica, go in and schedule the simulation
        if not rep.startswith('rep'):
            continue

        # check if there are multiple files in a single replica with the format prod*.alch
        files = glob.glob(os.path.join(lambda_dir, rep, 'prod*alch'))
        if len(files) == 0:
            print('Not a single prod file in this dir', lambda_dir, rep)
        elif len(files) == 1:
            shutil.copy(files[0], os.path.join(lambda_dir, rep, 'prod_merged.alch'))
        else:
            # take the prod.alch as teh reference, from every other remove the comments #
            lines = open(os.path.join(lambda_dir, rep, 'prod.alch')).readlines()
            for other_prod in files:
                if other_prod == 'prod.alch':
                    continue

                next_lines = open(other_prod).readlines()
                # remove the comments
                data = filter(lambda l:not l.startswith('#'), next_lines)
                lines.extend(data)
            # save the results
            open(os.path.join(lambda_dir, rep, 'prod_merged.alch'), 'w').writelines(lines)