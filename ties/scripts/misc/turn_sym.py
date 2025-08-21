"""
Walk through the generated structure and turn the files into symbolic links if possible.
"""

import os


def check(dirname):
    for lambda_dir in os.listdir(dirname):
        if not lambda_dir.startswith("lambda_"):
            continue

        # go into the replicas
        for rep in os.listdir(os.path.join(dirname, lambda_dir)):
            # for each replica, go in and schedule the simulation
            if not rep.startswith("rep"):
                continue

            base_dir = os.getcwd()
            os.chdir(os.path.join(dirname, lambda_dir, rep))

            # files:
            files = [
                "sys_solv.pdb",
                "sys_solv_fep.pdb",
                "sys_solv.top",
                "constraint_f1.pdb",
                "constraint_f2.pdb",
                "constraint_f3.pdb",
                "constraint_f4.pdb",
                "morph_solv.top",
                "morph_solv.pdb",
                "morph_solv_fep.pdb",
            ]

            # rm the file and replace it with a symbolic link
            for filename in files:
                if os.path.isfile(filename):
                    os.remove(filename)
                    os.symlink("../../" + filename, filename)
                    print("replaced")

            # return to the root directory
            os.chdir(base_dir)


check("lig")
check("complex")
