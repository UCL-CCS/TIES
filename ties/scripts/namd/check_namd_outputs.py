"""
Check if the NAMD simulations finished for ligand and complex
"""

import os


def check(dirname):
    finished_sims = []
    for lambda_dir in os.listdir(dirname):
        if not lambda_dir.startswith("lambda_"):
            continue

        # go into the replicas
        for rep in os.listdir(os.path.join(dirname, lambda_dir)):
            # for each replica, go in and schedule the simulation
            if not rep.startswith("rep"):
                continue

            rep_dir = os.path.join(dirname, lambda_dir, rep)
            try:
                content = open(os.path.join(rep_dir, "prod.log")).read()
                if not (
                    "End of program" in content
                    or "WRITING VELOCITIES TO OUTPUT FILE AT STEP 3000000" in content
                ):
                    print("Prod not finished: %s" % rep_dir)
                else:
                    finished_sims.append(os.path.join(rep_dir, "prod.log"))
            except Exception:
                print("prod.log files not found: %s" % rep_dir)

    print("Total number of finished simulations: ", len(finished_sims))
    print("Which are: ", finished_sims)


check("lig")
check("complex")
