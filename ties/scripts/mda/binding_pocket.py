import MDAnalysis as mda
import glob
import os
import numpy as np
from MDAnalysis.analysis.distances import distance_array
from collections import OrderedDict
import matplotlib.pyplot as plt

"""
Track the distance to the pocket from each atom. 
Note that in the case of mcl1 l2 l32 this did not mean match.
"""


# load the protein.pdb and the morph.mol2
# note that the morph is docked
sys = mda.Universe("complex/morph_solv.top", "complex/morph_solv.pdb")
print("sys", sys)

# establish the binding site (within 5A of the original position)
ligand = sys.select_atoms("resname ged")
binding_site = sys.select_atoms(
    "not resname WAT and not resname ged and around 5 (resname ged)"
)
print(f"Binding site has atoms: {len(binding_site)}")


def get_avg_min_dst(ligpos, sitepos):
    dstm = distance_array(ligpos, sitepos)
    # print('shape', dstm.shape)
    # shape: ligand atoms are rows, and for each there is an array of bindding site atoms
    by_row = np.min(dstm, axis=1)
    # print('min points', by_row.shape)
    # take the average position
    return np.mean(by_row)


avg_dst = get_avg_min_dst(ligand.positions, binding_site.positions)
print(f"Average distance from lig to binding site: {avg_dst} A")

# now use the binding site to track where the ligand resides
# for every lambda, combine the data from each replica
data = {}
for lamdir in glob.glob(r"complex/lambda_[0-1].[0-9]*"):
    lamdir_num = float(lamdir.split("_")[-1])
    print(f"Next lambda: {lamdir_num}")
    # record the means
    data[lamdir_num] = []

    for rep in glob.glob(os.path.join(lamdir, "rep[0-9]*")):
        print("rep", rep)
        # now open the trajectory with the simulation
        u = mda.Universe("complex/morph_solv.top", os.path.join(rep, "prod.dcd"))
        lig = u.select_atoms("resname ged", updating=True)
        binding_pocket = u.select_atoms(
            "not resname ged and around 5 resname ged", updating=True
        )

        # gather data after the first 2 ns (EQ time)
        avg_min_frame = []
        for ts in u.trajectory[100::10]:
            # calculate the distance between the ligand and the binding pocket
            avg_min_dst = get_avg_min_dst(lig.positions, binding_pocket.positions)
            avg_min_frame.append(avg_min_dst)
        data[lamdir_num].append(avg_min_frame)

# sort in the growing lambdas
ord_data = OrderedDict(sorted(data.items(), key=lambda x: x[0]))
print(ord_data.keys())

# for each lambda, display the distribution, but use the same x axis
# 13 subplot, 3 rows x 5 columns
counter = 1
for lam, reps in ord_data.items():
    plt.title(f"Lambda {lam}")
    plt.subplot(3, 5, counter)
    xfrom, xto = 2.3, 2.7
    # plt.xlim([xfrom, xto])
    for rep_points in reps:
        # plt.hist(rep_points, alpha=0.5, bins=np.linspace(xfrom, xto, 10))
        plt.hist(rep_points, alpha=0.5)  # , bins=np.linspace(xfrom, xto, 10))
    counter += 1

plt.tight_layout()
plt.show()
