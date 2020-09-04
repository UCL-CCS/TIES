"""
Find how far are the waters
"""
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np

root_dir = Path('/home/dresio/ucl/validation/resp_validation/mcl1/l1_l8/complex/lambda_0.00/rep1')
u = mda.Universe(root_dir / 'morph_solv.prmtop', str(root_dir / 'prod.dcd'))
print('Loaded', u.trajectory)

prot = u.select_atoms('protein')
print('Prot: ', len(prot))
lig = u.select_atoms('resname ged')
print('Lig', len(lig))
water_o = u.select_atoms('resname WAT')
water_dummies = u.select_atoms('resname WAT and name EPW')
print('Wat', len(water_o))

closest_wat_number = 10
# create the trajectory writers, which will
mda_writers = [mda.Writer(f'com-wat{watNum}.dcd', n_atoms=len(lig) + len(prot) + 3*watNum)
               for watNum in range(1, closest_wat_number + 1)]

cum_sums = np.zeros(len(water_o.residues))

every_nth_frame = 1
number_of_frames = 2
start_frame = 0 # 1 means the second frame
water_ids_per_frame = []
for ts in u.trajectory[start_frame:number_of_frames:every_nth_frame]:
    print("Time", ts.time)
    # get distances from water to X
    wat_prot = distance_array(water_o.positions, prot.positions, box=ts.dimensions)
    wat_lig = distance_array(water_o.positions, lig.positions, box=ts.dimensions)
    # take the shortest distances
    wat_prot_shortest = np.min(wat_prot, axis=1)
    wat_lig_shortest = np.min(wat_lig, axis=1)
    # get the minimum out of the 2 hydrogens and the oxygen
    wat_prot_res_min = np.min(wat_prot_shortest.reshape([7730, 3]), axis=1)
    wat_lig_res_min = np.min(wat_lig_shortest.reshape([7730, 3]), axis=1)
    # add up the shortest distances
    wat_others_sum = wat_prot_res_min + wat_lig_res_min
    cum_sums += wat_others_sum

    # PER FRAME SAVE
    # prepare for sorting
    chained = [(wid, dst) for wid, dst in zip(water_o.residues.resids, wat_others_sum)]
    # order by the distance
    chained_ord = sorted(chained, key=lambda x: x[1])
    # extract only the water IDs
    water_ids_per_frame.append([ts.frame] + [wid for wid, dst in chained_ord])

    # write the trajectories with the closest water molecules
    for watNum in range(0, closest_wat_number):
        # select the best waters, get their IDs and extract the indices
        com_wat_ids = [chained_ord[x][0] - 1 for x in range(watNum + 1)]
        com_wat = u.residues[com_wat_ids]
        # exclude the dummy atoms
        com_wat_no_dummies = com_wat.atoms.difference(water_dummies)
        mda_writers[watNum].write(prot + lig + com_wat_no_dummies.atoms)

# save the water molecules
np.savetxt('closest_water_ids.dat', water_ids_per_frame, fmt='%d', header='Frame closest_watID1 closest_watID2 ..')

# bind the two fields together for sorting
chained = [(wid, cum) for wid, cum in zip(water_o.residues.resids, cum_sums)]
# order by the cumulative sum
chained_ord = sorted(chained, key=lambda x:x[1])
to_save = np.array(chained_ord).T
print("Final shape", to_save.shape)
# sort according to the cumulative time
np.savetxt('cumulative_water_time.dat', to_save, fmt='%.1f')
