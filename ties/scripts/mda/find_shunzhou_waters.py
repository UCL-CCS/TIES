"""
Find how far are the waters
"""
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import time
import argparse


def find_waters(topo, traj, prot_sel='protein', lig_sel='resname ged', wat_sel='resname WAT',
                output_file='closest_water_ids.dat', output_summary_file='cumulative_water_time.dat',
                traj_dcd_prefix='com_wat'):
    start_time = time.time()
    u = mda.Universe(topo, traj)
    print('Loaded', u.trajectory)

    prot = u.select_atoms(prot_sel)
    print('Prot: ', len(prot))
    lig = u.select_atoms(lig_sel)
    print('Lig', len(lig))
    water_o = u.select_atoms(wat_sel)
    print('Wat', len(water_o))

    closest_wat_number = 10
    # create the trajectory writers, which will
    mda_writers = [mda.Writer(f'{traj_dcd_prefix}{watNum}.dcd', n_atoms=len(lig) + len(prot) + 3*watNum)
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
        wat_prot_res_min = np.min(wat_prot_shortest.reshape([len(water_o.residues), 3]), axis=1)
        wat_lig_res_min = np.min(wat_lig_shortest.reshape([len(water_o.residues), 3]), axis=1)
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
            mda_writers[watNum].write(prot + lig + com_wat.atoms)

    # save the water molecules
    np.savetxt(output_file, water_ids_per_frame, fmt='%d', header='Frame closest_watID1 closest_watID2 ..')

    # bind the two fields together for sorting
    chained = [(wid, cum) for wid, cum in zip(water_o.residues.resids, cum_sums)]
    # order by the cumulative sum
    chained_ord = sorted(chained, key=lambda x:x[1])
    to_save = np.array(chained_ord).T
    print("Final shape", to_save.shape)
    # sort according to the cumulative time
    np.savetxt(output_summary_file, to_save, fmt='%.1f')
    print(f'Altogether the analysis took {time.time() - start_time:.0f} seconds')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Shunzhou water finding')
    parser.add_argument('-traj', type=str, help='Trajectory file', dest='traj')
    parser.add_argument('-top', type=str, help='Topology file', dest='top')
    args = parser.parse_args()

    #root_dir = Path('/home/dresio/ucl/validation/resp_validation/mcl1/l1_l8/complex/lambda_0.00/rep1')
    find_waters(args.top, args.traj)
    # find_waters('morph_solv.prmtop', 'prod.dcd')