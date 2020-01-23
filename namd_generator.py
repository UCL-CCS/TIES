"""
Load two ligands, run the topology superimposer, and then
using the results, generate the NAMD input files.
"""
from os import path
import MDAnalysis as mda
from topology_superimposer import get_atoms_bonds_from_ac, superimpose_topologies, assign_coords_from_pdb
import os
import json

def getSuptop():
    # test case based on agastya_dataset/mcl1/l18-l39

    # load the PDB files
    mda_left_lig = mda.Universe('agastya_dataset/mcl1/l18-l39/hybrid_par/init_l18.pdb')
    mda_right_lig = mda.Universe('agastya_dataset/mcl1/l18-l39/hybrid_par/final_l39.pdb')

    # fetch the charge information from the antechamber .ac files
    leftlig_atoms, leftlig_bonds = get_atoms_bonds_from_ac('agastya_dataset/mcl1/l18-l39/hybrid_par/init_l18.ac')
    rightlig_atoms, rightlig_bonds = get_atoms_bonds_from_ac('agastya_dataset/mcl1/l18-l39/hybrid_par/final_l39.ac')

    assign_coords_from_pdb(leftlig_atoms, mda_left_lig)
    assign_coords_from_pdb(rightlig_atoms, mda_right_lig)

    # assign
    # fixme - there must be a better way to match these
    c11 = None
    ligand1_nodes = {}
    for atomNode in leftlig_atoms:
        ligand1_nodes[atomNode.get_id()] = atomNode
        if atomNode.atomName == 'C11':
            c11 = atomNode
    for nfrom, nto in leftlig_bonds:
        ligand1_nodes[nfrom].bindTo(ligand1_nodes[nto])

    c33 = None
    ligand2_nodes = {}
    for atomNode in rightlig_atoms:
        ligand2_nodes[atomNode.get_id()] = atomNode
        if atomNode.atomName == 'C33':
            c33 = atomNode
    for nfrom, nto in rightlig_bonds:
        ligand2_nodes[nfrom].bindTo(ligand2_nodes[nto])

    suptops = superimpose_topologies(ligand1_nodes.values(), ligand2_nodes.values(),
                                     starting_node_pairs=[(c11, c33)])
    assert len(suptops) == 1
    suptop = suptops[0]

    present = [('C18', 'C40'), ('O2', 'O5')]
    for atomName1, atomname2 in present:
        assert suptop.contains_atomNamePair(atomName1, atomname2)

    removed_pairs = [('C14', 'C33'), ('C15', 'C34'), ('C16', 'C35'), ('C17', 'C36')]
    for atomName1, atomname2 in removed_pairs:
        assert not suptop.contains_atomNamePair(atomName1, atomname2)

    return suptop

suptop = getSuptop()
# we have the sup top, so we know which atoms should stay the same
# and which atoms should disappear and appear
# we have the coordinates too,
# so we have to generate a file and a script that is going to add a solvent,
# and we kind of have to do it at this stage, because we have to then add the information about which
# atoms appear and which atoms disappear
# here we have an issue, because we have to do both: generate the waters, and then correct the B field
# about which atoms disappear and which appear
# for this reason, it is best to store the appearance and disappearnce in a separate file,
# we're going to store it this way for now
"""
File: l1_l2.fep     # this is directional, from l1 to l2 
Matching: C1-C11, C2-C12, ...
Appearing: C3, C4
Disappearing: C19, C20
"""
# the file above can be used to update other files (with MDAnalysis ideally), or a simple script
output_dir = '/home/dresio/code/BAC2020/namd_tests/l18-l39'
# save the results
summary_filename = os.path.join(output_dir, 'l18_l39.fep')
# fixme - check if the file exists
with open(summary_filename, 'w') as FOUT:
    # use json format, only use atomNames
    data = {
            'matching': [{str(n1): str(n2)} for n1, n2 in suptop.matched_pairs],
            'appearing': list(map(str, suptop.get_appearing_atoms())),
            'disappearing': list(map(str, suptop.get_disappearing_atoms()))
            }
    FOUT.write(json.dumps(data, indent=4))

# we have a file on which we can now rely to correct the solvated file,
# todo copy the solvation script, as well as the ligand1, ligand2,

# todo copy the little python script that will be used to udpate the solvated script

# todo copy the EQ script that will be run

# todo copy the lambda .fep file to

"""
what should be the overall format? 
# todo - generate completely different directories with scripts with namd for each lambda

# todo - use sqlite to synchronise the execuation and managing of all of the simulations? (ie one major script)
for example, imagine a script that says "do_ties" which knows that there is 13 x 5 different directories
which have to be run each, and what it does, it goes into each of them, and schedules them, but 
it first checks where the simulation is by looking up its little .db file, 
ie lambda1.1 simulation has a db which is empty, so it schedules it to run, but lambda 1.2 has already finished, 
and that's written in its DB, whereas lambda 1.3 has status "submitted" and therefore does not need to be 
submitted again, whereas lambda 1.5 has status "running" which also can be ignored. etc etc
"""

print('hi')