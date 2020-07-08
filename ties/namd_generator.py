"""
Load two ligands, run the topology superimposer, and then
using the results, generate the NAMD input files.

frcmod file format: http://ambermd.org/FileFormats.php#frcmod
"""
from os import path
import MDAnalysis as mda
from ties.topology_superimposer import get_atoms_bonds_from_mol2, superimpose_topologies, assign_coords_from_pdb
import os
import json

def getSuptop():
    # test case based on agastya_dataset/mcl1/l18-l39

    # fetch the charge information from the antechamber .ac files
    # ~/code/BAC2020/
    # fixme - you've switched to .mol2 format, so use mdanalysis for everything
    leftlig_atoms, leftlig_bonds, rightlig_atoms, rightlig_bonds, mda_l1, mda_l2 = \
        get_atoms_bonds_from_mol2('namd_tests/l18-l39/init_l18.mol2', 'namd_tests/l18-l39/final_l39.mol2')

    # assign
    # fixme - there must be a better way to match these
    c11 = None
    ligand1_nodes = {}
    for atomNode in leftlig_atoms:
        ligand1_nodes[atomNode.get_id()] = atomNode
        if atomNode.atomName == 'C11':
            c11 = atomNode
    for nfrom, nto, btype in leftlig_bonds:
        ligand1_nodes[nfrom].bindTo(ligand1_nodes[nto], btype)

    c33 = None
    ligand2_nodes = {}
    for atomNode in rightlig_atoms:
        ligand2_nodes[atomNode.get_id()] = atomNode
        if atomNode.atomName == 'C33':
            c33 = atomNode
    for nfrom, nto, btype in rightlig_bonds:
        ligand2_nodes[nfrom].bindTo(ligand2_nodes[nto], btype)

    suptops = superimpose_topologies(ligand1_nodes.values(), ligand2_nodes.values(),
                                     starting_node_pairs=[(c11, c33)])
    assert len(suptops) == 1
    suptop = suptops[0]

    present = [('C18', 'C40'), ('O2', 'O5')]
    for atomName1, atomname2 in present:
        assert suptop.contains_atomNamePair(atomName1, atomname2), f'({atomName1}, {atomname2}) not found'

    removed_pairs = [('C14', 'C33'), ('C15', 'C34'), ('C16', 'C35'), ('C17', 'C36')]
    for atomName1, atomname2 in removed_pairs:
        assert not suptop.contains_atomNamePair(atomName1, atomname2)

    return suptop, mda_l1, mda_l2

suptop, mda_l1, mda_l2 = getSuptop()
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
            'matching': {str(n1): str(n2) for n1, n2 in suptop.matched_pairs},
            'appearing': list(map(str, suptop.get_appearing_atoms())),
            'disappearing': list(map(str, suptop.get_disappearing_atoms()))
            }
    FOUT.write(json.dumps(data, indent=4))

# fixme - find another library that can handle writing to a PDB file
# todo
# we have a file on which we can now rely to correct the solvated file,
# save the ligand with all the appropriate atomic positions, write it using the pdb format
print('hi')
# pdb file format: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
# write a dual .pdb file
with open('namd_tests/l18-l39/l18l39.pdb', 'w') as FOUT:
    for atom in mda_l1.atoms:
        """
        There is only one forcefield which is shared across the two topologies. 
        Basically, we need to check whether the atom is in both topologies. 
        If that is the case, then the atom should have the same name, and therefore appear only once. 
        However, if there is a new atom, it should be specfically be outlined 
        that it is 1) new and 2) the right type
        """
        # write all the atoms if they are matched, that's the common part
        REMAINS = 0
        if suptop.contains_left_atomName(atom.name):
            line = f"ATOM  {atom.id:>5d} {atom.name:>4s} {atom.resname:>3s}  " \
                   f"{atom.resid:>4d}    " \
                   f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}" \
                   f"{1.0:>6.2f}{REMAINS:>6.2f}" + (' ' * 11) + \
                   '  ' + '  ' + '\n'
            FOUT.write(line)
        else:
            # this atom was not found, this means it disappears, so we should update the
            DISAPPEARING_ATOM = -1.0
            line = f"ATOM  {atom.id:>5d} {atom.name:>4s} {atom.resname:>3s}  " \
                   f"{atom.resid:>4d}    " \
                   f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}" \
                   f"{1.0:>6.2f}{DISAPPEARING_ATOM:>6.2f}" + \
                   (' ' * 11) + \
                   '  ' + '  ' + '\n'
            FOUT.write(line)
    # add atoms from the right topology,
    # which are going to be created
    for atom in mda_l2.atoms:
        if not suptop.contains_right_atomName(atom.name):
            APPEARING_ATOM = 1.0
            line = f"ATOM  {atom.id:>5d} {atom.name:>4s} {atom.resname:>3s}  " \
                   f"{atom.resid:>4d}    " \
                   f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}" \
                   f"{1.0:>6.2f}{APPEARING_ATOM:>6.2f}" + \
                   (' ' * 11) + \
                   '  ' + '  ' + '\n'
            FOUT.write(line)

# read the corresponding topology
ltop_file = '/home/dresio/code/BAC2020/namd_tests/l18-l39/init_l18.mol2'
rtop_file = '/home/dresio/code/BAC2020/namd_tests/l18-l39/final_l39.mol2'
top_merged = '/home/dresio/code/BAC2020/namd_tests/l18-l39/merged.mol2'

ltop = mda.Universe(ltop_file)
rtop = mda.Universe(rtop_file)
# todo check if the bonds information are there?
assert len(ltop.bonds) > 0 and len(rtop.bonds) > 0

# recreate the mol2 file that is merged and contains the correct atoms from both
def write_merged(suptop, merged_filename):
    # mol2 format: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
    with open(merged_filename, 'w') as FOUT:
        bonds = suptop.getDualTopologyBonds()

        FOUT.write('@<TRIPOS>MOLECULE ' + os.linesep)
        # name of the molecule
        FOUT.write('merged ' + os.linesep)
        # num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
        # fixme this is tricky
        FOUT.write(f'{suptop.getUnqiueAtomCount():d} '
                   f'{len(bonds):d}' + os.linesep)
        # mole type
        FOUT.write('SMALL ' + os.linesep)
        # charge_type
        FOUT.write('NO_CHARGES ' + os.linesep)
        FOUT.write(os.linesep)

        # write the atoms
        FOUT.write('@<TRIPOS>ATOM ' + os.linesep)
        # atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
        # e.g.
        #       1 O4           3.6010   -50.1310     7.2170 o          1 L39      -0.815300

        # so from the two topologies all the atoms are needed and they need to have a different atom_id
        # so we might need to name the atom_id for them, other details are however pretty much the same
        # the importance of atom_name is difficult to estimate

        # we are going to assign IDs in the superimposed topology in order to track which atoms have IDs
        # and which don't

        subst_id = 1    # resid basically
        # write all the atoms that were matched first with their IDs
        # prepare all the atoms, note that we use primarily the left ligand naming
        all_atoms = [left for left, right in suptop.matched_pairs] + suptop.get_unmatched_atoms()
        # reorder the list according to the ID
        all_atoms.sort(key=lambda atom: suptop.get_generated_atom_ID(atom))

        for atom in all_atoms:
            FOUT.write(f'{suptop.get_generated_atom_ID(atom)} {atom.atomName} '
                       f'{atom.position[0]:.4f} {atom.position[1]:.4f} {atom.position[2]:.4f} '
                       f'{atom.type.lower()} {subst_id} {atom.resname} {atom.charge:.6f} {os.linesep}')

        # for left_atom, _ in suptop.matched_pairs:
        #     # note that the atom id is the most important
        #     FOUT.write(f'{suptop.get_generated_atom_ID(left_atom)} {left_atom.atomName} '
        #                f'{left_atom.position[0]:.4f} {left_atom.position[1]:.4f} {left_atom.position[2]:.4f} '
        #                f'{left_atom.type} {subst_id} {left_atom.resname} {left_atom.charge} {os.linesep}')

        # write the IDs for the atoms which are appearing/disappearing
        # for unmatched in suptop.get_unmatched_atoms():
        #     FOUT.write(f'{suptop.get_generated_atom_ID(unmatched)} {unmatched.atomName} '
        #                f'{unmatched.position[0]:.4f} {unmatched.position[1]:.4f} {unmatched.position[2]:.4f} '
        #                f'{unmatched.type} {subst_id} {unmatched.resname} {unmatched.charge} {os.linesep}')

        FOUT.write(os.linesep)

        # write bonds
        FOUT.write('@<TRIPOS>BOND ' + os.linesep)

        # we have to list every bond:
        # 1) all the bonds between the paired atoms, so that part is easy
        # 2) bonds which link the disappearing atoms, and their connection to the paired atoms
        # 3) bonds which link the appearing atoms, and their connections to the paired atoms

        bond_counter = 1
        list(bonds)
        for bond_from_id, bond_to_id, bond_type in sorted(list(bonds)):
            # Bond Line Format:
            # bond_id origin_atom_id target_atom_id bond_type [status_bits]
            FOUT.write(f'{bond_counter} {bond_from_id} {bond_to_id} {bond_type}' + os.linesep)
            bond_counter += 1

write_merged(suptop, top_merged)


print('hi')
# todo write the topology file? example:
# /home/dresio/ucl/namd-ties-tutorial/FEP-tutorial-files/03.Mutating-tyrosine-into-alanine/tyr2ala.top

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