"""
Load two ligands, run the topology superimposer, and then
using the results, generate the NAMD input files.

todo - use ambertools python interface rather than ambertools directly,
       or find some other API bbblocks building blocks? or some other API,
       or create minimal your own
todo - adjust the input for surfsara for nwo
todo - estimate the size and legnth of the simulation, to adjust the number of cores/nodes
todo -

frcmod file format: `http://ambermd.org/FileFormats.php#frcmod`
"""
from topology_superimposer import get_atoms_bonds_from_mol2, superimpose_topologies, assign_coords_from_pdb
import os
import json
import numpy as np
import shutil
import sys
import subprocess
from collections import OrderedDict
import copy
import MDAnalysis as mda

def getSuptop(mol1, mol2):
    # use mdanalysis to load the files
    leftlig_atoms, leftlig_bonds, rightlig_atoms, rightlig_bonds, mda_l1, mda_l2 = \
        get_atoms_bonds_from_mol2(mol1, mol2)

    # assign
    # fixme - Ideally I would reuse the mdanalysis data for this
    startLeft = None
    ligand1_nodes = {}
    for atomNode in leftlig_atoms:
        ligand1_nodes[atomNode.get_id()] = atomNode
        if atomNode.atomName == 'O3':
            startLeft = atomNode
    for nfrom, nto, btype in leftlig_bonds:
        ligand1_nodes[nfrom].bindTo(ligand1_nodes[nto], btype)

    startRight = None
    ligand2_nodes = {}
    for atomNode in rightlig_atoms:
        ligand2_nodes[atomNode.get_id()] = atomNode
        if atomNode.atomName == 'O7':
            startRight = atomNode
    for nfrom, nto, btype in rightlig_bonds:
        ligand2_nodes[nfrom].bindTo(ligand2_nodes[nto], btype)

    suptops = superimpose_topologies(ligand1_nodes.values(), ligand2_nodes.values(),
                                     starting_node_pairs=[(startLeft, startRight)])
    assert len(suptops) == 1
    return suptops[0], mda_l1, mda_l2


def save_superimposition_results(filepath):
    # fixme - check if the file exists
    with open(filepath, 'w') as FOUT:
        # use json format, only use atomNames
        data = {
                'matching': {str(n1): str(n2) for n1, n2 in suptop.matched_pairs},
                'appearing': list(map(str, suptop.get_appearing_atoms())),
                'disappearing': list(map(str, suptop.get_disappearing_atoms()))
                }
        FOUT.write(json.dumps(data, indent=4))


def write_dual_top_pdb(filepath):
    # fixme - find another library that can handle writing to a PDB file, mdanalysis?
    # save the ligand with all the appropriate atomic positions, write it using the pdb format
    # pdb file format: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    # write a dual .pdb file
    with open(filepath, 'w') as FOUT:
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


def write_merged(suptop, merged_filename):
    # recreate the mol2 file that is merged and contains the correct atoms from both
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


def join_frcmod_files(f1, f2, output_filepath):
    # fixme - load f1 and f2

    def get_section(name, rlines):
        """
        Chips away from the lines until the section is ready

        fixme is there a .frcmod reader in ambertools?
        http://ambermd.org/FileFormats.php#frcmod
        """
        section_names = ['MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON']
        assert name in rlines.pop().strip()

        section = []
        while not (len(rlines) == 0 or any(rlines[-1].startswith(sname) for sname in section_names)):
            nextl = rlines.pop().strip()
            if nextl == '':
                continue
            # depending on the column name, parse differently
            if name == 'ANGLE':
                # e.g.
                # c -cc-na   86.700     123.270   same as c2-cc-na, penalty score=  2.6
                atom_types = nextl[:8]
                other = nextl[9:].split()[::-1]
                # The harmonic force constants for the angle "ITT"-"JTT"-
                #                     "KTT" in units of kcal/mol/(rad**2) (radians are the
                #                     traditional unit for angle parameters in force fields).
                harmonicForceConstant = float(other.pop())
                # TEQ        The equilibrium bond angle for the above angle in degrees.
                eq_bond_angle = float(other.pop())
                # the overall angle
                section.append([atom_types, harmonicForceConstant, eq_bond_angle])
            elif name == 'DIHE':
                # e.g.
                # ca-ca-cd-cc   1    0.505       180.000           2.000      same as c2-ce-ca-ca, penalty score=229.0
                atom_types = nextl[:11]
                other = nextl[11:].split()[::-1]
                """
                IDIVF      The factor by which the torsional barrier is divided.
                    Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
                    details. Basically, the actual torsional potential is

                           (PK/IDIVF) * (1 + cos(PN*phi - PHASE))

                 PK         The barrier height divided by a factor of 2.

                 PHASE      The phase shift angle in the torsional function.

                            The unit is degrees.

                 PN         The periodicity of the torsional barrier.
                            NOTE: If PN .lt. 0.0 then the torsional potential
                                  is assumed to have more than one term, and the
                                  values of the rest of the terms are read from the
                                  next cards until a positive PN is encountered.  The
                                  negative value of pn is used only for identifying
                                  the existence of the next term and only the
                                  absolute value of PN is kept.
                """
                IDIVF = float(other.pop())
                PK = float(other.pop())
                PHASE = float(other.pop())
                PN = float(other.pop())
                section.append([atom_types, IDIVF, PK, PHASE, PN])
            elif name == 'IMPROPER':
                # e.g.
                # cc-o -c -o          1.1          180.0         2.0          Using general improper torsional angle  X- o- c- o, penalty score=  3.0)
                # ...  IDIVF , PK , PHASE , PN
                atom_types = nextl[:11]
                other = nextl[11:].split()[::-1]
                # fixme - what is going on here? why is not generated this number?
                # IDIVF = float(other.pop())
                PK = float(other.pop())
                PHASE = float(other.pop())
                PN = float(other.pop())
                if PN < 0:
                    raise Exception('Unimplemented - ordering using with negative 0')
                section.append([atom_types, PK, PHASE, PN])
            else:
                section.append(nextl.split())
        return {name: section}

    def load_frcmod(filepath):
        # remark line
        rlines = open(filepath).readlines()[::-1]
        assert 'Remark' in rlines.pop()

        parsed = OrderedDict()
        for section_name in ['MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON']:
            parsed.update(get_section(section_name, rlines))

        return parsed

    def join_frcmod(left_frc, right_frc):
        joined = OrderedDict()
        for left, right in zip(left_frc.items(), right_frc.items()):
            lname, litems = left
            rname, ritems = right
            assert lname == rname

            joined[lname] = copy.copy(litems)

            if lname == 'MASS':
                if len(litems) > 0 or len(ritems) > 0:
                    raise Exception('Unimplemented')
            elif lname == 'BOND':
                if len(litems) > 0 or len(ritems) > 0:
                    raise Exception('Unimplemented')
            # ANGLE, e.g.
            # c -cc-na   86.700     123.270   same as c2-cc-na, penalty score=  2.6
            elif lname == 'ANGLE':
                for ritem in ritems:
                    # if the item is not in the litems, add it there
                    # extra the first three terms to determine if it is present
                    # fixme - note we are ignoring the "same as" note
                    if ritem not in joined[lname]:
                        joined[lname].append(ritem)
            elif lname == 'DIHE':
                for ritem in ritems:
                    if ritem not in joined[lname]:
                        joined[lname].append(ritem)
            elif lname == 'IMPROPER':
                for ritem in ritems:
                    if ritem not in joined[lname]:
                        joined[lname].append(ritem)
        return joined

    def write_frcmod(frcmod, filename):
        with open(filename, 'w') as FOUT:
            for sname, items in frcmod.items():
                FOUT.write(f'{sname}' + os.linesep)
                for item in items:
                    atom_types = item[0]
                    FOUT.write(atom_types)
                    numbers = ' \t'.join([str(n) for n in item[1:]])
                    FOUT.write(' \t' + numbers)
                    FOUT.write(os.linesep)
                # the ending line
                FOUT.write(os.linesep)
                print('hi')

    left_frc = load_frcmod(f1)
    right_frc = load_frcmod(f2)
    joined_frc = join_frcmod(left_frc, right_frc)
    write_frcmod(joined_frc, output_filepath)


def correct_fep_tempfactor(fep_json_filename, source_pdb_filename, new_pdb_filename):
    """
    fixme - this function does not need to use the file?
    we have the json information available here.
    """
    u = mda.Universe(source_pdb_filename)
    fep_meta_str = open(fep_json_filename).read()
    fep_meta = json.loads(fep_meta_str)

    left_matched = list(fep_meta['matching'].keys())
    appearing_atoms = fep_meta['appearing']
    disappearing_atoms = fep_meta['disappearing']

    # update the Temp column
    for atom in u.atoms:
        # ignore water and ions
        if atom.resname in ['WAT', 'Na+', 'Cl-']:
            continue

        # if the atom was "matched", meaning present in both ligands (left and right)
        # then ignore
        # note: we only use the left ligand
        if atom.name in left_matched:
            continue
        elif atom.name in appearing_atoms:
            # appearing atoms should
            atom.tempfactor = 1
        elif atom.name in disappearing_atoms:
            atom.tempfactor = -1

    u.atoms.write(new_pdb_filename)


workplace_root = '/home/dresio/code/BAC2020/namd_study/mcl_l8_l18'
# todo - check if there is left.pdb and right.pdb
if not os.path.isfile(os.path.join(workplace_root, 'left.pdb')):
    print('File left.pdb not found in', workplace_root)
    sys.exit(1)
elif not os.path.isfile(os.path.join(workplace_root, 'left.pdb')):
    print('File right.pdb not found in', workplace_root)
    sys.exit(1)
# copy the ambertools.sh for 1) creating .mol2 - antechamber, 2) optimising the structure with sqm
antechamber_sqm_script_name = 'assign_charge_parameters.sh'
script_dir = 'md_scripts'
# fixme - what are the net charges?
shutil.copy(os.path.join(script_dir, antechamber_sqm_script_name), workplace_root)
# execute the script (the script has to source amber.sh)
# do not do this if there is a .frcmod files
if not os.path.isfile(os.path.join(workplace_root, 'left.frcmod')):
    # fixme - add checks for other files
    output = subprocess.check_output(['sh', os.path.join(workplace_root, antechamber_sqm_script_name), 'left', 'right'])
    # todo - CHECK IF THE RESULTS ARE CORRECT

# load the files (.mol2) and superimpose the two topologies
# fixme - superimpose the molecules or stop relaying on RMSD info
# fixme - call any of the tools you have (antechamber, parmchk2)
suptop, mda_l1, mda_l2 = getSuptop(os.path.join(workplace_root, 'left.mol2'),
                                   os.path.join(workplace_root, 'right.mol2'))
# save the results of the topology superimposition as a json
top_sup_joint_meta = os.path.join(workplace_root, 'joint_meta_fep.json')
save_superimposition_results(top_sup_joint_meta)
write_dual_top_pdb(os.path.join(workplace_root, 'left_right.pdb'))
# save the merged topologies as a .mol2 file
top_merged_filename = os.path.join(workplace_root, 'merged.mol2')
write_merged(suptop, top_merged_filename)

# check if the .frcmod were generated
left_frcmod = os.path.join(workplace_root, 'left.frcmod')
right_frcmod = os.path.join(workplace_root, 'right.frcmod')
if not os.path.isfile(left_frcmod):
    sys.exit(5)
elif not os.path.isfile(right_frcmod):
    sys.exit(5)

# generate the joint .frcmod file
merged_frc_filename = os.path.join(workplace_root, 'merged.frcmod')
join_frcmod_files(left_frcmod, right_frcmod, merged_frc_filename)

# copy the solvate script for tleap
shutil.copy(os.path.join(script_dir, "run_tleap.sh"), workplace_root)
shutil.copy(os.path.join(script_dir, "leap.in"), workplace_root)
# solvate using AmberTools, copy leap.in and use tleap
output = subprocess.check_output(['sh', os.path.join(workplace_root, "run_tleap.sh")])
# this generates the "merged_solvated.pdb" which does not have .fep information in the .pdb tempfactor columns
merged_solvated = os.path.join(workplace_root, "merged_solvated.pdb")
merged_solvated_fep = os.path.join(workplace_root, "merged_solvated_fep.pdb")
correct_fep_tempfactor(top_sup_joint_meta, merged_solvated, merged_solvated_fep)


# Generate the directory structure for all the lambdas, and copy the files
for lambda_step in [0, 0.05] + list(np.linspace(0.1, 0.9, 9)) + [0.95, 1]:
    lambda_path = os.path.join(workplace_root, f'lambda_{lambda_step:.2f}')
    if not os.path.exists(lambda_path):
        os.makedirs(lambda_path)

    # for each lambda create 5 replicas
    for replica_no in range(1, 5 + 1):
        replica_dir = os.path.join(lambda_path, f'rep{replica_no}')
        if not os.path.exists(replica_dir):
            os.makedirs(replica_dir)

        # copy the files for each simulation
        # "coordinates" with the fep metadata
        shutil.copy(merged_solvated_fep, replica_dir)
        # fixme - where was this file generated
        shutil.copy(os.path.join(workplace_root, "merged_solvated.top"), replica_dir)

        # copy the NAMD files
        shutil.copy(os.path.join(script_dir, "eq.namd"), replica_dir)
        shutil.copy(os.path.join(script_dir, "prod.namd"), replica_dir)

        # copy the .tcl script that runs TIES
        shutil.copy(os.path.join(script_dir, "fep.tcl"), replica_dir)

        # set the lambda value for the directory
        with open(os.path.join(replica_dir, 'lambda'), 'w') as FOUT:
            FOUT.write(f'{lambda_step:.2f}')

        # copy the surfsara submit script - fixme - make this general
        shutil.copy(os.path.join(script_dir, "surfsara.sh"), os.path.join(replica_dir, 'submit.sh'))

        # copy the scheduler to the main directory
        shutil.copy(os.path.join(script_dir, "schedule_separately.py"), workplace_root)
        shutil.copy(os.path.join(script_dir, "check_namd_outputs.py"), workplace_root)

        # fixme States show the progress of the simulation.
print('hi')

# todo copy the solvation script, as well as the ligand1, ligand2,
# todo copy the little python script that will be used to udpate the solvated script
# todo copy the EQ script that will be run
# todo copy the lambda .fep file
"""
# todo - generate completely different directories with scripts with namd for each lambda
# todo - use sqlite to synchronise the execution and managing of all of the simulations? (ie one major script)
for example, imagine a script that says "do_ties" which knows that there is 13 x 5 different directories
which have to be run each, and what it does, it goes into each of them, and schedules them, but 
it first checks where the simulation is by looking up its little .db file, 
ie lambda1.1 simulation has a db which is empty, so it schedules it to run, but lambda 1.2 has already finished, 
and that's written in its DB, whereas lambda 1.3 has status "submitted" and therefore does not need to be 
submitted again, whereas lambda 1.5 has status "running" which also can be ignored. etc etc
"""

