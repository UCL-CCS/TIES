import tempfile
from topology_superimposer import get_atoms_bonds_from_mol2, superimpose_topologies, general_atom_types2
import os
import sys
import json
import numpy as np
from collections import OrderedDict
import copy
import MDAnalysis as mda
import topology_superimposer
import shutil


def getSuptop(mol1, mol2, reference_match=None, force_mismatch=None,
              no_disjoint_components=True, net_charge_filter=True,
              ignore_charges_completely=False,
              use_general_type=True,
              ignore_bond_types=False,
              ignore_coords=False,
              left_coords_are_ref=True,
              align_molecules=True,
              use_only_gentype=False,
              check_atom_names_unique=True):
    # use mdanalysis to load the files
    # fixme - should not squash all messsages. For example, wrong type file should not be squashed
    leftlig_atoms, leftlig_bonds, rightlig_atoms, rightlig_bonds, mda_l1, mda_l2 = \
        get_atoms_bonds_from_mol2(mol1, mol2, use_general_type=use_general_type)

    # assign
    # fixme - Ideally I would reuse the mdanalysis data for this
    ligand1_nodes = {}
    for atomNode in leftlig_atoms:
        ligand1_nodes[atomNode.get_id()] = atomNode
        if reference_match is not None and atomNode.atomName == reference_match[0]:
            startLeft = atomNode
    for nfrom, nto, btype in leftlig_bonds:
        ligand1_nodes[nfrom].bindTo(ligand1_nodes[nto], btype)

    ligand2_nodes = {}
    for atomNode in rightlig_atoms:
        ligand2_nodes[atomNode.get_id()] = atomNode
        if reference_match is not None and atomNode.atomName == reference_match[1]:
            startRight = atomNode
    for nfrom, nto, btype in rightlig_bonds:
        ligand2_nodes[nfrom].bindTo(ligand2_nodes[nto], btype)

    if reference_match is None:
        starting_node_pairs = None
    else:
        starting_node_pairs = [(startLeft, startRight)]

    # fixme - simplify to only take the mdanalysis as input
    suptops = superimpose_topologies(ligand1_nodes.values(), ligand2_nodes.values(),
                                     starting_node_pairs=starting_node_pairs,
                                     force_mismatch=force_mismatch,
                                     no_disjoint_components=no_disjoint_components,
                                     net_charge_filter=net_charge_filter,
                                     ligandLmda=mda_l1, ligandRmda=mda_l2,
                                     ignore_charges_completely=ignore_charges_completely,
                                     ignore_bond_types=ignore_bond_types,
                                     ignore_coords=ignore_coords,
                                     left_coords_are_ref=left_coords_are_ref,
                                     align_molecules=align_molecules,
                                     use_general_type=use_general_type,
                                     use_only_gentype=use_only_gentype,
                                     check_atom_names_unique=check_atom_names_unique)
    assert len(suptops) == 1
    return suptops[0], mda_l1, mda_l2


def ensureUniqueAtomNames(left_mol2, right_mol2):
    """
    Ensure that each has a name that is unique to both files.

    # rename the ligand to ensure that no atom has the same name
    # name atoms using the first letter (C, N, ..) and count them
    # keep the names if possible (ie if they are already different)
    """
    # load both ligands
    left = topology_superimposer.load_mol2_wrapper(left_mol2)
    right = topology_superimposer.load_mol2_wrapper(right_mol2)

    # first, ensure that all the atom names are unique
    L_atom_names = [a.name for a in left.atoms]
    L_names_unique = len(set(L_atom_names)) == len(L_atom_names)
    L_correct_format = is_correct_atom_name_format(L_atom_names)

    if not L_names_unique or not L_correct_format:
        print('Renaming Left Molecule Atom Names (Because it is needed)')
        name_counter_L_nodes = rename_ligand(left.atoms)
        L_atom_names = [a.name for a in left.atoms]
    else:
        name_counter_L_nodes = get_atom_names_counter(left.atoms)

    R_atom_names = [a.name for a in right.atoms]
    R_names_unique = len(set(R_atom_names)) == len(R_atom_names)
    R_correct_format = is_correct_atom_name_format(R_atom_names)
    L_R_overlap = len(set(R_atom_names).intersection(set(L_atom_names))) > 0

    if not R_names_unique or not R_correct_format or L_R_overlap:
        print('Renaming Right Molecule Atom Names (Because it is needed)')
        rename_ligand(right.atoms, name_counter=name_counter_L_nodes)
    # each atom name is unique, fixme - this check does not apply anymore
    # ie it is fine for a molecule to use general type
    # assert len(set(R_atom_names)) == len(R_atom_names)

    # make a copy of the original files
    shutil.copy(left_mol2, str(left_mol2) + '_before_atomNameChange.mol2')
    shutil.copy(right_mol2, str(right_mol2) + '_before_atomNameChange.mol2')

    # overwrite the files
    left.atoms.write(left_mol2)
    right.atoms.write(right_mol2)


def is_correct_atom_name_format(names):
    """
    check if the atom format is like "C15", ie atom type followed by a number
    input atoms: MDAnalysis atoms
    """
    for name in names:
        afterLetters = [i for i, l in enumerate(name) if l.isalpha()][-1] + 1

        atom_name = name[:afterLetters]
        if len(atom_name) == 0:
            return False

        atom_number = name[afterLetters:]
        try:
            int(atom_number)
        except:
            return False

    return True


def rename_ligand(atoms, name_counter=None):
    """
    name_counter: a dictionary with atom as the key such as 'N', 'C', etc,
    the counter keeps track of the last used counter for each name.
    Empty means that the counting will start from 1.
    input atoms: mdanalysis atoms
    """
    if name_counter is None:
        name_counter = {}

    for atom in atoms:
        # get the first letters that is not a character
        afterLetters = [i for i, l in enumerate(atom.name) if l.isalpha()][-1] + 1

        atom_name = atom.name[:afterLetters]
        last_used_counter = name_counter.get(atom_name, 0)

        # rename
        last_used_counter += 1
        newAtomName = atom_name + str(last_used_counter)
        print(f'Renaming {atom.name} to {newAtomName}')
        atom.name = newAtomName

        # update the counter
        name_counter[atom_name] = last_used_counter

    return name_counter


def get_atom_names_counter(atoms):
    """
    name_counter: a dictionary with atom as the key such as 'N', 'C', etc,
    the counter keeps track of the last used counter for each name.
    Ie if there are C1, C2, C3, this will return {'C':3} as the last counter.
    """
    name_counter = {}

    for atom in atoms:
        # get the first letters that is not a character
        afterLetters = [i for i, l in enumerate(atom.name) if l.isalpha()][-1] + 1

        atom_name = atom.name[:afterLetters]
        atom_number = int(atom.name[afterLetters:])
        last_used_counter = name_counter.get(atom_name, 0)

        # update the counter
        name_counter[atom_name] = max(last_used_counter, atom_number)

    return name_counter


def save_superimposition_results(filepath, suptop):
    # fixme - check if the file exists
    with open(filepath, 'w') as FOUT:
        # use json format, only use atomNames
        data = {
            'matching': {str(n1): str(n2) for n1, n2 in suptop.matched_pairs},
            'appearing': list(map(str, suptop.get_appearing_atoms())),
            'disappearing': list(map(str, suptop.get_disappearing_atoms()))
        }
        FOUT.write(json.dumps(data, indent=4))


def write_dual_top_pdb(filepath, mda_l1, mda_l2, suptop):
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


def write_merged(suptop, merged_filename, use_left_charges=True, use_left_coords=True):
    # fixme - make this as a method of suptop as well
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

        # fixme - for writing, modify things to achieve the desired output
        # note - we are modifying in place our atoms
        for left, right in suptop.matched_pairs:
            print(f'Aligned {left.originalAtomName} id {left.atomId} with {right.originalAtomName} id {right.atomId}')
            if not use_left_charges:
                left.charge = right.charge
            if not use_left_coords:
                left.position = right.position

        subst_id = 1  # resid basically
        # write all the atoms that were matched first with their IDs
        # prepare all the atoms, note that we use primarily the left ligand naming
        all_atoms = [left for left, right in suptop.matched_pairs] + suptop.get_unmatched_atoms()
        unmatched_atoms = suptop.get_unmatched_atoms()
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
        for bond_from_id, bond_to_id, bond_type in sorted(list(bonds), key=lambda b: b[:2]):
            # Bond Line Format:
            # bond_id origin_atom_id target_atom_id bond_type [status_bits]
            FOUT.write(f'{bond_counter} {bond_from_id} {bond_to_id} {bond_type}' + os.linesep)
            bond_counter += 1


def parse_frcmod_sections(filename):
    """
    Copied from the previous TIES. It's simpler and this approach must be fine then.
    """
    frcmod_info = {}
    section = 'REMARK'

    with open(filename) as F:
        for line in F:
            start_line = line[0:9].strip()

            if start_line in ['MASS', 'BOND', 'IMPROPER',
                              'NONBON', 'ANGLE', 'DIHE']:
                section = start_line
                frcmod_info[section] = []
            elif line.strip() and section != 'REMARK':
                frcmod_info[section].append(line)

    return frcmod_info


def _merge_frcmod_section(ref_lines, other_lines):
    """
    A helper function for merging lines in .frcmod files.
    Note that the order has to be kept. This is because some lines need to follow other lines.
    In this case, we exclude lines that are present in ref_lines.
    fixme - since there are duplicate lines, we need to check the duplicates and their presence,
    """
    merged_section = copy.copy(ref_lines)
    for line in other_lines:
        if line not in ref_lines:
            merged_section.append(line)

    return merged_section


def join_frcmod_files2(filename1, filename2, output_filename):
    """
    Copied from the previous TIES. It's simpler and this approach must be fine then.
    """
    frcmod_info1 = parse_frcmod_sections(filename1)
    frcmod_info2 = parse_frcmod_sections(filename2)

    with open(output_filename, 'w') as FOUT:
        FOUT.write('merged frcmod\n')

        for section in ['MASS', 'BOND', 'ANGLE',
                        'DIHE', 'IMPROPER', 'NONBON']:
            section_lines = frcmod_info1[section] + frcmod_info2[section]
            FOUT.write('{0:s}\n'.format(section))
            for line in section_lines:
                FOUT.write('{0:s}'.format(line))
            FOUT.write('\n')

        FOUT.write('\n\n')

    return


def join_frcmod_files(f1, f2, output_filepath):
    """
    This implementation should be used. Switch to join_frcmod_files2.
    This version might be removed if the simple approach is fine.
    """
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
                for ritem in ritems:
                    if len(litems) > 0 or len(ritems) > 0:
                        if ritem not in joined[lname]:
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
            elif lname == 'NONBON':
                # if they're empty
                if not litems and not ritems:
                    continue

                raise Exception('Unimplemented')
            else:
                raise Exception('Unimplemented')
        return joined

    def write_frcmod(frcmod, filename):
        with open(filename, 'w') as FOUT:
            FOUT.write('GENERATED .frcmod by joining two .frcmod files' + os.linesep)
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
        # ignore water and ions and non-ligand resname
        # we only modify the protein, so ignore the ligand resnames
        # fixme .. why is it called mer, is it tleap?
        if atom.resname != 'mer':
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

    u.atoms.write(new_pdb_filename)  # , file_format='PDB') - fixme?


def get_ligand_resname(filename):
    lig_resnames = mda.Universe(filename).residues.resnames
    assert len(lig_resnames) == 1
    return lig_resnames


def get_morphed_ligand_resnames(filename):
    lig_resnames = mda.Universe(filename).residues.resnames
    # assert len(lig_resnames) == 2
    return lig_resnames


def get_PBC_coords(pdb_file):
    """
    Return [x, y, z]
    """
    u = mda.Universe(pdb_file)
    x = np.abs(max(u.atoms.positions[:, 0]) - min(u.atoms.positions[:, 0]))
    y = np.abs(max(u.atoms.positions[:, 1]) - min(u.atoms.positions[:, 1]))
    z = np.abs(max(u.atoms.positions[:, 2]) - min(u.atoms.positions[:, 2]))
    return (x, y, z)


def extract_PBC_oct_from_tleap_log(leap_log):
    """
    http://ambermd.org/namd/namd_amber.html
    Return the 9 numbers for the truncated octahedron unit cell in namd
    cellBasisVector1  d         0.0            0.0
    cellBasisVector2  (-1/3)*d (2/3)sqrt(2)*d  0.0
    cellBasisVector3  (-1/3)*d (-1/3)sqrt(2)*d (-1/3)sqrt(6)*d
    """
    leapl_log_lines = open(leap_log).readlines()
    line_to_extract = "Total bounding box for atom centers:"
    line_of_interest = list(filter(lambda l: line_to_extract in l, leapl_log_lines))
    d1, d2, d3 = line_of_interest[-1].split(line_to_extract)[1].split()
    d1, d2, d3 = float(d1), float(d2), float(d3)
    assert d1 == d2 == d3
    # scale the d since after minimisation the system turns out to be much smaller
    d = d1 * 0.8
    return {
        'cbv1': d, 'cbv2': 0, 'cbv3': 0,
        'cbv4': (-1/3.0)*d, 'cbv5': (2/3.0)*np.sqrt(2)*d, 'cbv6': 0,
        'cbv7': (-1/3.0)*d, 'cbv8': (-1/3.0)*np.sqrt(2)*d, 'cbv9': (-1/3)*np.sqrt(6)*d,
    }


def prepare_antechamber_parmchk2(source_script, target_script, net_charge):
    """
    Prepare the ambertools scripts.
    Particularly, create the scritp so that it has the net charge
    # fixme - run antechamber directly with the right settings from here?
    # fixme - check if antechamber has a python interface?
    """
    net_charge_set = open(source_script).read().format(net_charge=net_charge)
    open(target_script, 'w').write(net_charge_set)


def set_charges_from_ac(mol2_filename, ac_ref_filename):
    # ! the mol file will be overwritten
    print('Overwriting the mol2 file with charges from the ac file')
    # load the charges from the .ac file
    ac_atoms, _ = topology_superimposer.get_atoms_bonds_from_ac(ac_ref_filename)
    # load the .mol2 files with MDAnalysis and correct the charges
    mol2 = topology_superimposer.load_mol2_wrapper(mol2_filename)

    for mol2_atom in mol2.atoms:
        found_match = False
        for ac_atom in ac_atoms:
            if mol2_atom.name.upper() == ac_atom.atomName.upper():
                found_match = True
                mol2_atom.charge = ac_atom.charge
                break
        assert found_match, "Could not find the following atom in the AC: " + mol2_atom.name

    # update the mol2 file
    mol2.atoms.write(mol2_filename)


def set_charges_from_mol2(mol2_filename, mol2_ref_filename, by_atom_name=False, by_index=False, by_general_atom_type=False):
    # mol2_filename will be overwritten!
    print(f'Overwriting {mol2_filename} mol2 file with charges from {mol2_ref_filename} file')
    # load the ref charges
    ref_mol2 = topology_superimposer.load_mol2_wrapper(mol2_ref_filename)
    # load the .mol2 files with MDAnalysis and correct the charges
    mol2 = topology_superimposer.load_mol2_wrapper(mol2_filename)

    if by_atom_name and by_index:
        raise ValueError('Cannot have both. They are exclusive')
    elif not by_atom_name and not by_index:
        raise ValueError('Either option has to be selected.')

    # save the sum of charges before
    ref_sum_q =  sum(a.charge for a in ref_mol2.atoms)

    if by_atom_name:
        for mol2_atom in mol2.atoms:
            found_match = False
            for ref_atom in ref_mol2.atoms:
                if mol2_atom.name.upper() == ref_atom.name.upper():
                    if found_match == True:
                        raise Exception('AtomNames are not unique or do not match')
                    found_match = True
                    mol2_atom.charge = ref_atom.charge
            assert found_match, "Could not find the following atom: " + mol2_atom.name
    elif by_general_atom_type:
        for mol2_atom in mol2.atoms:
            found_match = False
            for ref_atom in ref_mol2.atoms:
                if general_atom_types2[mol2_atom.type.upper()] == general_atom_types2[ref_atom.type.upper()]:
                    if found_match == True:
                        raise Exception('AtomNames are not unique or do not match')
                    found_match = True
                    mol2_atom.charge = ref_atom.charge
            assert found_match, "Could not find the following atom in the AC: " + mol2_atom.name
    elif by_index:
        for mol2_atom, ref_atom in zip(mol2.atoms, ref_mol2.atoms):
                atype = general_atom_types2[mol2_atom.type.upper()]
                reftype = general_atom_types2[ref_atom.type.upper()]
                if atype != reftype:
                    raise Exception(f"The found general type {atype} does not equal to the reference type {reftype} ")

                mol2_atom.charge = ref_atom.charge

    assert ref_sum_q == sum(a.charge for a in mol2.atoms)
    # update the mol2 file
    mol2.atoms.write(mol2_filename)


def set_coor_from_ref(mol2_filename, coor_ref_filename, by_atom_name=False, by_index=False, by_general_atom_type=False):
    """
    fixme - add additional checks that check before and after averages to ensure that the copying was done right
    """
    # mol2_filename will be overwritten!
    print(f'Overwriting {mol2_filename} mol2 file with coordinates from {coor_ref_filename} file')
    # load the ref coordinates
    ref_mol2 = topology_superimposer.load_mol2_wrapper(coor_ref_filename)
    # load the .mol2 files with MDAnalysis and correct the charges
    mol2 = topology_superimposer.load_mol2_wrapper(mol2_filename)

    ref_pos_sum = np.sum(ref_mol2.atoms.positions)

    if by_atom_name and by_index:
        raise ValueError('Cannot have both. They are exclusive')
    elif not by_atom_name and not by_index:
        raise ValueError('Either option has to be selected.')

    if by_general_atom_type:
        for mol2_atom in mol2.atoms:
            found_match = False
            for ref_atom in ref_mol2.atoms:
                if general_atom_types2[mol2_atom.type.upper()] == general_atom_types2[ref_atom.type.upper()]:
                    found_match = True
                    mol2_atom.position = ref_atom.position
                    break
            assert found_match, "Could not find the following atom in the AC: " + mol2_atom.name
    if by_atom_name:
        for mol2_atom in mol2.atoms:
            found_match = False
            for ref_atom in ref_mol2.atoms:
                if mol2_atom.name.upper() == ref_atom.name.upper():
                    found_match = True
                    mol2_atom.position = ref_atom.position
                    break
            assert found_match, "Could not find the following atom in the AC: " + mol2_atom.name
    elif by_index:
        for mol2_atom, ref_atom in zip(mol2.atoms, ref_mol2.atoms):
                atype = general_atom_types2[mol2_atom.type.upper()]
                reftype = general_atom_types2[ref_atom.type.upper()]
                if atype != reftype:
                    raise Exception(f"The found general type {atype} does not equal to the reference type {reftype} ")

                mol2_atom.position = ref_atom.position

    if ref_pos_sum != np.sum(mol2.atoms.positions):
        raise Exception('Copying of the coordinates did not work correctly')

    # update the mol2 file
    mol2.atoms.write(mol2_filename)


def set_coor_from_ref_by_named_pairs(mol2_filename, coor_ref_filename, output_filename, left_right_pairs_filename):
    """
    Set coordinates but use atom names provided by the user.

    Example of the left_right_pairs_filename content:
    # flip the first ring
    # move the first c and its h
    C32 C18
    H34 C19
    # second pair
    C33 C17
    # the actual matching pair
    C31 C16
    H28 H11
    # the second matching pair
    C30 C15
    H29 H12
    #
    C35 C14
    # flip the other ring with itself
    C39 C36
    C36 C39
    H33 H30
    H30 H33
    C37 C38
    C38 C37
    H31 H32
    H32 H31
    """
    # fixme - check if the names are unique

    # parse "left_right_pairs_filename
    # format per line: leftatom_name right_atom_name
    lines = open(left_right_pairs_filename).read().split(os.linesep)
    left_right_pairs = (l.split() for l in lines if not l.strip().startswith('#'))

    # load the ref coordinates
    ref_mol2 = topology_superimposer.load_mol2_wrapper(coor_ref_filename)
    # load the .mol2 files with MDAnalysis and correct the charges
    static_mol2 = topology_superimposer.load_mol2_wrapper(mol2_filename)
    # this is being modified
    mod_mol2 = topology_superimposer.load_mol2_wrapper(mol2_filename)


    for pair in left_right_pairs:
        print('find pair', pair)
        new_pos = False
        for mol2_atom in static_mol2.atoms:
            # check if we are assigning from another molecule
            for ref_atom in ref_mol2.atoms:
                if mol2_atom.name.upper() == pair[0] and ref_atom.name.upper() == pair[1]:
                    new_pos = ref_atom.position
            # check if we are trying to assing coords from the same molecule
            for another_atom in static_mol2.atoms:
                if mol2_atom.name.upper() == pair[0] and another_atom.name.upper() == pair[1]:
                    new_pos = another_atom.position

        if new_pos is False:
            raise Exception("Could not find this pair: " + str(pair))

        # assign the position to the right atom
        # find pair[0]
        found = False
        for atom in mod_mol2.atoms:
            if atom.name.upper() == pair[0]:
                atom.position = new_pos
                found = True
                break
        assert found


    # update the mol2 file
    mod_mol2.atoms.write(output_filename)


def update_PBC_in_namd_input(namd_filename, new_pbc_box, structure_filename, constraint_lines=''):
    """
    fixme - rename this file since it generates the .eq files
    These are the lines we modify:
    cellBasisVector1	{cell_x}  0.000  0.000
    cellBasisVector2	 0.000  {cell_y}  0.000
    cellBasisVector3	 0.000  0.000 {cell_z}

    With x/y/z replacing the 3 values
    """
    assert len(new_pbc_box) == 3

    reformatted_namd_in = open(namd_filename).read().format(
        cell_x=new_pbc_box[0], cell_y=new_pbc_box[1], cell_z=new_pbc_box[2],
        constraints=constraint_lines, output='test_output', structure=structure_filename)

    # write to the file
    open(namd_filename, 'w').write(reformatted_namd_in)


def create_4_constraint_files(original_pdb, location):
    # Generate 4 constraint files and return the filenames
    """
coordinates  complex.pdb
constraints  on
consexp  2
consref  complex.pdb ;#need all positions
conskfile  constraint_f4.pdb
conskcol  B
    """

    def create_constraint(mda_universe, output, constraint):
        # for each atom, give the B column the right value
        for atom in mda_universe.atoms:
            # ignore water
            if atom.resname in ['WAT', 'Na+', 'TIP3W', 'TIP3']:
                continue

            # set each atom depending on whether it is a H or not
            if atom.name.upper().startswith('H'):
                atom.tempfactor = 0
            else:
                # restrain the heavy atom
                atom.tempfactor = constraint

        mda_universe.atoms.write(output)

    # create the 4 constraint files
    filenames = []
    u = mda.Universe(original_pdb)
    for i in range(1, 4 + 1):
        next_constraint_filename = location / f'constraint_f{i:d}.pdb'
        create_constraint(u, next_constraint_filename, i)
        filenames.append(next_constraint_filename)

    return filenames


def init_namd_file_min(from_dir, to_dir, filename, structure_name, pbc_box):
    min_namd_initialised = open(os.path.join(from_dir, filename)).read() \
        .format(structure_name=structure_name, **pbc_box)
    open(os.path.join(to_dir, filename), 'w').write(min_namd_initialised)


def generate_namd_eq(namd_eq, dst_dir, structure_name='sys_solv'):
    input_data = open(namd_eq).read()
    eq_namd_filenames = []
    for i in range(4):
        constraints = f"""
constraints  on
consexp  2
# use the same file for the position reference and the B column
consref  constraint_f{4 - i}.pdb ;#need all positions
conskfile  constraint_f{4 - i}.pdb
conskcol  B
        """
        if i == 0:
            prev_output = 'min_out'
        else:
            # take the output from the previous run
            prev_output = 'eq_out_%d' % i

        reformatted_namd_in = input_data.format(
            constraints=constraints, output='eq_out_%d' % (i + 1),
            prev_output=prev_output, structure_name=structure_name)
        # write to the file, start eq files count from 1
        next_eq_step_filename = dst_dir / ("eq_step%d.namd" % (i + 1))
        open(next_eq_step_filename, 'w').write(reformatted_namd_in)
        eq_namd_filenames.append(next_eq_step_filename)
    return eq_namd_filenames


def redistribute_charges(mda):
    """
    Calculate the original charges in the matched component.
    """


    return
