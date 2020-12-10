#!/usr/bin/env python3
"""
Exposes a basic terminal interace to TIES 20.

Load two ligands, run the topology superimposer, and then using the results, generate the NAMD input files.
"""
import os
import shutil
import sys
import subprocess
import itertools
import csv
import argparse
from pathlib import Path, PurePosixPath

from ties.generator import *


def find_antechamber(args):
    # set up ambertools
    if args.ambertools_home is not None:
        # user manually specified the variable.
        ambertools_bin = args.ambertools_home / 'bin'
    elif os.getenv('AMBERHOME'):
        ambertools_bin = PurePosixPath(os.getenv('AMBERHOME')) / 'bin'
    elif os.getenv('AMBER_PREFIX'):
        ambertools_bin = PurePosixPath(os.getenv('AMBER_PREFIX')) / 'bin'
    else:
        print('Error: Cannot find ambertools. $AMBERHOME and $AMBER_PREFIX are empty')
        print('Option 1: source your ambertools script amber.sh')
        print('Option 2: specify manually the path to amberhome with -ambertools option')
        sys.exit()

    # fixme - test ambertools at this stage before proceeding
    return ambertools_bin


def convert_ac_to_mol2(filename, output_name, args, workplace_root, antechamber_dr='no'):
    """
    If the file is not a prepi file, this function does not do anything.
    Otherwise, antechamber is called to conver the .prepi file into a .mol2 file.

    @output_name: 'left' or 'right'

    Returns: the name of the original file, or of it was .prepi, a new filename with .mol2
    """
    base, ext = os.path.splitext(str(filename).lower())
    if ext != '.ac':
        return filename

    print('Will convert .ac to a .mol2 file')
    print('Searching for antechamber')
    ambertools_bin = find_antechamber(args)

    # a directory for the operation
    conversion_dir = f'converting_{output_name}_to_mol2'
    if not os.path.isdir(conversion_dir):
        os.makedirs(conversion_dir)

    # subprocess options for calling ambertools
    subprocess_kwargs = {
        "check": True, "text": True,
        "cwd": workplace_root / conversion_dir,
        "timeout": 30  # seconds
    }

    # prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
    print('Antechamber: converting the .ac to .mol2')
    new_filename = output_name + '_converted_from_ac.mol2'

    log_filename = workplace_root / conversion_dir / "antechamber_conversion.log"
    with open(log_filename, 'w') as LOG:
        try:
            subprocess.run([ambertools_bin / 'antechamber',
                                '-i', filename.resolve(), '-fi', 'ac',
                                '-o', new_filename, '-fo', 'mol2',
                                '-dr', antechamber_dr],
                                stdout=LOG, stderr=LOG,
                               **subprocess_kwargs)
        except subprocess.CalledProcessError as E:
            print('ERROR: An error occured during the antechamber conversion from .ac to .mol2 data type. ')
            print(f'ERROR: The output was saved in the directory: { workplace_root / conversion_dir}')
            print(f'ERROR: Please see the log file for the exact error information: {log_filename}')
            raise E

    # return a path to the new file as the input
    converted_file = workplace_root / conversion_dir / new_filename
    print(f'Converted .ac file to .mol2. The location of the new file: {converted_file}')
    return converted_file


def str2bool(v):
    "ArgumentParser tool to figure out the bool value"
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


class Ligand:
    """
    The ligand helper class

    TODO - use a general conf class/dict to find out about ambertools
    """
    LIG_COUNTER = 0

    UNIQ_ATOM_NAME_DIR = 'unique_atom_names'
    FRCMOD_DIR = 'ligand_frcmods'

    def __init__(self, ligand, workplace_root):
        # save workplace root
        self.workplace_root = Path(workplace_root)
        self.original_input = self.workplace_root / ligand

        # check if the input files exist
        if not self.original_input.is_file():
            print(f'ERROR: Ligand file {self.original_input} not found')
            sys.exit(1)

        # last used representative Path file
        self.current = self.original_input

        Ligand.LIG_COUNTER += 1
        # internal name without an extension
        self.internal_name = f'ligand{Ligand.LIG_COUNTER:d}'

        self.frcmod = None

    def make_atom_names_unique(self):
        """
        Ensure that each atom has a unique name.

        rename the ligand to ensure that no atom has the same atom name
        using the first letter (C, N, ..) and count them
        keep the names if possible (ie if they are already different)

        @parameter save_update: if the path is provided, the updated file
            will be saved with the unique names and a handle to the new file (MDAnalysis universe)
            will be returned.
        """

        # save the output here
        os.makedirs(self.workplace_root / Ligand.UNIQ_ATOM_NAME_DIR, exist_ok=True)

        # load the ligand with MDAnalysis
        ligand_universe = load_mol2_wrapper(self.current)

        # ensure that all the atom names are unique
        atom_names = [a.name for a in ligand_universe.atoms]
        names_unique = len(set(atom_names)) == len(atom_names)

        if not names_unique or not are_correct_names(atom_names):
            print(f'Atom names in your molecule ({self.original_input}) are either not unique '
                  f'or do not follow NameDigit format (e.g. C15). Renaming')
            rename_ligand(ligand_universe.atoms)

        ligand_with_uniq_atom_names = self.workplace_root / Ligand.UNIQ_ATOM_NAME_DIR / \
                                      (self.internal_name + self.current.suffix)
        ligand_universe.atoms.write(ligand_with_uniq_atom_names)

        self.ligand_with_uniq_atom_names = ligand_with_uniq_atom_names
        # this object is now represented by the updated ligand
        self.current = ligand_with_uniq_atom_names

    # make this into a python property
    def suffix(self):
        return self.current.suffix.lower()

    def antechamber_prepare_mol2(self, ambertools_bin, atom_type, net_charge, antechamber_dr, antechamber_charge_type):
        """
        Convert the files into .mol2 files. Generate BCC charges if needed.
        A helper function that calls antechamber and ensures that the log is kept.
        The default behaviour is to keep the results in the file.

        # antechamber note:
        # charge type -c is not used if user provided prefer to use their charges
        """
        print('Antechamber: converting to .mol2 and generating charges if necessary')
        if not antechamber_charge_type:
            print('Antechamber: Ignoring atom charges. The user-provided atom charges will be used. ')
        else:
            print('Antechamber: Generating BCC charges')

        mol2_cwd = self.workplace_root / 'mol2_prep' / self.internal_name

        # prepare the directory
        if not mol2_cwd.is_dir():
            mol2_cwd.mkdir(parents=True, exist_ok=True)

        subprocess_kwargs = {
            "check": True, "text": True,
            "timeout": 60 * 30  # 30 minutes
        }

        # fixme - does not throw errors? try to make it throw an error
        # do not redo if the target file exists
        mol2_target = mol2_cwd / f'{self.internal_name}.mol2'
        if not (mol2_target).is_file():
            log_filename = mol2_cwd / "antechamber.log"
            with open(log_filename, 'w') as LOG:
                try:
                    subprocess.run([ambertools_bin / 'antechamber',
                                    '-i', self.current, '-fi', self.current.suffix[1:],
                                    '-o', mol2_target, '-fo', 'mol2',
                                    '-at', atom_type, '-nc', str(net_charge),
                                    '-dr', antechamber_dr] + antechamber_charge_type,
                                   cwd=mol2_cwd,
                                   stdout=LOG, stderr=LOG,
                                   **subprocess_kwargs)
                except subprocess.CalledProcessError as E:
                    print('ERROR: occured when creating the input .mol2 file with antechamber. ')
                    print(f'ERROR: The output was saved in the directory: {mol2_cwd}')
                    print(f'ERROR: can be found in the file: {log_filename}')
                    raise E
            print(f'Converted {self.original_input} into .mol2, Log: {log_filename}')
        else:
            print(f'File {mol2_target} already exists. Skipping. ')

            self.antechamber_mol2 = mol2_target
            self.current = mol2_target

            # remove any DUMMY DU atoms in the .mol2 atoms
            self.removeDU_atoms()

    def removeDU_atoms(self):
        """
        Ambertools antechamber creates sometimes DU dummy atoms.
        These are not created when BCC charges are computed from scratch.
        They are only created if you reuse existing charges.
        They appear to be a side effect. We remove the dummy atoms therefore.
        """
        mol2_u = load_mol2_wrapper(self.current)
        # check if there are any DU atoms
        has_DU = any(a.type == 'DU' for a in mol2_u.atoms)
        if not has_DU:
            return

        # make a backup copy before (to simplify naming)
        shutil.move(self.current, self.current.parent / ('lig.beforeRemovingDU' + self.current.suffix))

        # remove DU type atoms and save the file
        mol2_u.select_atoms('not type DU').atoms.write(self.current)
        print('Removed dummy atoms with type "DU"')

    def generate_frcmod(self, parmchk2, atom_type):
        # fixme - work on the file handles instaed of the constant stitching
        print(f'Parmchk2: generate the .frcmod for {self.internal_name}.mol2')

        # prepare cwd
        cwd = self.workplace_root / Ligand.FRCMOD_DIR / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        log_filename = cwd / "parmchk2.log"
        target_frcmod = f'{self.internal_name}.frcmod'
        with open(log_filename, 'w') as LOG:
            try:
                subprocess.run([parmchk2,
                                '-i', self.current,
                                '-o', target_frcmod,
                                '-f', 'mol2',
                                '-s', atom_type],
                               stdout=LOG, stderr=LOG,
                               check= True, text=True,
                               cwd= cwd, timeout=20,  # 20 seconds
                                )
            except subprocess.CalledProcessError as E:
                print('ERROR: An error occured during the antechamber conversion from .ac to .mol2 data type. ')
                print(f'ERROR: The output was saved in the directory: {cwd}')
                print(f'ERROR: Please see the log file for the exact error information: {log_filename}')
                raise E

        print(f'Parmchk2: created .frcmod: {target_frcmod}')
        self.frcmod = cwd / target_frcmod

class Morph():
    """
    A convenience class to help organise hybrids/morphs.

    Note that Morph actually in a way duplicates the concept of a SupTop.
    This distinction could be removed easily in the future.
    However, way to simplify SupTop might be worth some effort.
    """
    FRCMOD_DIR = Path('morph_frcmods')
    FRCMOD_TEST_DIR = FRCMOD_DIR / 'tests'
    UNIQUE_ATOM_NAMES = Path('morph_unique_atom_names')

    def __init__(self, ligA, ligZ, workplace_root):
        self.ligA = ligA
        self.ligZ = ligZ
        self.workplace_root = workplace_root

        self.internal_name = f'{ligA.internal_name}_{ligZ.internal_name}'
        self.mol2 = None
        self.pdb = None
        self.summary = None
        self.suptop = None
        self.mda_l1 = None
        self.mda_l2 = None

    def set_suptop(self, suptop, mda_l1, mda_l2):
        self.suptop = suptop
        self.mda_l1 = mda_l1
        self.mda_l2 = mda_l2

    def unique_atomres_names(self):
        """
        Ensure that each has a name that is unique to both ligands.

        Rename the ligand to ensure that no atom has the same name
        name atoms using the first letter (C, N, ..) and count them
        keep the names if possible (ie if they are already different)

        Resnames are set to "INI" and "FIN", this is useful for the hybrid dual topology
        """
        self.ligA.make_atom_names_unique()
        self.ligZ.make_atom_names_unique()

        # load both ligands
        left = load_mol2_wrapper(self.ligA.current)
        right = load_mol2_wrapper(self.ligZ.current)

        ligands_atom_names_overlap = len(set(right.atoms.names).intersection(set(left.atoms.names))) > 0

        right_renamed = False
        if ligands_atom_names_overlap:
            print(f'Renaming right molecule ({self.ligZ.internal_name}) atom names because they are not unique')
            name_counter_L_nodes = get_atom_names_counter(left.atoms)
            rename_ligand(right.atoms, name_counter=name_counter_L_nodes)
            right_renamed = True

        # rename the resnames to INI and FIN
        left.residues.resnames = ['INI']
        right.residues.resnames = ['FIN']

        # prepare the destination directory
        cwd = self.workplace_root / Morph.UNIQUE_ATOM_NAMES / f'{self.ligA.internal_name}_{self.ligZ.internal_name}'
        cwd.mkdir(parents=True, exist_ok=True)

        # save the updated atom names
        self.renamed_ligA = cwd / (self.ligA.internal_name + '.mol2')
        left.atoms.write(self.renamed_ligA)
        self.renamed_ligZ = cwd / (self.ligZ.internal_name + '.mol2')
        right.atoms.write(self.renamed_ligZ)

    def write_superimposition_json(self):
        """
        Writes a simple .json file with a summary of which atoms are classified as appearing, disappearing.
        """

        # store at the root for now
        matching_json = self.workplace_root / f'meta_fep_{self.ligA.internal_name}_{self.ligZ.internal_name}.json'

        with open(matching_json, 'w') as FOUT:
            # use json format, only use atomNames
            app, dis = self.suptop.get_single_topology_app()
            summary = {
                # the dual topology information
                'matched': {str(n1): str(n2) for n1, n2 in self.suptop.matched_pairs},
                'appearing': list(map(str, self.suptop.get_appearing_atoms())),
                'disappearing': list(map(str, self.suptop.get_disappearing_atoms())),
                # single topology information
                'single_top_matched': {str(n1): str(n2) for n1, n2 in self.suptop.get_single_topology_region()},
                # NAMD hybrid single-dual topology info
                'single_top_appearing': list(map(str, app)),
                'single_top_disappearing': list(map(str, dis)),
            }
            FOUT.write(json.dumps(summary, indent=4))

        self.summary = summary

    def write_morph_pdb(self, hybrid_single_dual_top):
        morph_pdb_path = self.workplace_root / f'{self.ligA.internal_name}_{self.ligZ.internal_name}_morph.pdb'

        # def write_morph_top_pdb(filepath, mda_l1, mda_l2, suptop, hybrid_single_dual_top=False):
        if hybrid_single_dual_top:
            # the NAMD hybrid single dual topology
            # rename the ligand on the left to INI
            # and the ligand on the right to END

            # make a copy of the suptop here to ensure that the modifications won't affect it
            st = copy.copy(self.suptop)

            # first, set all the matched pairs to -2 and 2 (single topology)
            # regardless of how they were mismatched
            raise NotImplementedError('Cannot yet write hybrid single dual topology .pdb file')

            # then, set the different atoms to -1 and 1 (dual topology)

            # save in a single PDB file
            # Note that the atoms from left to right
            # in the single topology region have to
            # be separated by the same number
            # fixme - make a check for that
            return
        # fixme - find another library that can handle writing to a PDB file, MDAnalysis
        # save the ligand with all the appropriate atomic positions, write it using the pdb format
        # pdb file format: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        # write a dual .pdb file
        with open(morph_pdb_path, 'w') as FOUT:
            for atom in self.mda_l1.atoms:
                """
                There is only one forcefield which is shared across the two topologies. 
                Basically, we need to check whether the atom is in both topologies. 
                If that is the case, then the atom should have the same name, and therefore appear only once. 
                However, if there is a new atom, it should be specfically be outlined 
                that it is 1) new and 2) the right type
                """
                # write all the atoms if they are matched, that's the common part
                REMAINS = 0
                if self.suptop.contains_left_atom_name(atom.name):
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
            for atom in self.mda_l2.atoms:
                if not self.suptop.contains_right_atom_name(atom.name):
                    APPEARING_ATOM = 1.0
                    line = f"ATOM  {atom.id:>5d} {atom.name:>4s} {atom.resname:>3s}  " \
                           f"{atom.resid:>4d}    " \
                           f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}" \
                           f"{1.0:>6.2f}{APPEARING_ATOM:>6.2f}" + \
                           (' ' * 11) + \
                           '  ' + '  ' + '\n'
                    FOUT.write(line)

    def write_hybrid_mol2(self, use_left_charges=True, use_left_coords=True):
        hybrid_mol2 = self.workplace_root / f'{self.ligA.internal_name}_{self.ligZ.internal_name}_morph.mol2'

        # fixme - make this as a method of suptop as well
        # recreate the mol2 file that is merged and contains the correct atoms from both
        # mol2 format: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
        # fixme - build this molecule using the MDAnalysis builder instead of the current approach
        # where the formatting is done manually
        with open(hybrid_mol2, 'w') as FOUT:
            bonds = self.suptop.get_dual_topology_bonds()

            FOUT.write('@<TRIPOS>MOLECULE ' + os.linesep)
            # name of the molecule
            FOUT.write('merged ' + os.linesep)
            # num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
            # fixme this is tricky
            FOUT.write(f'{self.suptop.get_unique_atom_count():d} '
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
            for left, right in self.suptop.matched_pairs:
                print(
                    f'Aligned {left.originalAtomName} id {left.atomId} with {right.originalAtomName} id {right.atomId}')
                if not use_left_charges:
                    left.charge = right.charge
                if not use_left_coords:
                    left.position = right.position

            subst_id = 1  # resid basically
            # write all the atoms that were matched first with their IDs
            # prepare all the atoms, note that we use primarily the left ligand naming
            all_atoms = [left for left, right in self.suptop.matched_pairs] + self.suptop.get_unmatched_atoms()
            unmatched_atoms = self.suptop.get_unmatched_atoms()
            # reorder the list according to the ID
            all_atoms.sort(key=lambda atom: self.suptop.get_generated_atom_id(atom))

            for atom in all_atoms:
                FOUT.write(f'{self.suptop.get_generated_atom_id(atom)} {atom.name} '
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

        self.mol2 = hybrid_mol2

    def join_frcmod_files(self, ambertools_bin, ambertools_script_dir, protein_ff, ligand_ff):
        """
        Copied from the previous TIES. It's simpler and this appears to work fine.
        I tested the duplication and that seemed to have no effect on the final results.

        Note that we are also testing the .frcmod here with the specific user FF.
        """
        frcmod_info1 = parse_frcmod_sections(self.ligA.frcmod)
        frcmod_info2 = parse_frcmod_sections(self.ligZ.frcmod)

        # prep directory
        cwd = Path(Morph.FRCMOD_DIR)
        if not cwd.is_dir():
            cwd.mkdir(exist_ok=True)

        # fixme: use the provided cwd here, otherwise this will not work if the wrong cwd is used
        # have some conf module instead of this
        morph_frcmod = cwd / f'{self.ligA.internal_name}_{self.ligZ.internal_name}_morph.frcmod'
        with open(morph_frcmod, 'w') as FOUT:
            FOUT.write('merged frcmod\n')

            for section in ['MASS', 'BOND', 'ANGLE',
                            'DIHE', 'IMPROPER', 'NONBON']:
                section_lines = frcmod_info1[section] + frcmod_info2[section]
                FOUT.write('{0:s}\n'.format(section))
                for line in section_lines:
                    FOUT.write('{0:s}'.format(line))
                FOUT.write('\n')

            FOUT.write('\n\n')

        # this is our current frcmod file
        self.frcmod = morph_frcmod

        # as part of the .frcmod writing
        # insert dummy angles/dihedrals if a morph .frcmod requires
        # new terms between the appearing/disappearing atoms
        self._check_hybrid_frcmod(ambertools_bin, ambertools_script_dir, protein_ff, ligand_ff)

    def _check_hybrid_frcmod(self, ambertools_bin, ambertools_script_dir, protein_ff, ligand_ff):
        """
        Previous code: https://github.com/UCL-CCS/BacScratch/blob/master/agastya/ties_hybrid_topology_creator/output.py
        Check that the output library can be used to create a valid amber topology.
        Add missing terms with no force to pass the topology creation.
        Returns the corrected .frcmod content, otherwise throws an exception.
        """
        # prepare the working directory
        cwd = self.workplace_root / self.FRCMOD_TEST_DIR / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        modified_hybrid_frcmod = cwd / f'{self.internal_name}_corrected.frcmod'

        # prepare tleap input
        leap_in_test = 'leap_test_morph.in'
        leap_in_conf = open(ambertools_script_dir / leap_in_test).read()
        open(cwd / leap_in_test, 'w').write(leap_in_conf.format(
                                mol2=os.path.relpath(self.mol2, cwd),
                                frcmod=os.path.relpath(self.frcmod, cwd),
                                protein_ff=protein_ff, ligand_ff=ligand_ff))

        # attempt generating the .top
        print('Create amber7 topology .top')
        tleap_process = subprocess.run([ambertools_bin / 'tleap', '-s', '-f', leap_in_test],
                                       cwd=cwd, text=True, timeout=20,
                                       capture_output=True, check=True)

        # save stdout and stderr
        open(cwd / 'tleap_scan_check.log', 'w').write(tleap_process.stdout + tleap_process.stderr)

        if 'Errors = 0' in tleap_process.stdout:
            print('Test hybrid .frcmod: OK, no dummy angle/dihedrals inserted.')
            return

        # extract the missing angles/dihedrals
        missing_angles = []
        missing_dihedrals = []
        for line in tleap_process.stdout.splitlines():
            if "Could not find angle parameter:" in line:
                cols = line.split(':')
                angle = cols[-1].strip()
                if angle not in missing_angles:
                    missing_angles.append(angle)
            elif "No torsion terms for" in line:
                cols = line.split()
                torsion = cols[-1].strip()
                if torsion not in missing_dihedrals:
                    missing_dihedrals.append(torsion)

        if missing_angles or missing_dihedrals:
            print('WARNING: Adding dummy dihedrals/angles to frcmod to generate .top')
            frcmod_lines = open(modified_hybrid_frcmod).readlines()
            # overwriting the .frcmod with dummy angles/dihedrals
            with open(modified_hybrid_frcmod, 'w') as NEW_FRCMOD:
                for line in frcmod_lines:
                    NEW_FRCMOD.write(line)
                    if 'ANGLE' in line:
                        for angle in missing_angles:
                            dummy_angle = f'{angle:<14}0  120.010  \t\t# Dummy angle\n'
                            NEW_FRCMOD.write(dummy_angle)
                            print(f'Added dummy angle: "{dummy_angle}"')
                    if 'DIHE' in line:
                        for dihedral in missing_dihedrals:
                            dummy_dihedral = f'{dihedral:<14}1  0.00  180.000  2.000   \t\t# Dummy dihedrals\n'
                            NEW_FRCMOD.write(dummy_dihedral)
                            print(f'Added dummy dihedral: "{dummy_dihedral}"')

            # verify that adding the dummy angles/dihedrals worked
            tleap_process = subprocess.run([ambertools_bin / 'tleap', '-s', '-f', leap_in_test],
                                           cwd=cwd, text=True, timeout=60 * 10, capture_output=True, check=True)

            if not "Errors = 0" in tleap_process.stdout:
                raise Exception('ERROR: Could not generate the .top file after adding dummy angles/dihedrals')

        print('Morph .frcmod after the insertion of dummy angle/dihedrals: OK')
        # set this .frcmod as the correct one now,
        self.frcmod_before_correction = self.frcmod
        self.frcmod = modified_hybrid_frcmod


def command_line_script():
    parser = argparse.ArgumentParser(description='TIES 20')
    parser.add_argument('action', metavar='command', type=str,
                        help='Action to be performed. E.g. "ties rename, ties create, ties .." ')
    parser.add_argument('-l', '--ligands', metavar='Ligands ', dest='ligands',
                        nargs='+',
                        #type=Path, required=False,
                        help='A list of the ligands in the right order. When two ligands are given they would be '
                             'treated as Left and Right (Right-Left). '
                             'If more ligands are given, Lead Optmisation Mapping (like LOMAP) will be used.')
    parser.add_argument('-cwd', metavar='Current_Working_Directory', dest='cwd',
                        type=Path, required=False,
                        help='If not provided, current directory is used. ')
    parser.add_argument('-nc', '--net-charge', metavar='Right_Ligand_File', dest='net_charge',
                        type=int, required=False,
                        help='The right ligand filename')
    parser.add_argument('-p', '--protein', metavar='Protein_file', dest='protein',
                        type=Path, required=False,
                        help='The protein file')
    parser.add_argument('-qtol', '--q-atom-tolerance', metavar='Charge_tolerance', dest='qtol',
                        type=float, required=False, default=0.1,
                        help='The maximum difference in charge between any two atoms in electrons. ')
    parser.add_argument('-netqtol', '--q-net-tolerance', metavar='Net_Charge_tolerance', dest='net_charge_threshold',
                        type=float, required=False, default=0.1,
                        help='The maximum difference in charge between the two ligands. ')
    parser.add_argument('-use-provided-q', '--use-user-q', metavar='True_False', dest='ligands_have_q',
                        type=str2bool, required=False,
                        help='Use charges provided with the ligand files (with .mol2). '
                             'Default: If .mol2 is given, using the given charges will be attempted. '
                             'Default: If .pdb is given, then BCC charges are computed with antechamber. ')
    parser.add_argument('-align-mcs', '--align-ligands-mcs', metavar='Align_Ligands_MCS', dest='align_mcs',
                        type=str2bool, required=False, default=False,
                        help='Align the right ligand to the left using the determined common component')
    parser.add_argument('-ant-dr', '--antechamber-dr', metavar='AntechamberDrMode', dest='antechamber_dr',
                        type=str2bool, required=False, default=True,
                        help='Antechamber dr is turned off by default. It is playing up with .mol2 files. '
                             'Please ensure that you input is valid, or turn on the antechamber dr. ')
    parser.add_argument('-ambertools', '--ambertools-home', metavar='Ambertools_path', dest='ambertools_home',
                        type=Path, required=False,
                        help='Path to the home directory of ambertools '
                             '(the one that contains "bin" directory and others)')
    parser.add_argument('-match', '--manual-match', metavar='ManualMatch', dest='manual_match_file',
                        type=Path, required=False,
                        help='A path to a file that contains atom-pairs (one per line).'
                             'Each pair should be separated by dash "-" character.'
                             'For example a line with CL2-BR2 means that CL2 atom will be transformed to BR2 atom.')
    parser.add_argument('-mismatch', '--manual-mismatch', metavar='ManualMismatch', dest='manual_mismatch_file',
                        type=Path, required=False,
                        help='A path to a file that contains atom-pairs (one per line).'
                             'Each pair should be separated by dash "-" character.'
                             'For example a line with CL2-BR2 means that CL2 atom cannot be matched to BR2 atom.'
                             'Note that this option might have undesired consequences. ')
    parser.add_argument('-redist-q', '--redistribute-q-over-unmatched', dest='redistribute_charges_over_unmatched',
                        type=str2bool, required=False, default=True,
                        help='Averaging the charges in the matched areas changes the overall molecule charge slightly. '
                             'Redistribute the lost/gained charges over the unmatched area to make the two molecules '
                             'equal. ')
    parser.add_argument('-hybrid-top', '--hybrid-singe-dual-top', dest='hybrid_single_dual_top',
                        type=str2bool, required=False, default=False,
                        help='Use hybrid single-dual topology in NAMD. See NAMD manual and '
                             'https://pubs.acs.org/doi/10.1021/acs.jcim.9b00362')
    parser.add_argument('-disjoint-allowed', '--disjoint-appearing-disappearing', dest='allow_disjoint_components',
                        type=str2bool, required=False, default=False,
                        help='Allow the molecules to be divided by "disappearing/appearing" atoms.')
    # temporary
    parser.add_argument('-amberff', '--amberff-name', metavar='AmberFFName', dest='amber_tleap_forcefield',
                        type=str, required=False, default='leaprc.protein.ff14SB',
                        help='This is a temporary solution. Ambertools tleap filename of the amberforcefield to be used. '
                             'Another one is "leaprc.ff99SBildn" etc. ')
    parser.add_argument('-ligff', '--ligandff-name', metavar='LigandFF', dest='ligand_ff',
                        type=str, required=False, default='gaff',
                        help='Either "gaff" or "gaff2"')
    parser.add_argument('-namd_prod', '--namd-prod', metavar='NAMD_prod', dest='namd_prod',
                        type=str, required=False, default='prod.namd',
                        help='This is a temporary solution. The name of the file to be used for the production. ')
    args = parser.parse_args()

    # ------------------ Configuration

    # set the working directory
    if args.cwd:
        # user provided
        workplace_root = args.cwd
    else:
        workplace_root = Path(os.getcwd())
    print(f'Working Directory: {workplace_root}')

    # check if ligand files are fine
    if args.ligands is None or len(args.ligands) < 2:
        print('Please supply at least two ligand files with -l (--ligands). E.g. -l file1.pdb file2.pdb ')
        sys.exit()

    if args.antechamber_dr is True:
        antechamber_dr = 'yes'
    else:
        antechamber_dr = 'no'
    print(f'Antechamber dr: {antechamber_dr}')

    # create ligands
    ligands = [Ligand(lig, workplace_root) for lig in args.ligands]

    command = args.action

    if command == 'rename':
        print('Atom names will be renamed to ensure that the atom names are unique across the two molecules.')
        # this case assumes that there are only two ligands given
        if len(ligands) != 2:
            print('ERROR: to rename atoms to be unique across two ligands, '
                  'you have to provide exactly two ligands with -l.'
                  'E.g. ties rename -l init.mol2 final.mol2')
            sys.exit()
        renameAtomNamesUniqueAndResnames(ligands[0], ligands[1])
        sys.exit()

    # fixme
    # If .ac format (ambertools, similar to .pdb), convert it to .mol2,
    ligands_right_format = [convert_ac_to_mol2(lig, 'left', args, workplace_root, antechamber_dr)
                            for lig in ligands]

    # ensure that each atom name is unique (see #237)
    # this helps to track further problems later
    [lig.make_atom_names_unique() for lig in ligands]

    if command != 'create':
        print('please provide action (rename/create)')
        sys.exit()

    if not args.protein:
        protein_filename = None
        print('No protein was select. Generating only one delta G directory.')
    else:
        protein_filename = Path(args.protein)
        if not protein_filename.is_file():
            print(f'Protein file (-p) name/path does not seem to lead to a file: {protein_filename}')
            sys.exit(1)

    align_molecules = args.align_mcs
    # use the original coords because antechamber can change them slightly
    use_original_coor = False

    ambertools_bin = find_antechamber(args)

    if args.net_charge is None:
        print('Please supply the net charge of the ligands with -nc. '
              'For now they all have to have the same charges.')
        sys.exit()
    else:
        net_charge = args.net_charge

    # charge tolerance
    atom_pair_q_atol = args.qtol
    net_charge_threshold = args.net_charge_threshold

    if args.ligands_have_q:
        # fixme - check if this is .mol2 or pdb
        use_provided_charges = args.ligands_have_q
        file_input_type = 'mol2'
        for lig in ligands:
            if lig.suffix != '.mol2':
                print('ERROR: if q are provided, the filetype .mol2/.ac has to be used across the molecules.')
                sys.exit(1)
    else:
        # determine whether charges are provided using the file extensions
        if all(lig.suffix() == '.mol2' for lig in ligands):
            print('Assuming that charges are provided based on the filetype .ac/.mol2')
            use_provided_charges = True
            file_input_type = 'mol2'
        elif all(lig.suffix() == '.pdb' for lig in ligands):
            print('Assuming that charges are not provided based on the filetype .ac/.mol2')
            use_provided_charges = False
            file_input_type = 'pdb'
        else:
            raise ValueError('The requested ligand files have different formats. This is not yet supported.')

    if use_provided_charges:
        # ignore the charges
        antechamber_charge_type = []
    else:
        # compute the charges
        antechamber_charge_type = ['-c', 'bcc'] # 'AM1-BCC'

    allow_disjoint_components = args.allow_disjoint_components
    print(f'Allowing disjoint components: {allow_disjoint_components}')

    # A user-provided list of pairs that should be matched
    # fixme - this requires better parsing capabilities
    # fixme - this should be handled by a function and throw ArgumentException if something is off
    manually_matched = []
    if args.manual_match_file is not None:
        with open(args.manual_match_file) as IN:
            for left_atom, right_atom in csv.reader(IN, delimiter='-'):
                manually_matched.append((left_atom.strip(), right_atom.strip()))
    if len(manually_matched) > 1:
        raise NotImplementedError('Currently only one atom pair can be matched - others were not tested')
    force_mismatch = []
    if args.manual_mismatch_file is not None:
        with open(args.manual_mismatch_file) as IN:
            for left_atom, right_atom in csv.reader(IN, delimiter='-'):
                force_mismatch.append((left_atom, right_atom))
        # fixme - check that there is no overlap between match and mismatch lists

    # ambertools forcefield
    amber_forcefield = args.amber_tleap_forcefield
    print(f'Amber forcefield name: {amber_forcefield}')

    namd_prod = "prod_2017.namd"  # only different because uses Berendsen
    print(f'The NAMD production file used: {namd_prod}')

    redist_q_over_unmatched = args.redistribute_charges_over_unmatched
    print(f'Distribute of introduced q disparity in the alchemical region: {redist_q_over_unmatched}')

    # use NAMD hybrid single dual topology
    use_hybrid_single_dual_top = args.hybrid_single_dual_top
    if use_hybrid_single_dual_top:
        ignore_charges_completely = True
    else:
        ignore_charges_completely = False

    # used for naming atom types,
    # fixme - we have to make sure this is consistent across the files (and ff leap.in files)
    atom_type = args.ligand_ff
    if atom_type == 'gaff':
        # they both use the same ff
        ligand_ff = 'leaprc.gaff'
    elif atom_type == 'gaff2':
        ligand_ff = 'leaprc.gaff2'
    else:
        raise ValueError('Argument -ligff cannot be anything else but "gaff" or "gaff2" currently')

    # not configurable currently
    hpc_submit = None  # "hpc_hartree_hsp.sh"

    # INTERNAL CONFIGURATION
    # set the path to the scripts
    code_root = Path(os.path.dirname(__file__))
    # scripts/input files
    script_dir = code_root / PurePosixPath('scripts')
    namd_script_dir = script_dir / 'namd'
    ambertools_script_dir = script_dir / 'ambertools'

    # Start of TIES

    # subprocess options for calling ambertools
    subprocess_kwargs = {
        "check" : True, "text" : True,
        "cwd" : workplace_root,
        "timeout" : 60 * 60 # 60 minute timeout
    }

    # prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
    [lig.antechamber_prepare_mol2(ambertools_bin, atom_type, net_charge, antechamber_dr, antechamber_charge_type)
        for lig in ligands]

    # when the user provides charges, BCC and minimisation is not carried out, so the coordinates are correct
    if use_original_coor and not use_provided_charges:
        print(f'Copying coordinates from {left_ligand} and {right_ligand} since antechamber changes them slightly')
        # copy the files before applying the coordinates
        shutil.copy(workplace_root / 'left.mol2', workplace_root / 'left_before_COOR.mol2')
        shutil.copy(workplace_root / 'right.mol2', workplace_root / 'right_before_COOR.mol2')
        set_coor_from_ref(workplace_root / 'left.mol2', workplace_root / left_ligand, by_atom_name=True)
        set_coor_from_ref(workplace_root / 'right.mol2', workplace_root / right_ligand, by_atom_name=True)
        # set_coor_from_ref(workplace_root / 'left.mol2', workplace_root / left_ligand, by_index=True)
        # set_coor_from_ref(workplace_root / 'right.mol2', workplace_root / right_ligand, by_index=True)

    # generate all possible morphs
    morphs = [Morph(ligA, ligZ, workplace_root) for ligA, ligZ in itertools.combinations(ligands, r=2)]

    # superimpose the two topologies
    for morph in morphs:
        # rename the atom names to ensure they are unique across the two molecules
        # we need to execute our pipeline for every pair and create the directories
        # which means we'll end up with a set of pairs
        morph.unique_atomres_names()

        suptop, mda_l1, mda_l2 = getSuptop(morph.renamed_ligA, morph.renamed_ligZ,
                                           align_molecules=align_molecules,
                                           pair_charge_atol=atom_pair_q_atol,
                                           manual_match=manually_matched,
                                           force_mismatch=force_mismatch,
                                           net_charge_threshold=net_charge_threshold,
                                           redistribute_charges_over_unmatched=redist_q_over_unmatched,
                                           ignore_charges_completely=ignore_charges_completely,
                                           no_disjoint_components=not allow_disjoint_components)

        morph.set_suptop(suptop, mda_l1, mda_l2)
        # save meta data
        morph.write_superimposition_json()
        morph.write_morph_pdb(hybrid_single_dual_top=use_hybrid_single_dual_top)
        morph.write_hybrid_mol2()

        morphs.append(morph)

    def rewrite_mol2_hybrid_top(file, single_top_atom_names):
        # in the  case of the hybrid single-dual topology in NAMD
        # the .mol2 files have to be rewritten so that
        # the atoms dual-topology atoms that appear/disappear
        # are placed at the beginning of the molecule
        # (single topology atoms have to be separted by
        # the same distance)
        shutil.copy(file, os.path.splitext(file)[0] + '_before_sdtop_reordering.mol2' )
        u = mda.Universe(file)
        # select the single top area, use their original order
        single_top_area = u.select_atoms('name ' +  ' '.join(single_top_atom_names))
        # all others are mutating
        dual_top_area = u.select_atoms('not name ' + ' '.join(single_top_atom_names))
        new_order_u = single_top_area + dual_top_area
        new_order_u.atoms.write(file)

    if use_hybrid_single_dual_top:
        raise NotImplementedError('Hybrid single-dual not done with LOMAP feature. ')
        rewrite_mol2_hybrid_top('left.mol2', list(matching_info["single_top_matched"].keys()))
        rewrite_mol2_hybrid_top('right.mol2', list(matching_info["single_top_matched"].values()))

    # generate the frcmod for each ligand
    print('Ambertools parmchk2 generating .frcmod for ligands')
    [lig.generate_frcmod(ambertools_bin / 'parmchk2', atom_type) for lig in ligands]

    # join the .frcmod files for each pair
    print('Ambertools parmchk2 generating .frcmod for morphs')
    [morph.join_frcmod_files(ambertools_bin, ambertools_script_dir, amber_forcefield, ligand_ff) for morph in morphs]

    ##########################################################
    # ------------------   Ligand ----------------------------
    # pick the right tleap instructions
    if use_hybrid_single_dual_top:
        ligand_tleap_in = 'leap_ligand_sdtop.in'
    else:
        ligand_tleap_in = 'leap_ligand.in'

    # fixme - at this point you'd know which pairs to set up
    for morph in morphs:
        prepare_inputs(morph,
                       dir_prefix='lig',
                       protein=None,
                       namd_script_loc=namd_script_dir,
                       scripts_loc=script_dir,
                       tleap_in=ligand_tleap_in,
                       protein_ff=amber_forcefield,
                       ligand_ff=ligand_ff,
                       net_charge=net_charge,
                       ambertools_script_dir=ambertools_script_dir,
                       subprocess_kwargs=subprocess_kwargs, # drop this?
                       ambertools_bin=ambertools_bin,
                       namd_prod=namd_prod,
                       hybrid_topology=use_hybrid_single_dual_top
                       )
        print(f'Ligand {morph} directory populated successfully')

    ##########################################################
    # ------------------ Complex  ----------------------------
    # calculate the charges of the protein (using ambertools)
    if protein_filename is not None:
        protein_net_charge = get_protein_net_charge(workplace_root, protein_filename.resolve(),
                               ambertools_bin, ambertools_script_dir / 'solv_prot.in',
                               amber_forcefield)
        print(f'Protein net charge: {protein_net_charge}')

        # pick the right tleap instuctions
        if use_hybrid_single_dual_top:
            complex_tleap_in = 'leap_complex_sdtop.in'
        else:
            complex_tleap_in = 'leap_complex.in'

        prepare_inputs(workplace_root, directory='complex',
                       protein=protein_filename,
                       hybrid_mol2=hybrid_mol2,
                       hybrid_frc=hybrid_frcmod,
                       left_right_mapping=matching_json,
                       namd_script_loc=namd_script_dir,
                       scripts_loc=script_dir,
                       tleap_in=complex_tleap_in,
                       protein_ff=amber_forcefield,
                       ligand_ff=ligand_ff,
                       net_charge=net_charge + protein_net_charge,
                       ambertools_script_dir=ambertools_script_dir,
                       subprocess_kwargs=subprocess_kwargs,
                       ambertools_bin=ambertools_bin,
                       namd_prod=namd_prod,
                       hybrid_topology=use_hybrid_single_dual_top)

    # prepare the post-analysis scripts
    shutil.copy(namd_script_dir / "check_namd_outputs.py", workplace_root)
    shutil.copy(namd_script_dir / "ddg.py", workplace_root)

    print('TIES 20 Finished')

if __name__ == '__main__':
    command_line_script()