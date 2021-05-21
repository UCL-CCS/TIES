import sys
import os
import subprocess
import shutil
from pathlib import Path

import numpy as np

import ties.helpers
from ties.config import Config


class Ligand:
    """
    The ligand helper class. It tracks the different copies of the original input files.
    It also offers ligand-oriented operations.

    TODO - use a general conf class/dict to find out about ambertools
    """
    LIG_COUNTER = 0

    ROOT = Path('prep')
    UNIQ_ATOM_NAME_DIR = ROOT / 'unique_atom_names'
    FRCMOD_DIR = ROOT / 'ligand_frcmods'
    ACPREP_CONVERT = ROOT / 'acprep_to_mol2'
    MOL2 = ROOT / 'mol2'

    _USED_FILENAMES = set()

    def __init__(self, ligand, config=None, save=True):
        self.save = save
        # save workplace root
        self.config = Config() if config is None else config
        self.original_input = Path(ligand).absolute()

        # check if the input files exist
        if not self.original_input.is_file():
            print(f'ERROR: Ligand file {self.original_input} not found')
            sys.exit(1)

        # internal name without an extension
        self.internal_name = self.original_input.stem

        # ligand names have to be unique
        if self.internal_name in Ligand._USED_FILENAMES:
            print(f'ERROR: the ligand filename {self.internal_name} is not unique in the list of ligands. ')
            sys.exit(1)
        else:
            Ligand._USED_FILENAMES.add(self.internal_name)

        # last used representative Path file
        self.current = self.original_input

        # internal index
        self.index = Ligand.LIG_COUNTER
        Ligand.LIG_COUNTER += 1

        self.ligand_with_uniq_atom_names = None

        # If .ac format (ambertools, similar to .pdb), convert it to .mol2 using antechamber
        self.convert_acprep_to_mol2()

    def __repr__(self):
        # return self.original_input.stem
        return self.internal_name

    def convert_acprep_to_mol2(self):
        """
        If the file is not a prep/ac file, this function does not do anything.
        Antechamber is called to convert the .prepi/.prep/.ac file into a .mol2 file.

        Returns: the name of the original file, or of it was .prepi, a new filename with .mol2
        """

        if self.current.suffix.lower() not in ('.ac', '.prep'):
            return

        filetype = {'.ac': 'ac', '.prep': 'prepi'}[self.current.suffix.lower()]

        cwd = self.config.workdir / Ligand.ACPREP_CONVERT / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        # prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
        print(f'Antechamber: converting {filetype} to mol2')
        new_current = cwd / (self.internal_name + '.mol2')

        log_filename = cwd / "antechamber_conversion.log"
        with open(log_filename, 'w') as LOG:
            try:
                subprocess.run([self.config.ambertools_antechamber,
                                '-i', self.current, '-fi', filetype,
                                '-o', new_current, '-fo', 'mol2',
                                '-dr', self.config.antechamber_dr],
                               stdout=LOG, stderr=LOG,
                               check=True, text=True,
                               cwd=cwd, timeout=30)
            except subprocess.CalledProcessError as E:
                print('ERROR: An error occurred during the antechamber conversion from .ac to .mol2 data type. ')
                print(f'ERROR: The output was saved in the directory: {cwd}')
                print(f'ERROR: Please see the log file for the exact error information: {log_filename}')
                raise E

        # update
        self.original_ac = self.current
        self.current = new_current
        print(f'Converted .ac file to .mol2. The location of the new file: {self.current}')

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

        # todo - if they already unique, ignore

        # save the output here
        os.makedirs(self.config.workdir / Ligand.UNIQ_ATOM_NAME_DIR, exist_ok=True)

        # load the ligand with MDAnalysis
        ligand_universe = ties.helpers.load_MDAnalysis_atom_group(self.current)

        # ensure that all the atom names are unique
        atom_names = [a.name for a in ligand_universe.atoms]
        names_unique = len(set(atom_names)) == len(atom_names)

        if not names_unique or not ties.helpers.are_correct_names(atom_names):
            print(f'Atom names in the molecule ({self.original_input}/{self.internal_name}) are either not unique '
                  f'or do not follow NameDigit format (e.g. C15). Renaming')
            _, old_new_atomname_map = ties.helpers.rename_ligand(ligand_universe.atoms)
        else:
            # the atom names stay the same, create a map fro the same elements to the same?
            old_new_atomname_map = {a.name: a.name for a in ligand_universe.atoms}

        self.old_new_atomname_map = old_new_atomname_map

        ligand_with_uniq_atom_names = self.config.workdir / Ligand.UNIQ_ATOM_NAME_DIR / \
                                      (self.internal_name + self.current.suffix)
        if self.save:
            ligand_universe.atoms.write(ligand_with_uniq_atom_names)

        self.ligand_with_uniq_atom_names = ligand_with_uniq_atom_names
        self.universe = ligand_universe
        # this object is now represented by the updated ligand
        self.current = ligand_with_uniq_atom_names

    # make this into a python property
    def suffix(self):
        return self.current.suffix.lower()

    def antechamber_prepare_mol2(self, atom_type, net_charge, antechamber_dr, antechamber_charge_type):
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

        mol2_cwd = self.config.workdir / self.MOL2 / self.internal_name

        # prepare the directory
        if not mol2_cwd.is_dir():
            mol2_cwd.mkdir(parents=True, exist_ok=True)

        mol2_target = mol2_cwd / f'{self.internal_name}.mol2'

        # copy the existing file if the file is already .mol2
        if self.current.suffix == '.mol2':
            print(f'Already .mol2 used. Copying {self.current} to {mol2_cwd}. ')
            shutil.copy(self.current, mol2_target)

        # do not redo if the target file exists
        if not (mol2_target).is_file():
            log_filename = mol2_cwd / "antechamber.log"
            with open(log_filename, 'w') as LOG:
                try:
                    subprocess.run([self.config.ambertools_antechamber,
                                    '-i', self.current, '-fi', self.current.suffix[1:],
                                    '-o', mol2_target, '-fo', 'mol2',
                                    '-at', atom_type, '-nc', str(net_charge),
                                    '-dr', antechamber_dr] + antechamber_charge_type,
                                   cwd=mol2_cwd,
                                   stdout=LOG, stderr=LOG,
                                   check=True, text=True,
                                   timeout=60 * 30  # 30 minutes
                                   )
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
        mol2_u = ties.helpers.load_MDAnalysis_atom_group(self.current)
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
        print(f'INFO: frcmod for {self} was computed before. Not repeating.')
        if hasattr(self, 'frcmod'):
            return

        # fixme - work on the file handles instaed of the constant stitching
        print(f'Parmchk2: generate the .frcmod for {self.internal_name}.mol2')

        # prepare cwd
        cwd = self.config.workdir / Ligand.FRCMOD_DIR / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        target_frcmod = f'{self.internal_name}.frcmod'
        log_filename = cwd / "parmchk2.log"
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
                print('ERROR file content: ', open(log_filename).read())
                print('ERROR: An error occured during the antechamber conversion from .ac to .mol2 data type. ')
                print(f'ERROR: The output was saved in the directory: {cwd}')
                print(f'ERROR: Please see the log file for the exact error information: {log_filename}')
                raise E

        print(f'Parmchk2: created .frcmod: {target_frcmod}')
        self.frcmod = cwd / target_frcmod

    def overwrite_coordinates_with(self, file, output_file):
        """
        Load coordinates from another file and overwrite the coordinates in the current file.
        """

        # load the current atoms with MDAnalysis
        mda_template = ties.helpers.load_MDAnalysis_atom_group(self.current)

        # load the file with the coordinates we want to use
        coords = ties.helpers.load_MDAnalysis_atom_group(file)

        # fixme: use the atom names
        by_atom_name = True
        by_index = False
        by_general_atom_type = False

        # mol2_filename will be overwritten!
        print(f'Writing to {self.current} the coordinates from {file}. ')

        coords_sum = np.sum(coords.atoms.positions)

        if by_atom_name and by_index:
            raise ValueError('Cannot have both. They are exclusive')
        elif not by_atom_name and not by_index:
            raise ValueError('Either option has to be selected.')

        if by_general_atom_type:
            for mol2_atom in mda_template.atoms:
                found_match = False
                for ref_atom in coords.atoms:
                    if element_from_type[mol2_atom.type.upper()] == element_from_type[ref_atom.type.upper()]:
                        found_match = True
                        mol2_atom.position = ref_atom.position
                        break
                assert found_match, "Could not find the following atom in the original file: " + mol2_atom.name
        if by_atom_name:
            for mol2_atom in mda_template.atoms:
                found_match = False
                for ref_atom in coords.atoms:
                    if mol2_atom.name.upper() == ref_atom.name.upper():
                        found_match = True
                        mol2_atom.position = ref_atom.position
                        break
                assert found_match, "Could not find the following atom name across the two files: " + mol2_atom.name
        elif by_index:
            for mol2_atom, ref_atom in zip(mda_template.atoms, coords.atoms):
                atype = element_from_type[mol2_atom.type.upper()]
                reftype = element_from_type[ref_atom.type.upper()]
                if atype != reftype:
                    raise Exception(
                        f"The found general type {atype} does not equal to the reference type {reftype} ")

                mol2_atom.position = ref_atom.position

        if np.testing.assert_almost_equal(coords_sum, np.sum(mda_template.atoms.positions), decimal=2):
            print('Different positions sums:', coords_sum, np.sum(mda_template.atoms.positions))
            raise Exception('Copying of the coordinates did not work correctly')

        # save the output file
        mda_template.atoms.write(output_file)

        print('hi')
        pass
