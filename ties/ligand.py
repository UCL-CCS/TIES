import sys
import os
import subprocess
import shutil
from pathlib import Path

import numpy as np
import parmed

import ties.helpers
from ties.config import Config


class Ligand:
    """
    The ligand helper class. Helps to load and manage the different copies of the ligand file.
    Specifically, it tracks the different copies of the original input files as it is transformed (e.g. charge assignment).

    :param ligand: ligand filepath
    :type ligand: string
    :param config: Optional configuration from which the relevant ligand settings can be used
    :type config: :class:`Config`
    :param save: write a file with unique atom names for further inspection
    :type save: bool
    """
    LIG_COUNTER = 0

    _USED_FILENAMES = set()

    def __init__(self, ligand, config=None, save=True):
        """Constructor method
        """

        self.save = save
        # save workplace root
        self.config = Config() if config is None else config
        self.config.ligand_files = ligand

        self.original_input = Path(ligand).absolute()

        # internal name without an extension
        self.internal_name = self.original_input.stem

        # ligand names have to be unique
        if self.internal_name in Ligand._USED_FILENAMES and self.config.uses_cmd:
            print(f'ERROR: the ligand filename {self.internal_name} is not unique in the list of ligands. ')
            sys.exit(1)
        else:
            Ligand._USED_FILENAMES.add(self.internal_name)

        # last used representative Path file
        self.current = self.original_input

        # internal index
        # TODO - move to config
        self.index = Ligand.LIG_COUNTER
        Ligand.LIG_COUNTER += 1

        self._renaming_map = None
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

        cwd = self.config.lig_acprep_dir / self.internal_name
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

    def are_atom_names_correct(self):
        """
        Checks if atom names:
         - are unique
         - have a correct format "LettersNumbers" e.g. C17
        """
        ligand = parmed.load_file(str(self.current), structure=True)
        atom_names = [a.name for a in ligand.atoms]

        are_uniqe = len(set(atom_names)) == len(atom_names)

        return are_uniqe and self._do_atom_names_have_correct_format(atom_names)

    @staticmethod
    def _do_atom_names_have_correct_format(names):
        """
        Check if the atom name is followed by a number, e.g. "C15"
        Note that the full atom name cannot be more than 4 characters.
        This is because the PDB format does not allow for more
        characters which can lead to inconsistencies.

        :param names: a list of atom names
        :type names: list[str]
        :return True if they all follow the correct format.
        """
        for name in names:
            # cannot exceed 4 characters
            if len(name) > 4:
                return False

            # count letters before any digit
            letter_count = 0
            for letter in name:
                if not letter.isalpha():
                    break

                letter_count += 1

            # at least one character
            if letter_count == 0:
                return False

            # extrac the number suffix
            atom_number = name[letter_count:]
            try:
                int(atom_number)
            except:
                return False

        return True

    def correct_atom_names(self):
        """
        Ensure that each atom name:
         - is unique
         - has letter followed by digits
         - has max 4 characters
        E.g. C17, NX23

        :param self.save: if the path is provided, the updated file
            will be saved with the unique names and a handle to the new file (ParmEd) will be returned.
        """
        if self.are_atom_names_correct():
            return

        print(f'Ligand {self.internal_name} will have its atom names renamed. ')

        ligand = parmed.load_file(str(self.current), structure=True)

        print(f'Atom names in the molecule ({self.original_input}/{self.internal_name}) are either not unique '
              f'or do not follow NameDigit format (e.g. C15). Renaming')
        _, renaming_map = ties.helpers.get_new_atom_names(ligand.atoms)
        self._renaming_map = renaming_map
        print(f'Rename map: {renaming_map}')

        # save the output here
        os.makedirs(self.config.lig_unique_atom_names_dir, exist_ok=True)

        ligand_with_uniq_atom_names = self.config.lig_unique_atom_names_dir / (self.internal_name + self.current.suffix)
        if self.save:
            ligand.save(str(ligand_with_uniq_atom_names))

        self.ligand_with_uniq_atom_names = ligand_with_uniq_atom_names
        self.parmed = ligand
        # this object is now represented by the updated ligand
        self.current = ligand_with_uniq_atom_names

    @property
    def renaming_map(self):
        """
        Otherwise, key: newName, value: oldName.

        If None, means no renaming took place.
        """
        return self._renaming_map

    @property
    def rev_renaming_map(self):
        return {b: c for c, b in self._renaming_map.items()}

    @renaming_map.setter
    def renaming_map(self, dict):
        if self._renaming_map is None:
            self._renaming_map = dict
        else:
            # this ligand was already renamed before.
            # so B -> A, where A is the original value,
            # and B is the new value (guarantted to be unique, therefore the key)
            # Now B is renamed to C, so we need to have C -> A
            # For each renamed B value here, we have to find the A value
            # fixme: this works only if Bs are unique

            # dict: C -> B. We know that C and B are unique. Therefore, reverse for convenience
            rev_dict = {b: c for c, b in dict.items()}
            self._renaming_map = {rev_dict[b]: a for b, a in self._renaming_map.items()}

    # make this into a python property
    def suffix(self):
        return self.current.suffix.lower()

    def antechamber_prepare_mol2(self, **kwargs):
        """
        Converts the ligand into a .mol2 format.

        BCC charges are generated if missing or requested.
        It calls antechamber (the charge type -c is not used if user prefers to use their charges).
        Any DU atoms created in the antechamber call are removed.

        :param atom_type: Atom type bla bla
        :type atom_type:
        :param net_charge:
        :type net_charge: int
        """
        self.config.set_configs(**kwargs)

        print('Antechamber: converting to .mol2 and generating charges if necessary')
        if self.config.ligands_contain_q:
            print('Antechamber: Generating .mol2 file with BCC charges')
        if not self.config.antechamber_charge_type:
            print('Antechamber: Ignoring atom charges. The user-provided atom charges will be used. ')
        
        mol2_cwd = self.config.lig_dir / self.internal_name

        # prepare the directory
        mol2_cwd.mkdir(parents=True, exist_ok=True)
        mol2_target = mol2_cwd / f'{self.internal_name}.mol2'

        # copy the existing file if the file is already .mol2
        if self.current.suffix == '.mol2' and self.config.ligands_contain_q:
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
                        '-at', self.config.ligand_ff_name, '-nc', str(self.config.ligand_net_charge),
                        '-dr', str(self.config.antechamber_dr)] + self.config.antechamber_charge_type,
                       cwd=mol2_cwd,
                       stdout=LOG, stderr=LOG,
                       check=True, text=True,
                       timeout=60 * 30  # 30 minutes
                       )
                except subprocess.CalledProcessError as ProcessError:
                    raise Exception(f'Could not convert the input into .mol2 file with antechamber. '
                                    f'See the log and its directory: {log_filename}') from ProcessError
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
        mol2 = parmed.load_file(str(self.current), structure=True)
        # check if there are any DU atoms
        has_DU = any(a.type == 'DU' for a in mol2.atoms)
        if not has_DU:
            return

        # make a backup copy before (to simplify naming)
        shutil.move(self.current, self.current.parent / ('lig.beforeRemovingDU' + self.current.suffix))

        # remove DU type atoms and save the file
        for atom in mol2.atoms:
            if atom.name != 'DU':
                continue

            atom.residue.delete_atom(atom)
        # save the updated molecule
        mol2.save(str(self.current))
        print('Removed dummy atoms with type "DU"')

    def generate_frcmod(self, **kwargs):
        """
            params
             - parmchk2
             - atom_type
        """
        self.config.set_configs(**kwargs)

        print(f'INFO: frcmod for {self} was computed before. Not repeating.')
        if hasattr(self, 'frcmod'):
            return

        # fixme - work on the file handles instaed of the constant stitching
        print(f'Parmchk2: generate the .frcmod for {self.internal_name}.mol2')

        # prepare cwd
        cwd = self.config.lig_frcmod_dir / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        target_frcmod = f'{self.internal_name}.frcmod'
        log_filename = cwd / "parmchk2.log"
        with open(log_filename, 'w') as LOG:
            try:
                subprocess.run([self.config.ambertools_parmchk2,
                                '-i', self.current,
                                '-o', target_frcmod,
                                '-f', 'mol2',
                                '-s', self.config.ligand_ff_name],
                               stdout=LOG, stderr=LOG,
                               check= True, text=True,
                               cwd= cwd, timeout=20,  # 20 seconds
                                )
            except subprocess.CalledProcessError as E:
                raise Exception(f"GAFF Error: Could not generate FRCMOD for file: {self.current} . "
                                f'See more here: {log_filename}') from E

        print(f'Parmchk2: created .frcmod: {target_frcmod}')
        self.frcmod = cwd / target_frcmod

    def overwrite_coordinates_with(self, file, output_file):
        """
        Load coordinates from another file and overwrite the coordinates in the current file.
        """

        # load the current atoms with ParmEd
        template = parmed.load_file(str(self.current), structure=True)

        # load the file with the coordinates we want to use
        coords = parmed.load_file(str(file), structure=True)

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
            for mol2_atom in template.atoms:
                found_match = False
                for ref_atom in coords.atoms:
                    if element_from_type[mol2_atom.type.upper()] == element_from_type[ref_atom.type.upper()]:
                        found_match = True
                        mol2_atom.position = ref_atom.position
                        break
                assert found_match, "Could not find the following atom in the original file: " + mol2_atom.name
        if by_atom_name:
            for mol2_atom in template.atoms:
                found_match = False
                for ref_atom in coords.atoms:
                    if mol2_atom.name.upper() == ref_atom.name.upper():
                        found_match = True
                        mol2_atom.position = ref_atom.position
                        break
                assert found_match, "Could not find the following atom name across the two files: " + mol2_atom.name
        elif by_index:
            for mol2_atom, ref_atom in zip(template.atoms, coords.atoms):
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