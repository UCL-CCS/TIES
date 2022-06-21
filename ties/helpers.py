"""
A list of functions with a clear purpose that does
not belong specifically to any of the existing units.
"""
import os
import sys
import argparse
import subprocess
import warnings
import pathlib



def get_new_atom_names(atoms, name_counter=None):
    """
    todo - add unit tests

    @parameter/returns name_counter: a dictionary with atom as the key such as 'N', 'C', etc,
    the counter keeps track of the last used counter for each name.
    Empty means that the counting will start from 1.
    input atoms: mdanalysis atoms
    """
    if name_counter is None:
        name_counter = {}

    renaming_map = {}

    for atom in atoms:
        # get the first letters that is not a character
        afterLetters = [i for i, l in enumerate(atom.name) if l.isalpha()][-1] + 1

        atom_name = atom.name[:afterLetters]
        last_used_counter = name_counter.get(atom_name, 0)

        # rename
        last_used_counter += 1
        newAtomName = atom_name + str(last_used_counter)
        renaming_map[newAtomName] = atom.name

        atom.name = newAtomName

        # update the counter
        name_counter[atom_name] = last_used_counter

    return name_counter, renaming_map


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


class ArgparseChecker():

    @staticmethod
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

    @staticmethod
    def existing_file(v):
        # check ligand arguments
        if not pathlib.Path(v).is_file():
            raise argparse.ArgumentTypeError(f'The file {v} could not be found.')
        return pathlib.Path(v)

    @staticmethod
    def ambertools_home(path):
        # check if this path points to ambertools
        amber_home = pathlib.Path(path)
        if not amber_home.is_dir():
            raise argparse.ArgumentTypeError(f'The path to ambertools does not point towards a directory. ')
        # check if the bin directory, antechamber and antechamber
        if not (amber_home / 'bin').is_dir():
            raise argparse.ArgumentTypeError(f'The path to ambertools does not contain the "bin" directory. ')
        if not (amber_home / 'bin' / 'antechamber').is_dir():
            raise argparse.ArgumentTypeError(f'The path to ambertools does not contain the "bin/antechamber" file. ')

        return amber_home