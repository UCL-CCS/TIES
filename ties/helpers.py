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

    # {new_uniqe_name: old_atom_name}
    reverse_renaming_map = {}

    for atom in atoms:
        # count the letters before any digit
        letter_count = 0
        for letter in atom.name:
            if not letter.isalpha():
                break

            letter_count += 1

        # use max 3 letters from the atom name
        letter_count = min(letter_count, 3)

        letters = atom.name[:letter_count]

        # how many atoms do we have with these letters? ie C1, C2, C3 -> 3
        last_used_counter = name_counter.get(letters, 0) + 1

        # rename
        new_name = letters + str(last_used_counter)

        # if the name is longer than 4 character,
        # shorten the number of letters
        if len(new_name) > 4:
            # the name is too long, use only the first character
            new_name = letters[:4-len(str(last_used_counter))] + str(last_used_counter)

            # we assume that there is fewer than 1000 atoms with that name
            assert len(str(last_used_counter)) < 1000

        reverse_renaming_map[new_name] = atom.name

        atom.name = new_name

        # update the counter
        name_counter[letters] = last_used_counter

    return name_counter, reverse_renaming_map


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

        # we are starting the counter from 0 as we always add 1 later on
        last_used_counter = name_counter.get(atom_name, 0)

        # update the counter
        name_counter[atom_name] = max(last_used_counter + 1, atom_number)

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
    def ratio(v):
        if ':' not in v:
            argparse.ArgumentTypeError('The ratio has to be separated by ":".')
        mcs, rmsd = v.split(':')
        return [float(mcs), float(rmsd)]

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