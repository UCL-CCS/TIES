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

import MDAnalysis


def are_correct_names(names):
    """
    Check if the atom name is followed by a number, e.g. "C15"
    @parameter names: a list of atom names
    @returns True if they all follow the correct format.
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


def load_mol2_wrapper(filename):
    # Load a .mol2 file
    # fixme - fuse with the .pdb loading, no distinction needed
    # ignore the .mol2 warnings about the mass
    warnings.filterwarnings(action='ignore', category=UserWarning,
                            message='Failed to guess the mass for the following atom types: '  # warning to ignore
                            )
    # squash the internal warning about parsing .mol2 within MDAnalysis
    warnings.filterwarnings(action='ignore', category=UserWarning,
                            message='Creating an ndarray from ragged nested sequences '
                            )
    # squash the internal warning about missing "cell dimensions"
    warnings.filterwarnings(action='ignore', category=UserWarning,
                            message='Unit cell dimensions not found. CRYST1 record set to unitary values.'
                            )
    u = MDAnalysis.Universe(filename)
    # turn off the filter warning after?
    return u


def rename_ligand(atoms, name_counter=None):
    """
    todo - add unit tests

    @parameter/returns name_counter: a dictionary with atom as the key such as 'N', 'C', etc,
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


def find_antechamber(args):
    # set up ambertools
    if args.ambertools_home is not None:
        # user manually specified the variable.
        ambertools_bin = args.ambertools_home / 'bin'
    elif os.getenv('AMBERHOME'):
        ambertools_bin = pathlib.Path(os.getenv('AMBERHOME')) / 'bin'
    elif os.getenv('AMBER_PREFIX'):
        ambertools_bin = pathlib.Path(os.getenv('AMBER_PREFIX')) / 'bin'
    else:
        print('Error: Cannot find ambertools. $AMBERHOME and $AMBER_PREFIX are empty')
        print('Option 1: source your ambertools script amber.sh')
        print('Option 2: specify manually the path to amberhome with -ambertools option')
        sys.exit()

    # fixme - test ambertools at this stage before proceeding
    return ambertools_bin