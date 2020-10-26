#!/usr/bin/env python3
"""
Load two ligands, run the topology superimposer, and then using the results, generate the NAMD input files.
"""
import os
import shutil
import sys
import subprocess
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


def prepare_antechamber_mol2_input(molecule, workplace_root, ambertools_bin, left_ligand,
                                   file_input_type, atom_type, net_charge,
                                   antechamber_dr, antechamber_charge_type):
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

    storage_directory = f"initial_{molecule}_preparation"

    # prepare the directory
    if not os.path.isdir(storage_directory):
        os.makedirs(storage_directory)

    cwd = workplace_root / storage_directory

    subprocess_kwargs = {
        "check" : True, "text" : True,
        "timeout" : 60 * 30 # 30 minutes
    }

    # fixme - does not throw errors? try to make it throw an error
    # todo - save the output in a file
    # todo - ensure that the antechamber is run in the right directory, which is a different directory
    log_filename = cwd / "antechamber.log"
    with open(log_filename, 'w') as LOG:
        try:
            subprocess.run([ambertools_bin / 'antechamber',
                            '-i', left_ligand, '-fi', file_input_type,
                            '-o', f'{molecule}.mol2', '-fo', 'mol2',
                            '-at', atom_type, '-nc', str(net_charge),
                            '-dr', antechamber_dr] + antechamber_charge_type,
                           cwd=cwd,
                           stdout=LOG, stderr=LOG,
                           **subprocess_kwargs)
        except subprocess.CalledProcessError as E:
            print('ERROR: occured when creating the input .mol2 file with antechamber. ')
            print(f'ERROR: The output was saved in the directory: {cwd}')
            print(f'ERROR: can be found in the file: {log_filename}')
            raise E
    print(f'Converted {molecule} into {molecule}.mol2, Log: {log_filename}')
    # copy the output to the root directory, e.g.  root/left.mol2
    shutil.copy(cwd / f'{molecule}.mol2', workplace_root / f'{molecule}.mol2')


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


def generate_frcmod(parmchk2, workplace_root, molecule, atom_type):
    print(f'Parmchk2: generate the .frcmod for {molecule}.mol2')

    # prepare the directory
    cwd_dir = workplace_root / f'generating_frcmod_{molecule}'
    if not os.path.isdir(cwd_dir):
        os.makedirs(cwd_dir)

    subprocess_kwargs = {
        "check": True, "text": True,
        "cwd": cwd_dir,
        "timeout": 60 * 5  # 5 minutes
    }

    log_filename = cwd_dir / "parmchk2.log"
    with open(log_filename, 'w') as LOG:
        try:
            subprocess.run([parmchk2,
                                '-i', (workplace_root / f'{molecule}.mol2').resolve(),
                                '-o', f'{molecule}.frcmod',
                                '-f', 'mol2',
                                '-s', atom_type],
                               stdout=LOG, stderr=LOG,
                               **subprocess_kwargs)
        except subprocess.CalledProcessError as E:
            print('ERROR: An error occured during the antechamber conversion from .ac to .mol2 data type. ')
            print(f'ERROR: The output was saved in the directory: {cwd_dir}')
            print(f'ERROR: Please see the log file for the exact error information: {log_filename}')
            raise E

    frcmod = cwd_dir / (str(molecule) + ".frcmod")
    print(f'Parmchk2: created successfully the .frcmod. Location: {frcmod}')
    return frcmod


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


def command_line_script():
    parser = argparse.ArgumentParser(description='TIES 20')
    parser.add_argument('action', metavar='command', type=str,
                        help='Action to be performed. E.g. "ties rename, ties create, ties .." ')
    parser.add_argument('-l', '--left-ligand', metavar='Left_Ligand_File', dest='left_ligand',
                        type=Path, required=False,
                        help='The left ligand file')
    parser.add_argument('-r', '--right-ligand', metavar='Right_Ligand_File', dest='right_ligand',
                        type=Path, required=False,
                        help='The right ligand file')
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
    if args.left_ligand is None or args.right_ligand is None:
        print('Please supply ligand files with -l and -r')
        sys.exit()

    left_ligand = args.left_ligand
    right_ligand = args.right_ligand

    # check if the input files exist
    # todo - verify the files at this stage
    if not (workplace_root / left_ligand).is_file():
        print(f'File {left_ligand} not found in {workplace_root}')
        sys.exit(1)
    elif not (workplace_root / right_ligand).is_file():
        print(f'File {right_ligand} not found in {workplace_root}')
        sys.exit(1)

    if args.antechamber_dr is True:
        antechamber_dr = 'yes'
    else:
        antechamber_dr = 'no'
    print(f'Will use antechamber dr: {antechamber_dr}')

    # fixme
    # if it is .ac file format (ambertools, almost the same as .pdb), convert it to .mol2,
    left_ligand = convert_ac_to_mol2(left_ligand, 'left', args, workplace_root, antechamber_dr)
    right_ligand = convert_ac_to_mol2(right_ligand, 'right', args, workplace_root, antechamber_dr)

    if args.action == 'rename':
        print('Atom names will be renamed to ensure that the atom names are unique across the two molecules.')
        renameAtomNamesUniqueAndResnames(left_ligand, right_ligand)
        sys.exit()
    elif args.action == 'create':
        print('Main protocol is used. ')
    else:
        print('please provide action (rename/create)')
        sys.exit()

    if not args.protein:
        protein_filename = None
        print('No protein was select. Generating only one delta G directory.')
    else:
        protein_filename = args.protein
        if not os.path.isfile(protein_filename):
            print(f'Protein file (-p) name/path does not seem to lead to a file: {protein_filename}')
            sys.exit(1)

    align_molecules = args.align_mcs
    # use the original coords because antechamber can change them slightly
    use_original_coor = True

    ambertools_bin = find_antechamber(args)

    if args.net_charge is None:
        print('Please supply the net charge of the ligands with -nc')
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
    else:
        # determine automatically based on the file extension
        left_ext = os.path.splitext(left_ligand)[-1].lower()
        right_ext = os.path.splitext(right_ligand)[-1].lower()
        if left_ext == '.mol2' and right_ext == '.mol2':
            use_provided_charges = True
            file_input_type = 'mol2'
        elif left_ext == '.pdb' and right_ext == '.pdb':
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
    print(f'Configuration: allowing disjoint components {allow_disjoint_components}')

    # A list of pairs that should be matched
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
    print(f'Will use amber forcefield name {amber_forcefield}')

    namd_prod = "prod_2017.namd"  # only different because uses Berendsen
    print(f'The NAMD production file in use is {namd_prod}')

    redist_q_over_unmatched = args.redistribute_charges_over_unmatched
    print(f'Will distribute introduced q disparity in the unmatched region. {redist_q_over_unmatched}')

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

    ## Start of TIES

    # subprocess options for calling ambertools
    subprocess_kwargs = {
        "check" : True, "text" : True,
        "cwd" : workplace_root,
        "stdout" : sys.stdout, #subprocess.PIPE,
        "stderr" : sys.stderr, #subprocess.PIPE,
        "timeout" : 60 * 60 # 60 minute timeout
    }

    # prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
    if not (workplace_root / 'left.mol2').is_file():
        prepare_antechamber_mol2_input("left", workplace_root, ambertools_bin, left_ligand,
                                       file_input_type, atom_type,
                                       net_charge, antechamber_dr,
                                       antechamber_charge_type)
    if not (workplace_root / 'right.mol2').is_file():
        prepare_antechamber_mol2_input("right", workplace_root, ambertools_bin, right_ligand,
                                       file_input_type, atom_type,
                                       net_charge, antechamber_dr,
                                       antechamber_charge_type)

    # scan the files to see if it created DUMMY DU atoms in the .mol2, remove the atoms if that's the case
    removeDU_atoms(workplace_root / 'left.mol2')
    removeDU_atoms(workplace_root / 'right.mol2')

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

    # rename the atom names to ensure they are unique across the two molecules
    renameAtomNamesUniqueAndResnames(workplace_root / 'left.mol2', workplace_root / 'right.mol2')

    # superimpose the two topologies
    suptop, mda_l1, mda_l2 = getSuptop(workplace_root / 'left.mol2', workplace_root / 'right.mol2',
                                       align_molecules=align_molecules,
                                       pair_charge_atol=atom_pair_q_atol,
                                       manual_match=manually_matched, force_mismatch=force_mismatch,
                                       net_charge_threshold=net_charge_threshold,
                                       redistribute_charges_over_unmatched=redist_q_over_unmatched,
                                       ignore_charges_completely=ignore_charges_completely,
                                       no_disjoint_components=not allow_disjoint_components)

    # save the superimposition results
    left_right_matching_json = workplace_root / 'joint_meta_fep.json'
    left_right_matching = save_superimposition_results(left_right_matching_json, suptop)
    # hybrid .pdb
    write_morph_top_pdb(workplace_root / 'morph.pdb', mda_l1, mda_l2, suptop,
                        hybrid_single_dual_top=use_hybrid_single_dual_top)
    # hybrid .mol2
    hybrid_mol2 = workplace_root / 'morph.mol2'
    # we save the merged .mol2 regardless of whether a dual or hybrid is used
    write_merged_mol2(suptop, hybrid_mol2)

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
        rewrite_mol2_hybrid_top('left.mol2', list(left_right_matching["single_top_matched"].keys()))
        rewrite_mol2_hybrid_top('right.mol2', list(left_right_matching["single_top_matched"].values()))


    # generate the functional forms
    print('Ambertools parmchk2 generating .frcmod')
    left_frcmod = generate_frcmod(ambertools_bin / 'parmchk2', workplace_root, 'left', atom_type)
    right_frcmod = generate_frcmod(ambertools_bin / 'parmchk2', workplace_root, 'right', atom_type)

    # join the .frcmod files
    hybrid_frcmod = workplace_root / 'morph.frcmod'
    join_frcmod_files2(left_frcmod, right_frcmod, hybrid_frcmod)

    # if the hybrid .frcmod needs new terms between the appearing/disappearing atoms, insert dummy ones
    # fixme - this should be a part of the "joining above"
    updated_frcmod_content = check_hybrid_frcmod(hybrid_mol2, hybrid_frcmod, amber_forcefield, ligand_ff,
                                                 ambertools_bin, ambertools_script_dir, cwd=workplace_root)
    with open(hybrid_frcmod, 'w') as FOUT:
        FOUT.write(updated_frcmod_content)

    ##########################################################
    # ------------------   Ligand ----------------------------
    # pick the right tleap instuctions
    if use_hybrid_single_dual_top:
        ligand_tleap_in = 'leap_ligand_sdtop.in'
    else:
        ligand_tleap_in = 'leap_ligand.in'

    prepare_inputs(workplace_root, directory='lig',
                   protein=None,
                   hybrid_mol2=hybrid_mol2,
                   hybrid_frc=hybrid_frcmod,
                   left_right_mapping=left_right_matching_json,
                   namd_script_loc=namd_script_dir,
                   scripts_loc=script_dir,
                   tleap_in=ligand_tleap_in,
                   protein_ff=amber_forcefield,
                   ligand_ff=ligand_ff,
                   net_charge=net_charge,
                   ambertools_script_dir=ambertools_script_dir,
                   subprocess_kwargs=subprocess_kwargs,
                   ambertools_bin=ambertools_bin,
                   namd_prod=namd_prod,
                   hybrid_topology=use_hybrid_single_dual_top
                   )

    ##########################################################
    # ------------------ Complex  ----------------------------
    # calculate the charges of the protein (using ambertools)
    if protein_filename is not None:
        protein_net_charge = get_protein_net_charge(workplace_root, protein_filename,
                               ambertools_bin, ambertools_script_dir / 'solv_prot.in',
                               subprocess_kwargs, amber_forcefield)

        # pick the right tleap instuctions
        if use_hybrid_single_dual_top:
            complex_tleap_in = 'leap_complex_sdtop.in'
        else:
            complex_tleap_in = 'leap_complex.in'

        prepare_inputs(workplace_root, directory='complex',
                       protein=protein_filename,
                       hybrid_mol2=hybrid_mol2,
                       hybrid_frc=hybrid_frcmod,
                       left_right_mapping=left_right_matching_json,
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