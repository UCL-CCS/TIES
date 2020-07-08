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

def command_line_script():
    parser = argparse.ArgumentParser(description='TIES 20')
    parser.add_argument('action', metavar='command', type=str,
                        help='Action to be performed. E.g. "ties rename .." ')
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
                        type=bool, required=False,
                        help='Use charges provided with the ligand files (with .mol2). '
                             'Default: If .mol2 is given, using the given charges will be attempted. '
                             'Default: If .pdb is given, then BCC charges are computed with antechamber. ')
    parser.add_argument('-align-mcs', '--align-ligands-mcs', metavar='Align_Ligands_MCS', dest='align_mcs',
                        type=bool, required=False, default=False,
                        help='Align the right ligand to the left using the determined common component')
    parser.add_argument('-ant-dr', '--antechamber-dr', metavar='AntechamberDrMode', dest='antechamber_dr',
                        type=bool, required=False, default=True,
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
    # temporary
    parser.add_argument('-amberff', '--amberff-name', metavar='AmberFFName', dest='amber_tleap_forcefield',
                        type=str, required=False, default='leaprc.protein.ff14SB',
                        help='This is a temporary solution. Ambertools tleap filename of the amberforcefield to be used. '
                             'Another one is "leaprc.ff99SBildn" etc. ')
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

    if args.action == 'rename':
        print('Atom names will be renamed to ensure that the atom names are unique across the two molecules.')
        renameAtomNamesUnique(left_ligand, right_ligand)
        sys.exit()
    elif args.action == 'create':
        print('Main protocol is used. ')
    else:
        print('please provide action (rename/create)')
        sys.exit()

    if not args.protein:
        print('Please supply the protein file with -p')
        sys.exit()
    else:
        protein_filename = args.protein
        if not os.path.isfile(protein_filename):
            print(f'Protein file (-p) name/path does not seem to lead to a file: {protein_filename}')
            sys.exit(1)

    align_molecules = args.align_mcs
    # use the original coords because antechamber can change them slightly
    use_original_coor = True

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
        print('Option 1: source your ambertools script.')
        print('Option 2: specify manually the path to amberhome with -ambertools')
        sys.exit()
    # fixme - test ambertools at this stage before proceeding

    if not args.net_charge:
        print('Please supply the net charge of the ligands with -nc')
        sys.exit()
    else:
        net_charge = args.net_charge

    # charge tolerance
    atom_pair_q_atol = args.qtol
    net_charge_threshold = args.net_charge_threshold

    if args.ligands_have_q is None:
        # determine automatically based on the file extension
        left_ext = os.path.splitext(left_ligand)[-1].lower()
        right_ext = os.path.splitext(right_ligand)[-1].lower()
        if left_ext == '.mol2' and right_ext == '.mol2':
            use_provided_charges = True
        elif left_ext == '.pdb' and right_ext == '.pdb':
            use_provided_charges = False
        else:
            raise ValueError('The requested ligand files have different formats. This is not yet supported.')
    else:
        # fixme - check if this is .mol2 or pdb
        use_provided_charges = args.ligands_have_q

    if use_provided_charges:
        # fixme - you cannot delete charges if that's what we want to keep
        # ignore the charges,
        charge_type = 'wc' # 'dc'  # 'Delete Charge'
        raise NotImplementedError('Does not deal with supplied charges at this moment')
    else:
        # compute the charges
        # fixme - test with .pdb and with .mol2
        charge_type = 'bcc'  # 'AM1-BCC'

    # A list of pairs that should be matched
    # fixme - this requires better parsing capabilities
    manually_matched = []
    if args.manual_match_file is not None:
        with open(args.manual_match_file) as IN:
            for left_atom, right_atom in csv.reader(IN, delimiter='-'):
                manually_matched.append((left_atom.strip(), right_atom.strip()))
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

    if args.antechamber_dr:
        antechamber_dr = 'on'
    else:
        antechamber_dr = 'off'
    print(f'Will use antechamber dr: {antechamber_dr}')

    # used for naming atom types,
    # fixme - we have to make sure this is consistent across the files (and ff leap.in files)
    atom_type = 'gaff'

    # not configurable currently
    hpc_submit = None  # "hpc_hartree_hsp.sh"

    # INTERNAL CONFIGURATION
    # set the path to the scripts
    code_root = Path(os.path.dirname(__file__))
    # scripts/input files
    script_dir = code_root / PurePosixPath('scripts')
    namd_script_dir = script_dir / 'namd'
    ambertools_script_dir = script_dir / 'ambertools'



    # subprocess options for calling ambertools
    subprocess_kwargs = {
        "check" : True, "text" : True,
        "cwd" : workplace_root,
        "stdout" : sys.stdout, #subprocess.PIPE,
        "stderr" : sys.stderr, #subprocess.PIPE,
        "timeout" : 60 * 10 # 10 minute timeout
    }

    # prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
    if not (workplace_root / 'left.mol2').is_file():
        print('Ambertools antechamber stage: converting to .mol2 and generating charges')
        # fixme - does not throw errors? try to make it throw an error
        subprocess.run([ambertools_bin / 'antechamber', '-i', left_ligand, '-fi', 'pdb',
                                            '-o', 'left.mol2', '-fo', 'mol2',
                                            '-c', charge_type, '-at', atom_type,
                                            '-nc', str(net_charge), '-dr', antechamber_dr],
                                   **subprocess_kwargs)

    # prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
    if not (workplace_root / 'right.mol2').is_file():
        print('Ambertools antechamber stage: converting to .mol2 and generating charges')
        subprocess.run([ambertools_bin / 'antechamber', '-i', right_ligand, '-fi', 'pdb',
                                    '-o', 'right.mol2', '-fo', 'mol2',
                                    '-c', charge_type, '-at', atom_type,
                                    '-nc', str(net_charge), '-dr', antechamber_dr], # FIXME - revisit the doctor
                                   **subprocess_kwargs)

    if use_original_coor:
        print(f'Copying coordinates from {left_ligand} and {right_ligand} since antechamber changes them slightly')
        # copy the files before applying the coordinates
        shutil.copy(workplace_root / 'left.mol2', workplace_root / 'left_before_COOR.mol2')
        shutil.copy(workplace_root / 'right.mol2', workplace_root / 'right_before_COOR.mol2')
        set_coor_from_ref(workplace_root / 'left.mol2', workplace_root / left_ligand, by_atom_name=True)
        set_coor_from_ref(workplace_root / 'right.mol2', workplace_root / right_ligand, by_atom_name=True)
        # set_coor_from_ref(workplace_root / 'left.mol2', workplace_root / left_ligand, by_index=True)
        # set_coor_from_ref(workplace_root / 'right.mol2', workplace_root / right_ligand, by_index=True)

    # rename the atom names to ensure they are unique across the two molecules
    renameAtomNamesUnique(workplace_root / 'left.mol2', workplace_root / 'right.mol2')

    # superimpose the two topologies
    suptop, mda_l1, mda_l2 = getSuptop(workplace_root / 'left.mol2', workplace_root / 'right.mol2',
                                       align_molecules=align_molecules,
                                       pair_charge_atol=atom_pair_q_atol,
                                       manual_match=manually_matched, force_mismatch=force_mismatch,
                                       net_charge_threshold=net_charge_threshold)

    # save the superimposition results
    left_right_matching_json = workplace_root / 'joint_meta_fep.json'
    save_superimposition_results(left_right_matching_json, suptop)
    # hybrid .pdb
    write_dual_top_pdb(workplace_root / 'morph.pdb', mda_l1, mda_l2, suptop)
    # hybrid .mol2
    hybrid_mol2 = workplace_root / 'morph.mol2'
    write_merged(suptop, hybrid_mol2)

    # generate the functional forms
    print('Ambertools parmchk2 generating .frcmod')
    left_chk2 = subprocess.run([ambertools_bin / 'parmchk2', '-i', 'left.mol2', '-o', 'left.frcmod',
                                 '-f', 'mol2', '-s', atom_type], **subprocess_kwargs)
    right_chk2 = subprocess.run([ambertools_bin / 'parmchk2', '-i', 'right.mol2', '-o', 'right.frcmod',
                                 '-f', 'mol2', '-s', atom_type], **subprocess_kwargs)

    # join the .frcmod files
    left_frcmod = workplace_root / 'left.frcmod'
    right_frcmod = workplace_root / 'right.frcmod'
    hybrid_frcmod = workplace_root / 'morph.frcmod'
    join_frcmod_files2(left_frcmod, right_frcmod, hybrid_frcmod)

    # if the hybrid .frcmod needs new terms between the appearing/disappearing atoms, insert dummy ones
    updated_frcmod_content = check_hybrid_frcmod(hybrid_mol2, hybrid_frcmod, amber_forcefield,
                                                 ambertools_bin, ambertools_script_dir, cwd=workplace_root)
    with open(hybrid_frcmod, 'w') as FOUT:
        FOUT.write(updated_frcmod_content)

    # prepare
    shutil.copy(namd_script_dir / "check_namd_outputs.py", workplace_root)
    shutil.copy(namd_script_dir / "ddg.py", workplace_root)


    ##########################################################
    # ------------------   Ligand ----------------------------
    prepare_inputs(workplace_root, directory='lig',
                   protein=None,
                   hybrid_mol2=hybrid_mol2,
                   hybrid_frc=hybrid_frcmod,
                   left_right_mapping=left_right_matching_json,
                   namd_script_loc=namd_script_dir,
                   scripts_loc=script_dir,
                   tleap_in='leap_ligand.in',
                   protein_ff=amber_forcefield,
                   net_charge=net_charge,
                   ambertools_script_dir=ambertools_script_dir,
                   subprocess_kwargs=subprocess_kwargs,
                   ambertools_bin=ambertools_bin,
                   namd_prod=namd_prod
                   )

    ##########################################################
    # ------------------ Complex  ----------------------------
    # calculate the charges of the protein (using ambertools)
    protein_net_charge = get_protein_net_charge(workplace_root, protein_filename,
                           ambertools_bin, ambertools_script_dir / 'solv_prot.in',
                           subprocess_kwargs, amber_forcefield)

    prepare_inputs(workplace_root, directory='complex',
                   protein=protein_filename,
                   hybrid_mol2=hybrid_mol2,
                   hybrid_frc=hybrid_frcmod,
                   left_right_mapping=left_right_matching_json,
                   namd_script_loc=namd_script_dir,
                   scripts_loc=script_dir,
                   tleap_in='leap_complex.in',
                   protein_ff=amber_forcefield,
                   net_charge=net_charge + protein_net_charge,
                   ambertools_script_dir=ambertools_script_dir,
                   subprocess_kwargs=subprocess_kwargs,
                   ambertools_bin=ambertools_bin,
                   namd_prod=namd_prod)

    print('TIES 20 Finished')

if __name__ == '__main__':
    command_line_script()