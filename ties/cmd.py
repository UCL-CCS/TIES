#!/usr/bin/env python3
"""
Exposes a terminal interface to TIES 20.
"""
import time
import itertools
import pathlib
import shutil
import sys

import ties.generator
from ties.helpers import *
import ties.config
import ties.ligand
import ties.pair
import ties.ligandmap


def command_line_script():
    parser = argparse.ArgumentParser(description='TIES 20')
    parser.add_argument('command', metavar='str', type=str,
                        help='Action to be performed. E.g. "rename, ties create, .." ')
    parser.add_argument('-l', '--ligands', metavar='files ', dest='ligands',
                        nargs='+',
                        type=ArgparseChecker.existing_file,
                        help='A list of uniquely named ligand files. '
                             'If more than 2 ligands are provided, '
                             'Lead Optimisation Mapping (TIES MAP) will be used.')
    parser.add_argument('-dir', '--ties-output-dir', metavar='directory', dest='workdir',
                        type=pathlib.Path, required=False,
                        help='If not provided, "ties20" directory will be created in the current directory. ')
    parser.add_argument('-nc', '--ligand-net-charge', metavar='integer', dest='ligand_net_charge',
                        type=int, required=False,
                        help='An integer representing the net charge of each ligand. '
                             'All ligands must have the same charge.')
    parser.add_argument('-p', '--protein', metavar='file', dest='protein',
                        type=ArgparseChecker.existing_file, required=False,
                        help='The protein file')
    parser.add_argument('-qtol', '--q-pair-tolerance', metavar='decimal', dest='atom_pair_q_atol',
                        type=float, required=False, default=0.1,
                        help='The maximum difference in charge between any two paired atoms (electron unit). '
                             'Default 0.1e')
    parser.add_argument('-netqtol', '--q-net-tolerance', metavar='decimal', dest='net_charge_threshold',
                        type=float, required=False, default=0.1,
                        help='The maximum difference in charge between two ligands. '
                             'Default 0.1e')
    parser.add_argument('-use-input-q', '--use-ligand-charges', metavar='boolean', dest='ligands_contain_q',
                        type=ArgparseChecker.str2bool, required=False, default=None,
                        help='Use charges provided with the ligand files (with .mol2). '
                             'Default: If .mol2 is given, using the given charges will be attempted. ' # fixme - test if this agrees with the files
                             'Default: If .pdb is given, then BCC charges are computed with antechamber. ')
    parser.add_argument('-align-mcs', '--align-ligands-mcs', metavar='boolean', dest='align_molecules_using_mcs',
                        type=ArgparseChecker.str2bool, required=False, default=False,
                        help='Align the coordinates in the "END" ligand to the "INITIAL" ligand using '
                             'the generated maximum common substructure (MCS).')
    parser.add_argument('-ant-dr', '--antechamber-dr', metavar='boolean', dest='antechamber_dr',
                        type=ArgparseChecker.str2bool, required=False, default=False,
                        help='Antechamber dr is turned off by default. It is playing up with .mol2 files. '
                             'Please ensure that your input is valid when you turn off the antechamber dr. ')
    parser.add_argument('-ambertools', '--ambertools-home', metavar='path', dest='ambertools_home',
                        type=ArgparseChecker.ambertools_home, required=False,
                        help='Path to the home directory of ambertools '
                             '(the one that contains "bin" directory which contains "antechamber" file)')
    parser.add_argument('-match', '--manual-match', metavar='file or pair', dest='manually_matched_atom_pairs',
                        type=pathlib.Path, required=False,
                        help='A path to a file that contains atom-pairs (one per line).'
                             'Each pair should be separated by dash "-" character.'
                             'For example a line with CL2-BR2 means that CL2 atom will be transformed to BR2 atom.')
    parser.add_argument('-supseed', '--superimposition-starting-pair', metavar='pair',
                        dest='superimposition_starting_pair',
                        type=str, required=False,
                        help='A starting seed for the superimposition. '
                             'This is the pair from which the superimposition algorithm will start the '
                             'traversal to find the superimposition. '
                             'Example: "C1-C34" where C1 is in the disappearing molecule, '
                             'and C34 is in the appearing molecule. ')
    parser.add_argument('-mismatch', '--manual-mismatch', metavar='file or pair', dest='manually_mismatched_pairs',
                        type=pathlib.Path, required=False,  # fixme - implement
                        help='A path to a file that contains atom-pairs (one per line).'
                             'Each pair should be separated by dash "-" character.'
                             'For example a line with CL2-BR2 means that CL2 atom cannot be matched to BR2 atom.'
                             'Note that this option might have undesired consequences. ')
    parser.add_argument('-redist-q', '--redistribute-q-over-unmatched', metavar='boolean',
                        dest='redistribute_q_over_unmatched',
                        type=ArgparseChecker.str2bool, required=False, default=True,
                        help='Averaging the charges in a pair of matched atoms (see --q-pair-tolerance) '
                             'changes the overall molecule charge slightly. '
                             'Redistribute the lost/gained charges over the unmatched area '
                             'to make the charges equal. ')
    parser.add_argument('-hybrid-top', '--hybrid-singe-dual-top', metavar='boolean',
                        dest='hybrid_single_dual_top',
                        type=ArgparseChecker.str2bool, required=False, default=False,
                        help='Use hybrid single-dual topology in NAMD. See NAMD manual and '
                             'https://pubs.acs.org/doi/10.1021/acs.jcim.9b00362')
    parser.add_argument('-disjoint-allowed', '--disjoint-appearing-disappearing', metavar='boolean',
                        dest='allow_disjoint_components',
                        type=ArgparseChecker.str2bool, required=False, default=False,
                        help='Allow the molecules to be divided by "disappearing/appearing" atoms.')
    # temporary # fixme - list the FFs here
    parser.add_argument('-pff', '--protein-ff', metavar='str', dest='protein_ff',
                        type=str, required=False, default='leaprc.protein.ff14SB',
                        help='This is a temporary solution. E.g. "leaprc.protein.ff14SB", '
                             ', "leaprc.ff99SBildn", etc.')
    parser.add_argument('-lff', '--ligand-ff', metavar='str', dest='ligand_ff_name',
                        type=str, required=False, default='gaff',
                        help='Either "gaff" or "gaff2"')
    parser.add_argument('-md', '--md-engine', metavar='str', dest='md_engine',
                        type=str, required=False, default='namd',
                        help='Generate input files for the selected engine. '
                             'Use value "no" if you do not want the MD input to be generated. '
                             'Default is "namd". ')
    parser.add_argument('-dirtree', '--engine-tree', metavar='bool', dest='lambda_rep_dir_tree',
                        type=ArgparseChecker.str2bool, required=False, default=False,
                        help='Generate the directory tree structure for each lambda/replica directory. ')
    # allow to overwrite the coordinates
    parser.add_argument('-crd', '--coordinates', metavar='file', dest='coordinates_file',
                        type=ArgparseChecker.existing_file, required=False,
                        help='The protein file')
    parser.add_argument('-o', '--output-file', metavar='str', dest='output_filename',
                        type=str, required=False,
                        help='Where to save the output file. The extension is necessary. ')
    # dev tools
    parser.add_argument('-noq', '--ignore-charges', metavar='boolean', dest='ignore_charges_completely',
                        type=ArgparseChecker.str2bool, required=False, default=False,
                        help='Ignore charges throughout. This is mainly for debugging. ')
    parser.add_argument('-elements', '--compare-elements', metavar='boolean', dest='use_element',
                        type=ArgparseChecker.str2bool, required=False, default=False,
                        help='Ignore the specific atom types in the superimposition. Use only the elements. ')

    # initialise the config class
    args = parser.parse_args()
    command = args.command.lower()
    delattr(args, 'command')
    # assign all the parsed arguments to the Config class
    config = ties.config.Config(**args.__dict__)

    # create ligands
    ligands = [ties.ligand.Ligand(lig, config) for lig in args.ligands]

    if command == 'rename':
        # this case assumes that there are only two ligands given
        if len(ligands) != 2:
            print('ERROR: to rename atoms to be unique across two ligands, '
                  'you have to provide exactly two ligands with -l.'
                  'E.g. ties rename -l init.mol2 final.mol2')
            sys.exit()
        print('Atom names will be renamed to ensure that the atom names are unique across the two molecules.')
        pair = ties.pair.Pair(ligands[0], ligands[1], config)
        pair.make_atom_names_unique()
        sys.exit()
    elif command == 'mergecrd':
        print('Merging coordinates. Will add coordinates from an external file. ')
        if len(ligands) > 1:
            print('ERROR: too many ligands ligand (-l). '
                  'Use one ligand when assigning coordinates from another file.')
            sys.exit(1)
        elif config.coordinates_file is None:
            print('ERROR: no file with the coordinates found. Please add a file with coordinates (-crd). ')
            sys.exit(1)
        elif args.output_filename is None:
            print('ERROR: no output file provided for the mergecrd. Please use (-o), e.g. -o output.mol2 ')
            sys.exit(1)

        # assign coordinates
        ligands[0].overwrite_coordinates_with(config.coordinates_file, args.output_filename)
        sys.exit(0)
    elif command == 'analyse':
        raise NotImplementedError('analysis is not yet implemented in the main TIES')
    elif command != 'create':
        print('Please provide action (rename/create)')
        sys.exit()

    # prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
    [lig.antechamber_prepare_mol2(config.ligand_ff_name, config.ligand_net_charge,
                                  config.antechamber_dr, config.antechamber_charge_type)
        for lig in ligands]

    # generate all pairings
    pairs = [ties.pair.Pair(ligA, ligZ, config) for ligA, ligZ in itertools.combinations(ligands, r=2)]

    # superimpose the paired topologies
    start_time = time.time()
    for pair in pairs:
        print(f'Next ligand pair: {pair.internal_name}')
        # rename the atom names to ensure they are unique across the two molecules
        pair.make_atom_names_unique()

        # superimpose the topologies
        hybrid = pair.compute_suptop()

        # save meta data
        hybrid.write_summary_json()
        hybrid.write_pdb()
        hybrid.write_mol2()
    print(f'Compared ligands to each other in: {time.time() - start_time:.1f} s')

    if config.use_hybrid_single_dual_top:
        raise NotImplementedError('Hybrid single-dual not done with LOMAP feature. ')
        rewrite_mol2_hybrid_top('left.mol2', list(matching_info["single_top_matched"].keys()))
        rewrite_mol2_hybrid_top('right.mol2', list(matching_info["single_top_matched"].values()))

    print('Ambertools parmchk2 generating .frcmod for ligands')
    [lig.generate_frcmod(config.ambertools_parmchk2, config.ligand_ff_name) for lig in ligands]

    # join the .frcmod files for each pair
    print('Ambertools parmchk2 generating .frcmod for pairs')
    [pair.merge_frcmod_files(config.ambertools_tleap, config.ambertools_script_dir,
                              config.protein_ff, config.ligand_ff) for pair in pairs]

    # transformation hunter
    if len(ligands) == 2:
        selected_pairs = pairs
    else:
        lm = LigandMap(ligands, pairs)
        lm.generate_map()
        lm.print_map()
        lm.generate_graph()
        # selected_pairs = lm.traveling_salesmen()
        selected_pairs = lm.kruskal()

    ##########################################################
    # ------------------   Ligand ----------------------------
    # fixme - switch to pair.prepare_inputs instead.
    # this way we could reuse a lot of this information
    for pair in selected_pairs:
        ties.generator.prepare_inputs(pair,
                       config.workdir,
                       protein=None,
                       namd_script_loc=config.namd_script_dir,
                       scripts_loc=config.script_dir,
                       tleap_in=config.ligand_tleap_in,
                       protein_ff=config.protein_ff,
                       ligand_ff=config.ligand_ff,
                       net_charge=config.ligand_net_charge,
                       ambertools_script_dir=config.ambertools_script_dir,
                       ambertools_tleap=config.ambertools_tleap,
                       hybrid_topology=config.use_hybrid_single_dual_top,
                       vmd_vis_script=config.vmd_vis_script,
                       md_engine=config.md_engine,
                       lambda_rep_dir_tree=config.lambda_rep_dir_tree,
                       )
        print(f'Ligand {pair} directory populated successfully')

    ##########################################################
    # ------------------ Complex  ----------------------------
    if config.protein is not None:
        # calculate the charges of the protein (using ambertools)
        protein_net_charge = ties.generator.get_protein_net_charge(config.workdir, config.protein.absolute(),
                               config.ambertools_tleap, config.tleap_check_protein,
                               config.protein_ff)
        print(f'Protein net charge: {protein_net_charge}')

        for pair in selected_pairs:
            ties.generator.prepare_inputs(pair,
                           config.workdir / 'com',
                           protein=config.protein,
                           namd_script_loc=config.namd_script_dir,
                           scripts_loc=config.script_dir,
                           tleap_in=config.complex_tleap_in,
                           protein_ff=config.protein_ff,
                           ligand_ff=config.ligand_ff,
                           net_charge=config.ligand_net_charge + protein_net_charge,
                           ambertools_script_dir=config.ambertools_script_dir,
                           ambertools_tleap=config.ambertools_tleap,
                           hybrid_topology=config.use_hybrid_single_dual_top,
                           vmd_vis_script=config.vmd_vis_script,
                           md_engine=config.md_engine,
                           lambda_rep_dir_tree=config.lambda_rep_dir_tree,
                           )

    # prepare the post-analysis scripts
    shutil.copy(config.namd_script_dir / "check_namd_outputs.py", config.workdir)
    shutil.copy(config.namd_script_dir / "ddg.py", config.workdir)

    print('TIES 20 Finished')

if __name__ == '__main__':
    command_line_script()