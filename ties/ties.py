#!/usr/bin/env python3
"""
Exposes a basic terminal interface to TIES 20.

Load two ligands, run the topology superimposer, and then using the results, generate the NAMD input files.
"""
import time
import itertools

from ties.generator import *
from ties.helpers import *
from ties.ligand import Ligand
from ties.morph import Morph
from ties.ligandmap import LigandMap
from ties.config import Config


def command_line_script():
    parser = argparse.ArgumentParser(description='TIES 20')
    parser.add_argument('action', metavar='command', type=str,
                        help='Action to be performed. E.g. "ties rename, ties create, .." ')
    parser.add_argument('-l', '--ligands', metavar='files ', dest='ligands',
                        nargs='+',
                        type=ArgparseChecker.existing_file, #required=True,
                        help='A list of uniquely named ligand files. '
                             'If more than 2 ligands are provided, '
                             'Lead Optimisation Mapping (TIES MAP) will be used.')
    parser.add_argument('-dir', '--ties-output-dir', metavar='directory', dest='tiesdir',
                        type=Path, required=False,
                        help='If not provided, "ties20" directory will be created in the current directory. ')
    parser.add_argument('-nc', '--ligand-net-charge', metavar='integer', dest='ligand_net_charge',
                        type=int, required=False,
                        help='An integer representing the net charge of each ligand. '
                             'All ligands must have the same charge.')
    parser.add_argument('-p', '--protein', metavar='file', dest='protein',
                        type=ArgparseChecker.existing_file, required=False,
                        help='The protein file')
    parser.add_argument('-qtol', '--q-pair-tolerance', metavar='decimal', dest='qtol',
                        type=float, required=False, default=0.1,
                        help='The maximum difference in charge between any two paired atoms (electron unit). '
                             'Default 0.1e')
    parser.add_argument('-netqtol', '--q-net-tolerance', metavar='decimal', dest='net_charge_threshold',
                        type=float, required=False, default=0.1,
                        help='The maximum difference in charge between two ligands. '
                             'Default 0.1e')
    parser.add_argument('-use-input-q', '--use-ligand-charges', metavar='boolean', dest='ligands_have_q',
                        type=ArgparseChecker.str2bool, required=False, default=None,
                        help='Use charges provided with the ligand files (with .mol2). '
                             'Default: If .mol2 is given, using the given charges will be attempted. ' # fixme - test if this agrees with the files
                             'Default: If .pdb is given, then BCC charges are computed with antechamber. ')
    parser.add_argument('-align-mcs', '--align-ligands-mcs', metavar='boolean', dest='align_mcs',
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
    parser.add_argument('-match', '--manual-match', metavar='file or pair', dest='manually_matched_pairs',
                        type=Path, required=False, # fixme - how multiple pairs affect it?
                        help='A path to a file that contains atom-pairs (one per line).'
                             'Each pair should be separated by dash "-" character.'
                             'For example a line with CL2-BR2 means that CL2 atom will be transformed to BR2 atom.')
    parser.add_argument('-mismatch', '--manual-mismatch', metavar='file or pair', dest='manually_mismatched_pairs',
                        type=Path, required=False, # fixme - implement
                        help='A path to a file that contains atom-pairs (one per line).'
                             'Each pair should be separated by dash "-" character.'
                             'For example a line with CL2-BR2 means that CL2 atom cannot be matched to BR2 atom.'
                             'Note that this option might have undesired consequences. ')
    parser.add_argument('-redist-q', '--redistribute-q-over-unmatched', metavar='boolean',
                        dest='redistribute_charges_over_unmatched',
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
    parser.add_argument('-pff', '--protein-ff', metavar='str', dest='protein_forcefield',
                        type=str, required=False, default='leaprc.protein.ff14SB',
                        help='This is a temporary solution. E.g. "leaprc.protein.ff14SB", '
                             ', "leaprc.ff99SBildn", etc.')
    parser.add_argument('-lff', '--ligand-ff', metavar='str', dest='ligand_ff_name',
                        type=str, required=False, default='gaff',
                        help='Either "gaff" or "gaff2"')
    parser.add_argument('-namd_prod', '--namd-prod', metavar='file', dest='namd_prod',
                        type=str, required=False, default='prod.namd',
                        help='This is a temporary solution. The name of the file to be used for the production. ')
    parser.add_argument('-md', '--md-engine', metavar='bool', dest='md_engine',
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

    args = parser.parse_args()

    # ------------------ Configuration and checks
    config = Config()
    config.workdir = args.tiesdir
    # ambertools
    config.ambertools_home = args.ambertools_home
    config.antechamber_dr = args.antechamber_dr
    # TIES 20 settings
    # charges
    config.ligand_net_charge = args.ligand_net_charge
    config.atom_pair_q_atol = args.qtol
    config.net_charge_threshold = args.net_charge_threshold
    config.redistribute_q_over_unmatched = args.redistribute_charges_over_unmatched
    config.ignore_charges_completely = args.ignore_charges_completely
    # coordinates
    config.align_molecules = args.align_mcs
    config.coordinates_file = args.coordinates_file
    # superimposition
    config.allow_disjoint_components = args.allow_disjoint_components
    config.use_element_in_superimposition = args.use_element
    # check the protein and the ligands
    config.protein = args.protein
    config.set_ligand_files(args.ligands, args.ligands_have_q,
                            args.manually_matched_pairs, args.manually_mismatched_pairs)
    config.protein_ff = args.protein_forcefield
    # ligand_ff_name is also referred to as "atom type", e.g. "gaff"
    config.set_ligand_ff(args.ligand_ff_name)
    config.use_hybrid_single_dual_top = args.hybrid_single_dual_top
    config.md_engine = args.md_engine
    config.lambda_rep_dir_tree = args.lambda_rep_dir_tree

    # TIES
    # create ligands
    ligands = [Ligand(lig, config) for lig in args.ligands]

    command = args.action
    if command == 'rename':
        # fixme - redo renaming?
        # this case assumes that there are only two ligands given
        if len(ligands) != 2:
            print('ERROR: to rename atoms to be unique across two ligands, '
                  'you have to provide exactly two ligands with -l.'
                  'E.g. ties rename -l init.mol2 final.mol2')
            sys.exit()
        print('Atom names will be renamed to ensure that the atom names are unique across the two molecules.')
        morph = Morph(ligands[0], ligands[1], config)
        morph.unique_atomres_names()
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

    # TIES
    # ensure that each atom name is unique
    [lig.make_atom_names_unique() for lig in ligands]

    # prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
    [lig.antechamber_prepare_mol2(config.ligand_ff_name, config.ligand_net_charge,
                                  config.antechamber_dr, config.antechamber_charge_type)
        for lig in ligands]

    # when the user provides charges, BCC and minimisation is not carried out, so the coordinates are correct
    if config.use_original_coor and not config.use_provided_charges:
        # todo - make this part of the lig.antechamber_prepare_mol2
        print(f'Copying coordinates from {ligands[0]} and {ligands[1]} since antechamber changes them slightly')
        # copy the files before applying the coordinates
        ligands[0].set_coor_from_ref(by_atom_name=True)

    # generate all pairings
    morphs = [Morph(ligA, ligZ, config) for ligA, ligZ in itertools.combinations(ligands, r=2)]

    # superimpose the two topologies
    start_time = time.time()
    for morph in morphs:
        print(f'Next ligand pair: {morph.internal_name}')
        # rename the atom names to ensure they are unique across the two molecules
        # we need to execute our pipeline for every pair and create the directories
        # which means we'll end up with a set of pairs
        morph.unique_atomres_names()

        # optimisation: check if the .json was created before, and use that to create the starting pair
        manually_matched = morph.check_json_file()
        if manually_matched is not None:
            print(f'Optimisation: reusing the preexisting .json file for {morph} as the starting point. ')

        morph.get_suptop(morph,
                         pair_charge_atol=config.atom_pair_q_atol,
                         manual_match=config.manually_matched_atom_pairs,
                         force_mismatch=None,
                         net_charge_threshold=config.net_charge_threshold,
                         redistribute_charges_over_unmatched=config.redistribute_q_over_unmatched,
                         ignore_charges_completely=config.ignore_charges_completely,
                         disjoint_components=config.allow_disjoint_components,
                         use_only_gentype=config.use_element_in_superimposition,
                         )

        # save meta data
        morph.write_summary_json()
        morph.write_pdb(hybrid_single_dual_top=config.use_hybrid_single_dual_top)
        morph.write_hybrid_mol2()
    print(f'Compared ligands to each other in: {time.time() - start_time:.1f} s')

    if config.use_hybrid_single_dual_top:
        raise NotImplementedError('Hybrid single-dual not done with LOMAP feature. ')
        rewrite_mol2_hybrid_top('left.mol2', list(matching_info["single_top_matched"].keys()))
        rewrite_mol2_hybrid_top('right.mol2', list(matching_info["single_top_matched"].values()))

    print('Ambertools parmchk2 generating .frcmod for ligands')
    [lig.generate_frcmod(config.ambertools_parmchk2, config.ligand_ff_name) for lig in ligands]

    # join the .frcmod files for each pair
    print('Ambertools parmchk2 generating .frcmod for morphs')
    [morph.join_frcmod_files(config.ambertools_tleap, config.ambertools_script_dir,
                             config.protein_ff, config.ligand_ff) for morph in morphs]

    # transformation hunter
    if len(ligands) == 2:
        selected_morphs = morphs
    else:
        lm = LigandMap(ligands, morphs)
        lm.generate_map()
        lm.print_map()
        lm.generate_graph()
        # selected_morphs = lm.traveling_salesmen()
        selected_morphs = lm.kruskal()

    ##########################################################
    # ------------------   Ligand ----------------------------
    # fixme - switch to morph.prepare_inputs instead.
    # this way we could reuse a lot of this information
    for morph in selected_morphs:
        prepare_inputs(morph,
                       config.workdir / 'lig',
                       protein=None,
                       namd_script_loc=config.namd_script_dir,
                       scripts_loc=config.script_dir,
                       tleap_in=config.ligand_tleap_in,
                       protein_ff=config.protein_ff,
                       ligand_ff=config.ligand_ff,
                       net_charge=config.ligand_net_charge,
                       ambertools_script_dir=config.ambertools_script_dir,
                       ambertools_tleap=config.ambertools_tleap,
                       namd_prod=config.namd_prod,
                       hybrid_topology=config.use_hybrid_single_dual_top,
                       vmd_vis_script=config.vmd_vis_script,
                       md_engine=config.md_engine,
                       lambda_rep_dir_tree=config.lambda_rep_dir_tree,
                       )
        print(f'Ligand {morph} directory populated successfully')

    ##########################################################
    # ------------------ Complex  ----------------------------
    # calculate the charges of the protein (using ambertools)
    if config.protein is not None:
        protein_net_charge = get_protein_net_charge(config.workdir, config.protein.absolute(),
                               config.ambertools_tleap, config.tleap_check_protein,
                               config.protein_ff)
        print(f'Protein net charge: {protein_net_charge}')

        for morph in selected_morphs:
            prepare_inputs(morph,
                           config.workdir / 'complex',
                           protein=config.protein,
                           namd_script_loc=config.namd_script_dir,
                           scripts_loc=config.script_dir,
                           tleap_in=config.complex_tleap_in,
                           protein_ff=config.protein_ff,
                           ligand_ff=config.ligand_ff,
                           net_charge=config.ligand_net_charge + protein_net_charge,
                           ambertools_script_dir=config.ambertools_script_dir,
                           ambertools_tleap=config.ambertools_tleap,
                           namd_prod=config.namd_prod,
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