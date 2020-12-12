#!/usr/bin/env python3
"""
Exposes a basic terminal interace to TIES 20.

Load two ligands, run the topology superimposer, and then using the results, generate the NAMD input files.
"""
import itertools
import csv

from ties.generator import *
from ties.helpers import *
from ties.Ligand import Ligand
from ties.Morph import Morph
from ties.LigandMap import LigandMap


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
    script_dir = code_root / Path('scripts')
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

    # generate all pairings
    morphs = [Morph(ligA, ligZ, workplace_root) for ligA, ligZ in itertools.combinations(ligands, r=2)]

    # superimpose the two topologies
    for morph in morphs:
        print(f'Next ligand pair: {morph.internal_name}')
        # rename the atom names to ensure they are unique across the two molecules
        # we need to execute our pipeline for every pair and create the directories
        # which means we'll end up with a set of pairs
        morph.unique_atomres_names()

        # TODO optimisation: check if the .json was created before, and use that to create the starting pair
        manually_matched = morph.check_json_file()

        suptop, mda_l1, mda_l2 = get_suptop(morph.renamed_ligA, morph.renamed_ligZ,
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

    # decide on which pairs to compare in order to obtain the full ranking
    if len(ligands) == 2:
        selected_morphs = morphs
    else:
        lm = LigandMap(ligands, morphs)
        lm.generate_map()
        lm.print_map()
        lm.generate_graph()
        selected_morphs = lm.traveling_salesmen()
        selected_morphs = lm.kruskal()

    ##########################################################
    # ------------------   Ligand ----------------------------
    # pick the right tleap instructions
    if use_hybrid_single_dual_top:
        ligand_tleap_in = 'leap_ligand_sdtop.in'
    else:
        ligand_tleap_in = 'leap_ligand.in'

    # fixme - at this point you'd know which pairs to set up
    # fixme - rather than using this, we should be able to have morph.prepare_inputs instead.
    # this way we could reuse a lot of this information
    for morph in selected_morphs:
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

        for morph in selected_morphs:
            prepare_inputs(morph, dir_prefix='complex',
                           protein=protein_filename,
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