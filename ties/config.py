import os
import sys
import pathlib
import csv

import MDAnalysis
from ties.helpers import load_MDAnalysis_atom_group, ArgparseChecker

class Config:

    def __init__(self, **kwargs):
        # set the path to the scripts
        self.code_root = pathlib.Path(os.path.dirname(__file__))

        # scripts/input files
        self.script_dir = self.code_root / 'scripts'
        self.namd_script_dir = self.script_dir / 'namd'
        self.ambertools_script_dir = self.script_dir / 'ambertools'
        self.tleap_check_protein = self.ambertools_script_dir / 'check_prot.in'
        self.vmd_vis_script = self.script_dir / 'vmd' / 'vis_morph.vmd'

        self._workdir = None
        self._antechamber_dr = None
        self._ambertools_home = None

        self._protein = None

        self._ligand_net_charge = None
        self._atom_pair_q_atol = 0.1
        self._net_charge_threshold = 0.1
        self._redistribute_q_over_unmatched = True
        self._allow_disjoint_components = False
        # use only the element in the superimposition rather than the specific atom type
        self._use_element = False
        self._use_element_in_superimposition = True

        # coordinates
        self._align_molecules_using_mcs = False
        self._use_original_coor = False
        self._coordinates_file = None

        self._ligand_files = None
        self._manually_matched_atom_pairs = None
        self._ligands_contain_q = None

        self._ligand_tleap_in = None

        self._superimposition_starting_pair = None

        self._protein_ff = None
        self._ligand_ff = 'leaprc.gaff'
        self._ligand_ff_name = 'gaff'

        # MD/NAMD production input file
        self._md_engine = 'namd'
        self._lambda_rep_dir_tree = False

        # experimental
        self._use_hybrid_single_dual_top = False
        self._ignore_charges_completely = False

        # assign all the initial configuration values
        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    def workdir(self):
        if self._workdir is None:
            self._workdir = pathlib.Path(os.getcwd()) / 'ties'
            self._workdir.mkdir(exist_ok=True)

        return self._workdir

    @workdir.setter
    def workdir(self, cwd):
        if cwd is not None:
            # user provided
            self._workdir = cwd.absolute()
        else:
            # current directory
            self._workdir = pathlib.Path(os.getcwd()) / 'ties20'
            self._workdir.mkdir(exist_ok=True)
        print(f'Working Directory: {self._workdir}')

    # --------------- general
    @property
    def protein(self):
        return self._protein

    @protein.setter
    def protein(self, path):
        if path is None:
            print('No protein was select. Skipping protein dG.')
            return

        #  can be loaded by MDAnalysis?
        print(f'Trying to open the protein file {path} with MDAnalysis..')
        load_MDAnalysis_atom_group(path)
        self._protein = path

    @property
    def ligand_files(self):
        return self._ligand_files

    @ligand_files.setter
    def ligand_files(self, files):
        if len(files) < 1:
            print('Please supply at least one ligand file with -l (--ligands). E.g. -l file1.pdb file2.pdb')
            sys.exit(1)

        # check if the ligands have the same extension
        if len({l.suffix for l in files}) != 1:
            # fixme - does it actually matter? you parse it anyway?
            print('Error: ligands (-l) have different extensions. '
                  'Please ensure all ligands have the same extension')
            sys.exit()

        # files must have unique names
        filenames = [l.stem for l in files]
        if len(filenames) != len(set(filenames)):
            print('Error: some ligand (-l) names are the same. '
                  'Error: ensure your ligands have unique names. ')
            sys.exit()

        # verify the files with MDAnalysis if possible
        for ligand in files:
            if ligand.suffix.lower() in ('.ac', '.prep'):
                print(f'Cannot verify .ac/.prep {ligand} with MDAnalysis. Skipping. ')
            else:
                print(f'Trying to open the ligand {ligand} with MDAnalysis..')
                ligand_universe = load_MDAnalysis_atom_group(ligand)
                # there should be one residue
                if len(ligand_universe.residues.resnames) > 1:
                    print(f'Warning: more than one residue name detected in the ligand {ligand}')

        # TODO - ensure that it is a connected component and there is no extra atoms

        self._ligand_files = files

    # --------------- ambertools
    @property
    def ambertools_home(self):
        return self._ambertools_home

    @property
    def ambertools_antechamber(self):
        return self._ambertools_home / 'bin' / 'antechamber'

    @property
    def ambertools_parmchk2(self):
        return self._ambertools_home / 'bin' / 'parmchk2'

    @property
    def ambertools_tleap(self):
        return self._ambertools_home / 'bin' / 'tleap'

    @ambertools_home.setter
    def ambertools_home(self, path):
        # ambertools path given?
        if path is None:
            # otherwise check env paths
            if os.getenv('AMBERHOME'):
                path = pathlib.Path(os.getenv('AMBERHOME'))
            elif os.getenv('AMBER_PREFIX'):
                path = pathlib.Path(os.getenv('AMBER_PREFIX'))
            else:
                print('Error: Cannot find ambertools. $AMBERHOME and $AMBER_PREFIX are empty')
                print('Option 1: source your ambertools script amber.sh')
                print('Option 2: specify manually the path to amberhome with -ambertools option')
                sys.exit()

        self._ambertools_home = path

    @property
    def antechamber_dr(self):
        return self._antechamber_dr

    @antechamber_dr.setter
    def antechamber_dr(self, bool):
        # using ambertools antechamber dr mode?
        if bool is True:
            self._antechamber_dr = 'yes'
        else:
            self._antechamber_dr = 'no'
        print(f'Antechamber dr: {self._antechamber_dr}')

    @property
    def ligand_net_charge(self):
        if self._ligand_net_charge is None:
            print(f'WARNING: Ligand net charge not provided (-nc) by the user. Assuming 0. ')
            self._ligand_net_charge = 0

        return self._ligand_net_charge

    @ligand_net_charge.setter
    def ligand_net_charge(self, net_charge):
        if net_charge is None:
            print('Ligand net charge was not supplied (-nc). Only limited functionalities will be available. ')

        self._ligand_net_charge = net_charge
        print(f'Ligand net charge (-nc): {net_charge}')

    @property
    def coordinates_file(self):
        return self._coordinates_file

    @coordinates_file.setter
    def coordinates_file(self, file):
        if file is None:
            return

        print(f'MDAnalysis: verifying the coordinate file {file}')
        load_MDAnalysis_atom_group(file)
        # fixme - warn if the atom names are not uniq, warn if there is more than one residue, no water, etc
        self._coordinates_file = file

    @property
    def atom_pair_q_atol(self):
        return self._atom_pair_q_atol

    @atom_pair_q_atol.setter
    def atom_pair_q_atol(self, atol):
        self._atom_pair_q_atol = atol
        print(f'The maximum acceptable difference in charge between two paired atoms: {atol:.2f}')

    @property
    def net_charge_threshold(self):
        return self._net_charge_threshold

    @net_charge_threshold.setter
    def net_charge_threshold(self, threshold):
        self._net_charge_threshold = threshold
        print(f'Using MCS net charge difference threshold of {threshold}')

    @property
    def ignore_charges_completely(self):
        return self._ignore_charges_completely

    @ignore_charges_completely.setter
    def ignore_charges_completely(self, bool):
        self._ignore_charges_completely = bool
        if bool:
            print('Ignoring the charges. ')
            self.redistribute_q_over_unmatched = False

    # superimposition
    @property
    def allow_disjoint_components(self):
        return self._allow_disjoint_components

    @allow_disjoint_components.setter
    def allow_disjoint_components(self, boolean):
        self._allow_disjoint_components = boolean
        print(f'Allowing disjoint components: {self._allow_disjoint_components}')

    @property
    def use_element_in_superimposition(self):
        return self._use_element_in_superimposition

    @use_element_in_superimposition.setter
    def use_element_in_superimposition(self, bool):
        self._use_element_in_superimposition = bool
        if bool:
            print('Warning. Ignoring the specific atom types during superimposition. '
                  'The results should not be used in production simulations.')

    @property
    def align_molecules_using_mcs(self):
        return self._align_molecules_using_mcs

    @align_molecules_using_mcs.setter
    def align_molecules_using_mcs(self, boolean):
        # align the coordinates in ligZ to the ligA using the MCS
        self._align_molecules_using_mcs = boolean
        # fixme - should be using the MCS before charges change
        print(f'Will align the coordinates using the final MCS: {boolean}')

    @property
    def use_original_coor(self):
        return self._use_original_coor

    @use_original_coor.setter
    def use_original_coor(self, boolean):
        # use the original coords because antechamber can change them slightly
        self._use_original_coor = boolean

    @property
    def ligands_contain_q(self):
        if self._ligand_files is None:
            raise ValueError('The variable ligands_contain_q is used but _ligand_files are not configured. ')

        # does the file type contain charges?
        ligand_ext = self._ligand_files[0].suffix.lower()
        if self._ligands_contain_q is True:
            if ligand_ext in {'.mol2', '.ac'}:
                self._ligands_contain_q = True
            else:
                print('ERROR: If charges are provided with the ligands, '
                      'the filetypes .mol2 or .ac have to be used.')
                sys.exit(1)
        elif self._ligands_contain_q is False:
            self._ligands_contain_q = False
        elif self._ligands_contain_q is None:
            # determine whether charges are provided using the file extensions
            if ligand_ext in {'.mol2', '.ac', '.prep'}:
                self._ligands_contain_q = True
                print('Assuming that charges are provided based on the filetype .ac/.mol2')
            elif ligand_ext == '.pdb':
                self._ligands_contain_q = False
                print('Assuming that charges are not provided in the given .pdb ligand files. ')
            else:
                print(f'Error: unrecognized file type {ligand_ext}. ')
                sys.exit(1)

        # fixme?
        if self._ligands_contain_q:
            # leave charges the way they are in the files
            # TODO - ensure that when antechamber_charge_type is accessed ,this function is used? implement in another
            self.antechamber_charge_type = []

        print(f'Ligand files already contain charges: {self._ligands_contain_q}')

        return self._ligands_contain_q

    @ligands_contain_q.setter
    def ligands_contain_q(self, boolean):
        self._ligands_contain_q = boolean

    @property
    def superimposition_starting_pair(self):
        return self._superimposition_starting_pair

    @superimposition_starting_pair.setter
    def superimposition_starting_pair(self, value):
        if value == None:
            self._superimposition_starting_pair = None
        else:
            atom_name_disappearing, atom_name_appearing = value.split('-')
            self. _superimposition_starting_pair = (atom_name_disappearing, atom_name_appearing)

    @property
    def manually_matched_atom_pairs(self):
        return self._manually_matched_atom_pairs

    @manually_matched_atom_pairs.setter
    def manually_matched_atom_pairs(self, file_or_pairs):
        # A user-provided list of pairs that should be matched
        # This functionality works only if there are two ligands
        if file_or_pairs is None:
            self._manually_matched_atom_pairs = []
            return
        if self._ligand_files is None:
            raise ValueError('Wrong use of Config class. Please set the ._ligand_files attribute first. ')
        elif len(self._ligand_files) != 2:
            raise ValueError('Error: Too many ligands. '
                             'Manually matching atom pairs works only with 2 ligands. ')

        # TODO
        # pair_format: C1-C7 or C1-C7,C2-C8
        # if different format, it is likely a file, check file

        manually_matched = []
        if file_or_pairs is not None and pathlib.Path(file_or_pairs).is_file():
            with open(file_or_pairs) as IN:
                for left_atom, right_atom in csv.reader(IN, delimiter='-'):
                    manually_matched.append((left_atom.strip(), right_atom.strip()))

        if len(manually_matched) > 1:
            raise NotImplementedError('Currently only one atom pair can be matched - others were not tested')

        self._manually_matched_atom_pairs = manually_matched

    @property
    def protein_ff(self):
        return self._protein_ff

    @protein_ff.setter
    def protein_ff(self, ff):
        self._protein_ff = ff
        print(f'Protein force field name: {self._protein_ff}')

    @property
    def md_engine(self):
        return self._md_engine

    @md_engine.setter
    def md_engine(self, value):
        supported = ['NAMD(2)', 'NAMD3', 'OpenMM']
        if type(value) == str and value.lower() == 'namd':
            self._md_engine = value
        elif type(value) == str and value.lower() == 'namd3':
            self._md_engine = value
        elif type(value) == str and value.lower() == 'openmm':
            self._md_engine = value
        else:
            print('Unknown engine {}. Supported engines {}'.format(value, supported))
            # check if it is a bool value
            response = ArgparseChecker.str2bool(value)
            print(f'Generating files for an MD engine: {response}')
            self._md_engine = response
        print(f'MD Engine: {value}')

    @property
    def lambda_rep_dir_tree(self):
        return self._lambda_rep_dir_tree

    @lambda_rep_dir_tree.setter
    def lambda_rep_dir_tree(self, value):
        self._lambda_rep_dir_tree = value

    @property
    def ligand_ff(self):
        return self._ligand_ff

    @property
    def ligand_ff_name(self):
        return self._ligand_ff_name

    @ligand_ff_name.setter
    def ligand_ff_name(self, atom_type):
        # save also the atom type
        if atom_type == 'gaff':
            # they both use the same ff
            self._ligand_ff = 'leaprc.gaff'
        elif atom_type == 'gaff2':
            self._ligand_ff = 'leaprc.gaff2'
        else:
            raise ValueError('Argument -lff cannot be anything else but "gaff" or "gaff2". ')

        self._ligand_ff_name = atom_type

    @property
    def redistribute_q_over_unmatched(self):
        return self._redistribute_q_over_unmatched

    @redistribute_q_over_unmatched.setter
    def redistribute_q_over_unmatched(self, boolean):
        # Redistribute charge imbalances created due to atom-pair averaging
        self._redistribute_q_over_unmatched = boolean
        print(f'Distribute the introduced charge disparity in the alchemical region: '
              f'{self._redistribute_q_over_unmatched}')

    @property
    def use_hybrid_single_dual_top(self):
        # hybrid single dual topology
        return self._use_hybrid_single_dual_top

    @use_hybrid_single_dual_top.setter
    def use_hybrid_single_dual_top(self, boolean):
        if boolean:
            self._use_hybrid_single_dual_top = True
            self._complex_tleap_in = 'leap_complex_sdtop.in'
            if self._ignore_charges_completely != True:
                raise Exception('Charges have to be ignored completely when using hybrid single-dual topology.')
        else:
            self._complex_tleap_in = 'leap_complex.in'

        self._use_hybrid_single_dual_top = boolean

    @property
    def ligand_tleap_in(self):
        if self._ligand_tleap_in is not None:
            # return the user provided filename
            return self._ligand_tleap_in
        if self.use_hybrid_single_dual_top:
            # return the default option for the hybrid
            return 'leap_ligand_sdtop.in'

        # return the default
        return 'leap_ligand.in'

    @property
    def complex_tleap_in(self):
        if self._complex_tleap_in is None:
            # assume that hybrid single-dual topology is not used
            # this will initiate the standard leap_ligand.in
            self.use_hybrid_single_dual_top = False
        return self._complex_tleap_in

    @staticmethod
    def get_element_map():
        # Get the mapping of atom types to elements
        element_map_filename = pathlib.Path(os.path.dirname(__file__)) / 'data' / 'element_atom_type_map.txt'
        # remove the comments lines with #
        lines = filter(lambda l: not l.strip().startswith('#') and not l.strip() == '', open(element_map_filename).readlines())
        # convert into a dictionary

        element_map = {}
        for line in lines:
            element, atom_types = line.split('=')

            for atom_type in atom_types.split():
                element_map[atom_type.strip()] = element.strip()

        return element_map

    # fixme - this should be determined at the location where it is relevant rather than here in the conf
    # antechamber parameters, by default compute AM1-BCC charges
    antechamber_charge_type = ['-c', 'bcc']
