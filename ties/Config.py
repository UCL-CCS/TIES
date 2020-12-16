import os
import sys
import pathlib
import csv

import MDAnalysis
from ties.helpers import load_MDAnalysis_atom_group

class Config:

    def __init__(self):
        # set the path to the scripts
        code_root = pathlib.Path(os.path.dirname(__file__))
        # scripts/input files
        self.script_dir = code_root / 'scripts'
        self.namd_script_dir = self.script_dir / 'namd'
        self.ambertools_script_dir = self.script_dir / 'ambertools'
        self.tleap_check_protein = self.ambertools_script_dir / 'check_prot.in'

        self._workplace_dir = None
        self._antechamber_dr = None
        self._ambertools_home = None

        self._protein = None

        self._ligand_net_charge = None
        self._atom_pair_q_atol = 0.1
        self._net_charge_threshold = 0.1
        self._redistribute_q_over_unmatched = True
        self._allow_disjoint_components = False

        # coordinates
        self._align_molecules = False
        self._use_original_coor = False

        self._ligand_files = None
        self._manually_matched_atom_pairs = None
        self._ligands_contain_q = None

        self._protein_ff = None
        self._ligand_ff = 'leaprc.gaff'
        self._ligand_ff_name = 'gaff'

        # NAMD production input file
        self._namd_prod = 'prod_2017.namd'  # Berendsen barostat ..
        print(f'The NAMD production file: {self._namd_prod}')

        # experimental
        self._use_hybrid_single_dual_top = False
        self._ignore_charges_completely = False

    @property
    def workdir(self):
        return self._workplace_dir

    @workdir.setter
    def workdir(self, cwd):
        if cwd is not None:
            # user provided
            self._workplace_dir = cwd
        else:
            # current directory
            self._workplace_dir = pathlib.Path(os.getcwd()) / 'ties20'
            self._workplace_dir.mkdir(exist_ok=True)
        print(f'Working Directory: {self._workplace_dir}')

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

    def set_ligand_files(self, files,
                         ligands_have_q=None,
                         manually_matched_pairs=None,
                         manually_mismatched_pairs=None):
        if len(files) < 2:
            print('Please supply at least two ligand files with -l (--ligands). E.g. -l file1.pdb file2.pdb')
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
            if ligand.suffix.lower() == '.ac':
                print(f'Cannot verify initially .ac file {ligand} with MDAnalysis. Skipping. ')
            else:
                print(f'Trying to open the ligand {ligand} with MDAnalysis..')
                ligand_universe = MDAnalysis.Universe(ligand)
                # there should be one residue
                if len(ligand_universe.residues.resnames) > 1:
                    print(f'Warning: more than one residue name detected in the ligand {ligand}')

        # TODO - ensure that it is a connected component and there is no extra atoms

        self._ligand_files = files

        # whether the ligands have their own charges
        self.ligands_contain_q = ligands_have_q
        # the manually matched pairs must be delivered after
        self.manually_matched_atom_pairs = manually_matched_pairs

        # fixme
        # todo
        # force_mismatch = []
        # if manually_mismatched_pairs is not None:
        #     with open(manually_mismatched_pairs) as IN:
        #         for left_atom, right_atom in csv.reader(IN, delimiter='-'):
        #             force_mismatch.append((left_atom, right_atom))
            # fixme - check that there is no overlap between match and mismatch lists

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
            raise ValueError('Ligand net charge (-nc) was not configured but is needed. ')
        return self._ligand_net_charge

    @ligand_net_charge.setter
    def ligand_net_charge(self, net_charge):
        if net_charge is None:
            print('Ligand net charge was not supplied (-nc). Only limited functionalities will be available. ')

        self._ligand_net_charge = net_charge
        print(f'Ligand net charge (-nc): {net_charge}')

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

    # superimposition
    @property
    def allow_disjoint_components(self):
        return self._allow_disjoint_components

    @allow_disjoint_components.setter
    def allow_disjoint_components(self, boolean):
        self._allow_disjoint_components = boolean
        print(f'Allowing disjoint components: {self._allow_disjoint_components}')

    @property
    def align_molecules(self):
        return self._align_molecules

    @align_molecules.setter
    def align_molecules(self, boolean):
        # align the coordinates in ligZ to the ligA using the MCS
        self._align_molecules = boolean
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
        return self._ligands_contain_q

    @ligands_contain_q.setter
    def ligands_contain_q(self, boolean):
        if self._ligand_files is None:
            raise ValueError('Wrong use of the Config class. Please set the ._ligand_files attribute first. ')

        # does the file type contain charges?
        ligand_ext = self._ligand_files[0].suffix.lower()
        if boolean is True:
            if ligand_ext in {'.mol2', '.ac'}:
                self._ligands_contain_q = True
            else:
                print('ERROR: If charges are provided with the ligands, '
                      'the filetypes .mol2 or .ac have to be used.')
                sys.exit(1)
        elif boolean is False:
            self._ligands_contain_q = False
        elif boolean is None:
            # determine whether charges are provided using the file extensions
            if ligand_ext in {'.mol2', '.ac'}:
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
            self.antechamber_charge_type = []

        print(f'Ligand files already contain charges: {self._ligands_contain_q}')

    @property
    def manually_matched_atom_pairs(self):
        return self._manually_matched_atom_pairs

    @manually_matched_atom_pairs.setter
    def manually_matched_atom_pairs(self, file_or_pairs):
        # A user-provided list of pairs that should be matched
        # This functionality works only if there are two ligands
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

    @property
    def protein_ff(self):
        return self._protein_ff

    @protein_ff.setter
    def protein_ff(self, ff):
        self._protein_ff = ff
        print(f'Protein force field name: {self._protein_ff}')

    # fixme - allow providing another .namd file
    @property
    def namd_prod(self):
        # fixme - add user support
        return self._namd_prod

    @property
    def ligand_ff(self):
        return self._ligand_ff

    @property
    def ligand_ff_name(self):
        return self._ligand_ff_name

    def set_ligand_ff(self, atom_type):
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
            self._ligand_tleap_in = 'leap_ligand_sdtop.in'
            self._complex_tleap_in = 'leap_complex_sdtop.in'
            self._ignore_charges_completely = True
        else:
            self._ignore_charges_completely = False
            self._ligand_tleap_in = 'leap_ligand.in'
            self._complex_tleap_in = 'leap_complex.in'

        self._use_hybrid_single_dual_top = boolean

    @property
    def ligand_tleap_in(self):
        if self._ligand_tleap_in is None:
            # assume that hybrid single-dual topology is not used
            # this will initiate the standard leap_ligand.in
            self.use_hybrid_single_dual_top = False
        return self._ligand_tleap_in

    @property
    def complex_tleap_in(self):
        if self._complex_tleap_in is None:
            # assume that hybrid single-dual topology is not used
            # this will initiate the standard leap_ligand.in
            self.use_hybrid_single_dual_top = False
        return self._complex_tleap_in

    # fixme - this should be determined at the location where it is relevant rather than here in the conf
    # antechamber parameters, by default compute AM1-BCC charges
    antechamber_charge_type = ['-c', 'bcc']
