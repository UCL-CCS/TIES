import os
import json
import copy
import subprocess
from pathlib import Path

import ties.helpers
from ties.generator import get_atoms_bonds_from_mol2
from ties.topology_superimposer import superimpose_topologies
import ties.config
import ties.ligand


class Pair():
    """
    A convenience class to help organise morphs.
    It offers functionality related to a pair of ligands (a transformation).

    Note that Morph actually in a way duplicates the concept of a SupTop.
    This distinction could be removed easily in the future.
    However, way to simplify SupTop might be worth some effort.
    """
    ROOT = Path('prep')
    FRCMOD_DIR_name = ROOT / Path('morph_frcmods')
    FRCMOD_TEST_DIR_name = FRCMOD_DIR_name / 'tests'
    UNIQUE_ATOM_NAMES_name = ROOT / Path('morph_unique_atom_names')

    def __init__(self, ligA, ligZ, config=None):
        # create a new config if it is not provided
        self.config = ties.config.Config() if config is None else config

        # tell Config about the ligands if necessary
        if self.config.ligands is None:
            self.config.ligands = [ligA, ligZ]

        # create ligands if they're just paths
        if isinstance(ligA, ties.ligand.Ligand):
            self.ligA = ligA
        else:
            self.ligA = ties.ligand.Ligand(ligA)
        if isinstance(ligZ, ties.ligand.Ligand):
            self.ligZ = ligZ
        else:
            self.ligZ = ties.ligand.Ligand(ligZ)

        # initialise the handles to the molecules that morph
        self.current_ligA = self.ligA.current
        self.current_ligZ = self.ligZ.current

        self.UNIQUE_ATOM_NAMES = self.config.workdir / Pair.UNIQUE_ATOM_NAMES_name
        self.FRCMOD_DIR = self.config.workdir / Pair.FRCMOD_DIR_name
        self.FRCMOD_TEST_DIR = self.config.workdir / Pair.FRCMOD_TEST_DIR_name

        self.internal_name = f'{self.ligA.internal_name}_{self.ligZ.internal_name}'
        self.mol2 = None
        self.pdb = None
        self.summary = None
        self.suptop = None
        self.mda_l1 = None
        self.mda_l2 = None
        self.distance = None

    def __repr__(self):
        if self.distance is None:
            return self.internal_name

        return f'{self.internal_name} ({self.distance:.2f})'

    def set_distance(self, value):
        self.distance = value

    def superimpose(self, **kwargs):

        # use MDAnalysis to load the files
        # fixme - move this to the Morph class instead of this place,
        # fixme - should not squash all messsages. For example, wrong type file should not be squashed
        leftlig_atoms, leftlig_bonds, rightlig_atoms, rightlig_bonds, mda_l1, mda_l2 = \
            get_atoms_bonds_from_mol2(self.current_ligA, self.current_ligZ,
                                      use_general_type=self.config.use_element_in_superimposition)
        # fixme - manual match should be improved here and allow for a sensible format.

        # empty lists count as None
        # if not force_mismatch:
        force_mismatch = None

        # assign
        # fixme - Ideally I would reuse the mdanalysis data for this,
        # MDAnalysis can use bonds if they are present - fixme
        # map atom IDs to their objects
        ligand1_nodes = {}
        for atomNode in leftlig_atoms:
            ligand1_nodes[atomNode.id] = atomNode
        # link them together
        for nfrom, nto, btype in leftlig_bonds:
            ligand1_nodes[nfrom].bind_to(ligand1_nodes[nto], btype)

        ligand2_nodes = {}
        for atomNode in rightlig_atoms:
            ligand2_nodes[atomNode.id] = atomNode
        for nfrom, nto, btype in rightlig_bonds:
            ligand2_nodes[nfrom].bind_to(ligand2_nodes[nto], btype)

        # fixme - this should be moved out of here,
        #  ideally there would be a function in the main interface for this
        manual_match = [] if self.config.manually_matched_atom_pairs is None else self.config.manually_matched_atom_pairs
        starting_node_pairs = []
        for l_aname, r_aname in manual_match:
            # find the starting node pairs, ie the manually matched pair(s)
            found_left_node = None
            for id, ln in ligand1_nodes.items():
                if l_aname == ln.name:
                    found_left_node = ln
            if found_left_node is None:
                raise ValueError(f'Manual Matching: could not find an atom name: "{l_aname}" in the left molecule')

            found_right_node = None
            for id, ln in ligand2_nodes.items():
                if r_aname == ln.name:
                    found_right_node = ln
            if found_right_node is None:
                raise ValueError(f'Manual Matching: could not find an atom name: "{r_aname}" in the right molecule')

            starting_node_pairs.append([found_left_node, found_right_node])

        if starting_node_pairs:
            print('Starting nodes will be used:', starting_node_pairs)

        # fixme - simplify to only take the mdanalysis as input
        suptops = superimpose_topologies(ligand1_nodes.values(), ligand2_nodes.values(),
                                         disjoint_components=self.config.allow_disjoint_components,
                                         net_charge_filter=True,
                                         pair_charge_atol=self.config.atom_pair_q_atol,
                                         net_charge_threshold=self.config.net_charge_threshold,
                                         redistribute_charges_over_unmatched=self.config.redistribute_q_over_unmatched,
                                         ignore_charges_completely=self.config.ignore_charges_completely,
                                         ignore_bond_types=True,
                                         ignore_coords=False,
                                         left_coords_are_ref=True,
                                         align_molecules=self.config.align_molecules_using_mcs,
                                         use_general_type=self.config.use_element_in_superimposition,
                                         # fixme - not the same ... use_element_in_superimposition,
                                         use_only_element=False,
                                         check_atom_names_unique=True, # fixme - remove?
                                         starting_pairs_heuristics=True, # fixme - add to config
                                         force_mismatch=self.config.manually_mismatched_pairs, # fixme - add to config
                                         starting_node_pairs=starting_node_pairs,
                                         ligand_l_mda=mda_l1, ligand_r_mda=mda_l2,
                                         starting_pair_seed=self.config.superimposition_starting_pair)

        assert len(suptops) == 1
        self.set_suptop(suptops[0], mda_l1, mda_l2)
        # attach the used config to the suptop
        suptops[0].config = self.config
        # attach the morph to the suptop
        suptops[0].morph = self

        return suptops[0]

    def set_suptop(self, suptop, mda_l1, mda_l2):
        self.suptop = suptop
        self.mda_l1 = mda_l1
        self.mda_l2 = mda_l2

    def make_atom_names_unique(self):
        """
        Ensure that each has a name that is unique to both ligands.

        Rename the ligand to ensure that no atom has the same name
        name atoms using the first letter (C, N, ..) and count them
        keep the names if possible (ie if they are already different)

        Resnames are set to "INI" and "FIN", this is useful for the hybrid dual topology

        fixme - rewrite this to be LigA.unique_names(LigZ) so that they can do it internally themselves?
        """
        # load both ligands
        left = ties.helpers.load_MDAnalysis_atom_group(self.ligA.current)
        right = ties.helpers.load_MDAnalysis_atom_group(self.ligZ.current)

        ligands_atom_names_overlap = len(set(right.atoms.names).intersection(set(left.atoms.names))) > 0

        if ligands_atom_names_overlap:
            print(f'Renaming right molecule ({self.ligZ.internal_name}) atom names because they are not unique')
            name_counter_L_nodes = ties.helpers.get_atom_names_counter(left.atoms)
            ties.helpers.rename_ligand(right.atoms, name_counter=name_counter_L_nodes)

        # rename the residue names to INI and FIN
        left.residues.resnames = ['INI']
        right.residues.resnames = ['FIN']

        # prepare the destination directory
        cwd = self.config.workdir / self.UNIQUE_ATOM_NAMES / f'{self.ligA.internal_name}_{self.ligZ.internal_name}'
        cwd.mkdir(parents=True, exist_ok=True)

        # save the updated atom names
        # fixme - write only when requested
        self.current_ligA = cwd / (self.ligA.internal_name + '.mol2')
        left.atoms.write(self.current_ligA)
        self.current_ligZ = cwd / (self.ligZ.internal_name + '.mol2')
        right.atoms.write(self.current_ligZ)

        # update the Ligands to use the renamed version

    def check_json_file(self):
        # An optimisation for testing
        # Check if the json was produced before, use it to
        # find the starting point to speed up the testing
        matching_json = self.config.workdir / f'fep_{self.ligA.internal_name}_{self.ligZ.internal_name}.json'
        if not matching_json.is_file():
            return None

        return [list(json.load(matching_json.open())['matched'].items())[0]]

    def merge_frcmod_files(self):
        """
        I tested the duplication and that had no effect on the final results (the .top file was identical).

        Note that we are also testing the .frcmod here with the specific user FF.
        """
        ambertools_tleap = self.config.ambertools_tleap
        ambertools_script_dir = self.config.ambertools_script_dir
        if self.config.protein is None:
            protein_ff = None
        else:
            protein_ff = self.config.protein_ff

        ligand_ff = self.config.ligand_ff

        frcmod_info1 = ties.helpers.parse_frcmod_sections(self.ligA.frcmod)
        frcmod_info2 = ties.helpers.parse_frcmod_sections(self.ligZ.frcmod)

        # prep directory
        cwd = self.FRCMOD_DIR
        if not cwd.is_dir():
            cwd.mkdir(exist_ok=True)

        # fixme: use the provided cwd here, otherwise this will not work if the wrong cwd is used
        # have some conf module instead of this
        morph_frcmod = cwd / f'{self.ligA.internal_name}_{self.ligZ.internal_name}_morph.frcmod'
        with open(morph_frcmod, 'w') as FOUT:
            FOUT.write('merged frcmod\n')

            for section in ['MASS', 'BOND', 'ANGLE',
                            'DIHE', 'IMPROPER', 'NONBON']:
                section_lines = frcmod_info1[section] + frcmod_info2[section]
                FOUT.write('{0:s}\n'.format(section))
                for line in section_lines:
                    FOUT.write('{0:s}'.format(line))
                FOUT.write('\n')

            FOUT.write('\n\n')

        # this is our current frcmod file
        self.frcmod = morph_frcmod

        # as part of the .frcmod writing
        # insert dummy angles/dihedrals if a morph .frcmod requires
        # new terms between the appearing/disappearing atoms
        self._check_hybrid_frcmod(ambertools_tleap, ambertools_script_dir, protein_ff, ligand_ff)

    def _check_hybrid_frcmod(self, ambertools_tleap, ambertools_script_dir, protein_ff, ligand_ff):
        """
        Previous code: https://github.com/UCL-CCS/BacScratch/blob/master/agastya/ties_hybrid_topology_creator/output.py
        Check that the output library can be used to create a valid amber topology.
        Add missing terms with no force to pass the topology creation.
        Returns the corrected .frcmod content, otherwise throws an exception.
        """
        # prepare the working directory
        cwd = self.config.workdir / self.FRCMOD_TEST_DIR / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        if protein_ff is None:
            protein_ff = '# no protein ff needed'
        else:
            protein_ff = 'source ' + protein_ff

        # prepare tleap input
        leap_in_test = 'leap_test_morph.in'
        leap_in_conf = open(ambertools_script_dir / leap_in_test).read()
        open(cwd / leap_in_test, 'w').write(leap_in_conf.format(
                                mol2=os.path.relpath(self.suptop.mol2, cwd),
                                frcmod=os.path.relpath(self.frcmod, cwd),
                                protein_ff=protein_ff, ligand_ff=ligand_ff))

        # attempt generating the .top
        print('Create amber7 topology .top')
        try:
            tleap_process = subprocess.run([ambertools_tleap, '-s', '-f', leap_in_test],
                                           cwd=cwd, text=True, timeout=20,
                                           capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            raise Exception(
                f'ERROR: Testing the topology with tleap broke. Return code: {err.returncode} '
                f'ERROR: Ambertools output: {err.stdout}') from err

        # save stdout and stderr
        open(cwd / 'tleap_scan_check.log', 'w').write(tleap_process.stdout + tleap_process.stderr)

        if 'Errors = 0' in tleap_process.stdout:
            print('Test hybrid .frcmod: OK, no dummy angle/dihedrals inserted.')
            return

        # extract the missing angles/dihedrals
        missing_bonds = set()
        missing_angles = []
        missing_dihedrals = []
        for line in tleap_process.stdout.splitlines():
            if "Could not find bond parameter for:" in line:
                bond = line.split(':')[-1].strip()
                missing_bonds.add(bond)
            elif "Could not find angle parameter:" in line:
                cols = line.split(':')
                angle = cols[-1].strip()
                if angle not in missing_angles:
                    missing_angles.append(angle)
            elif "No torsion terms for" in line:
                cols = line.split()
                torsion = cols[-1].strip()
                if torsion not in missing_dihedrals:
                    missing_dihedrals.append(torsion)

        modified_hybrid_frcmod = cwd / f'{self.internal_name}_corrected.frcmod'
        if missing_angles or missing_dihedrals:
            print('WARNING: Adding dummy dihedrals/angles to frcmod to generate .top')
            # read the original frcmod
            frcmod_lines = open(self.frcmod).readlines()
            # overwriting the .frcmod with dummy angles/dihedrals
            with open(modified_hybrid_frcmod, 'w') as NEW_FRCMOD:
                for line in frcmod_lines:
                    NEW_FRCMOD.write(line)
                    if 'BOND' in line:
                        for bond  in missing_bonds:
                            dummy_bond = f'{bond:<14}0  180  \t\t# Dummy bond\n'
                            NEW_FRCMOD.write(dummy_bond)
                            print(f'Added dummy bond: "{dummy_bond}"')
                    if 'ANGLE' in line:
                        for angle in missing_angles:
                            dummy_angle = f'{angle:<14}0  120.010  \t\t# Dummy angle\n'
                            NEW_FRCMOD.write(dummy_angle)
                            print(f'Added dummy angle: "{dummy_angle}"')
                    if 'DIHE' in line:
                        for dihedral in missing_dihedrals:
                            dummy_dihedral = f'{dihedral:<14}1  0.00  180.000  2.000   \t\t# Dummy dihedrals\n'
                            NEW_FRCMOD.write(dummy_dihedral)
                            print(f'Added dummy dihedral: "{dummy_dihedral}"')

            # update our tleap input test to use the corrected file
            leap_in_test_corrected = cwd / 'leap_test_morph_corrected.in'
            open(leap_in_test_corrected, 'w').write(leap_in_conf.format(
                                mol2=os.path.relpath(self.suptop.mol2, cwd),
                                frcmod=os.path.relpath(modified_hybrid_frcmod, cwd),
                                protein_ff=protein_ff, ligand_ff=ligand_ff))

            # verify that adding the dummy angles/dihedrals worked
            tleap_process = subprocess.run([ambertools_tleap, '-s', '-f', leap_in_test_corrected],
                                           cwd=cwd, text=True, timeout=60 * 10, capture_output=True, check=True)

            if not "Errors = 0" in tleap_process.stdout:
                raise Exception('ERROR: Could not generate the .top file after adding dummy angles/dihedrals')


        print('Morph .frcmod after the insertion of dummy angle/dihedrals: OK')
        # set this .frcmod as the correct one now,
        self.frcmod_before_correction = self.frcmod
        self.frcmod = modified_hybrid_frcmod

    def overlap_fractions(self):
        matched_fraction_left = len(self.suptop.matched_pairs) / float(len(self.suptop.top1))
        matched_fraction_right = len(self.suptop.matched_pairs) / float(len(self.suptop.top2))
        disappearing_atoms_fraction = (len(self.suptop.top1) - len(self.suptop.matched_pairs)) \
                                   / float(len(self.suptop.top1)) * 100
        appearing_atoms_fraction = (len(self.suptop.top2) - len(self.suptop.matched_pairs)) \
                                   / float(len(self.suptop.top2)) * 100

        return matched_fraction_left, matched_fraction_right, disappearing_atoms_fraction, appearing_atoms_fraction
