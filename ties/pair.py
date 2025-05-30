import os
import json
import copy
import shutil
import subprocess
import logging
from pathlib import Path

import parmed

import ties.helpers
from ties.generator import get_atoms_bonds_from_mol2
from ties.topology_superimposer import superimpose_topologies
import ties.config
import ties.ligand


logger = logging.getLogger(__name__)


class Pair():
    """
    Facilitates the creation of morphs.
    It offers functionality related to a pair of ligands (a transformation).

    :param ligA: The ligand to be used as the starting state for the transformation.
    :type ligA: :class:`Ligand` or string
    :param ligZ: The ligand to be used as the ending point of the transformation.
    :type ligZ: :class:`Ligand` or string
    :param config: The configuration object holding all settings.
    :type config: :class:`Config`
    """

    def __init__(self, ligA, ligZ, config=None, **kwargs):
        """
        Please use the Config class for the documentation of the possible kwargs.
        Each kwarg is passed to the config class.

        fixme - list all relevant kwargs here

            param ligand_net_charge: integer, net charge of each ligand (has to be the same)
        """

        # create a new config if it is not provided
        self.config = ties.config.Config() if config is None else config

        # channel all config variables to the config class
        self.config.set_configs(**kwargs)

        # tell Config about the ligands if necessary
        if self.config.ligands is None:
            self.config.ligands = [ligA, ligZ]

        # create ligands if they're just paths
        if isinstance(ligA, ties.ligand.Ligand):
            self.ligA = ligA
        else:
            self.ligA = ties.ligand.Ligand(ligA, self.config)

        if isinstance(ligZ, ties.ligand.Ligand):
            self.ligZ = ligZ
        else:
            self.ligZ = ties.ligand.Ligand(ligZ, self.config)

        # initialise the handles to the molecules that morph
        self.current_ligA = self.ligA.current
        self.current_ligZ = self.ligZ.current

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
        """
        Please see :class:`Config` class for the documentation of kwargs. The passed kwargs overwrite the config
        object passed in the constructor.

        fixme - list all relevant kwargs here

        :param use_element_in_superimposition: bool whether the superimposition should rely on the element initially,
            before refining the results with a more specific check of the atom type.
        :param manually_matched_atom_pairs:
        :param manually_mismatched_pairs:
        :param redistribute_q_over_unmatched:
        """
        self.config.set_configs(**kwargs)

        # use ParmEd to load the files
        # fixme - move this to the Morph class instead of this place,
        # fixme - should not squash all messsages. For example, wrong type file should not be squashed
        leftlig_atoms, leftlig_bonds, rightlig_atoms, rightlig_bonds, parmed_ligA, parmed_ligZ = \
            get_atoms_bonds_from_mol2(self.current_ligA, self.current_ligZ,
                                      use_general_type=self.config.use_element_in_superimposition)
        # fixme - manual match should be improved here and allow for a sensible format.

        # in case the atoms were renamed, pass the names via the map renaming map
        # TODO
        # ligZ_old_new_atomname_map
        new_mismatch_names = []
        for a, z in self.config.manually_mismatched_pairs:
            new_names = (self.ligA.rev_renaming_map[a], self.ligZ.rev_renaming_map[z])
            logger.debug(f'Selecting mismatching atoms. The mismatch {(a, z)}) was renamed to {new_names}')
            new_mismatch_names.append(new_names)

        # assign
        # fixme - Ideally I would reuse the ParmEd data for this,
        # ParmEd can use bonds if they are present - fixme
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
            logger.debug(f'Starting nodes will be used: {starting_node_pairs}')

        logging_key = str(self)

        # fixme - simplify to only take the ParmEd as input
        suptop = superimpose_topologies(ligand1_nodes.values(), ligand2_nodes.values(),
                                        disjoint_components=self.config.allow_disjoint_components,
                                        net_charge_filter=True,
                                        pair_charge_atol=self.config.atom_pair_q_atol,
                                        net_charge_threshold=self.config.net_charge_threshold,
                                        redistribute_charges_over_unmatched=self.config.redistribute_q_over_unmatched,
                                        ignore_charges_completely=self.config.ignore_charges_completely,
                                        ignore_bond_types=True,
                                        ignore_coords=False,
                                        align_molecules=self.config.align_molecules_using_mcs,
                                        use_general_type=self.config.use_element_in_superimposition,
                                        # fixme - not the same ... use_element_in_superimposition,
                                        use_only_element=False,
                                        check_atom_names_unique=True,  # fixme - remove?
                                        starting_pairs_heuristics=self.config.superimposition_starting_heuristic,  # fixme - add to config
                                        force_mismatch=new_mismatch_names,
                                        starting_node_pairs=starting_node_pairs,
                                        parmed_ligA=parmed_ligA, parmed_ligZ=parmed_ligZ,
                                        starting_pair_seed=self.config.superimposition_starting_pairs,
                                        logging_key=logging_key,
                                        config=self.config)

        self.set_suptop(suptop, parmed_ligA, parmed_ligZ)
        # attach the used config to the suptop

        if suptop is not None:
            suptop.config = self.config
            # attach the morph to the suptop
            suptop.morph = self

        return suptop

    def set_suptop(self, suptop, parmed_ligA, parmed_ligZ):
        """
        Attach a SuperimposedTopology object along with the ParmEd objects for the ligA and ligZ.

        :param suptop: :class:`SuperimposedTopology`
        :param parmed_ligA: An ParmEd for the ligA
        :param parmed_ligZ: An ParmEd for the ligZ
        """
        self.suptop = suptop
        self.parmed_ligA = parmed_ligA
        self.parmed_ligZ = parmed_ligZ

    def make_atom_names_unique(self, out_ligA_filename=None, out_ligZ_filename=None, save=True):
        """
        Ensure that each that atoms across the two ligands have unique names.

        While renaming atoms, start with the element (C, N, ..) followed by
         the count so far (e.g. C1, C2, N1).

        Resnames are set to "INI" and "FIN", this is useful for the hybrid dual topology.

        :param out_ligA_filename: The new filenames for the ligands with renamed atoms. If None, the default
            naming convention is used.
        :type out_ligA_filename: string or bool
        :param out_ligZ_filename: The new filenames for the ligands with renamed atoms. If None, the default
            naming convention is used.
        :type out_ligZ_filename: string or bool
        :param save: Whether to save to the disk the ligands after renaming the atoms
        :type save: bool
        """

        # The A ligand is a template for the renaming
        self.ligA.correct_atom_names()

        # load both ligands
        left = parmed.load_file(str(self.ligA.current), structure=True)
        right = parmed.load_file(str(self.ligZ.current), structure=True)

        common_atom_names = {a.name for a in right.atoms}.intersection({a.name for a in left.atoms})
        atom_names_overlap = len(common_atom_names) > 0

        if atom_names_overlap or not self.ligZ.are_atom_names_correct():
            logger.debug(f'Renaming ({self.ligA.internal_name}) molecule ({self.ligZ.internal_name}) atom names are either reused or do not follow the correct format. ')
            if atom_names_overlap:
                logger.debug(f'Common atom names: {common_atom_names}')
            name_counter_L_nodes = ties.helpers.get_atom_names_counter(left.atoms)
            _, renaming_map = ties.helpers.get_new_atom_names(right.atoms, name_counter=name_counter_L_nodes)
            self.ligZ.renaming_map = renaming_map

        # rename the residue names to INI and FIN
        for atom in left.atoms:
            atom.residue = 'INI'
        for atom in right.atoms:
            atom.residue = 'FIN'

        # fixme - instead of using the save parameter, have a method pair.save(filename1, filename2) and
        #  call it when necessary.
        # prepare the destination directory
        if not save:
            return

        if out_ligA_filename is None:
            cwd = self.config.pair_unique_atom_names_dir / f'{self.ligA.internal_name}_{self.ligZ.internal_name}'
            cwd.mkdir(parents=True, exist_ok=True)

            self.current_ligA = cwd / (self.ligA.internal_name + '.mol2')
            self.current_ligZ = cwd / (self.ligZ.internal_name + '.mol2')
        else:
            self.current_ligA = out_ligA_filename
            self.current_ligZ = out_ligZ_filename

        # save the updated atom names
        left.save(str(self.current_ligA))
        right.save(str(self.current_ligZ))

    def check_json_file(self):
        """
        Performance optimisation in case TIES is rerun again. Return the first matched atoms which
        can be used as a seed for the superimposition.

        :return: If the superimposition was computed before, and the .json file is available,
            gets one of the matched atoms.
        :rtype: [(ligA_atom, ligZ_atom)]
        """
        matching_json = self.config.workdir / f'fep_{self.ligA.internal_name}_{self.ligZ.internal_name}.json'
        if not matching_json.is_file():
            return None

        return [list(json.load(matching_json.open())['matched'].items())[0]]

    def merge_frcmod_files(self, ligcom=None):
        """
        Merges the .frcmod files generated for each ligand separately, simply by adding them together.

        The duplication has no effect on the final generated topology parm7 top file.

        We are also testing the .frcmod here with the user's force field in order to check if
        the merge works correctly.

        :param ligcom: Either "lig" if only ligands are present, or "com" if the complex is present.
            Helps with the directory structure.
        :type ligcom: string "lig" or "com"
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

        cwd = self.config.workdir

        # fixme: use the provided cwd here, otherwise this will not work if the wrong cwd is used
        # have some conf module instead of this
        if ligcom:
            morph_frcmod = cwd / f'ties-{self.ligA.internal_name}-{self.ligZ.internal_name}' / ligcom / 'build' / 'hybrid.frcmod'
        else:
            # fixme - clean up
            morph_frcmod = cwd / f'ties-{self.ligA.internal_name}-{self.ligZ.internal_name}' / 'build' / 'hybrid.frcmod'
        morph_frcmod.parent.mkdir(parents=True, exist_ok=True)
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
        # this is a trick to make sure tleap has everything it needs to generate the .top file
        correction_introduced = self._check_hybrid_frcmod(ambertools_tleap, ambertools_script_dir, protein_ff, ligand_ff)
        if correction_introduced:
            # move the .frcmod which turned out to be insufficient according to the test
            shutil.move(morph_frcmod, str(self.frcmod) + '.uncorrected' )
            # now copy in place the corrected version
            shutil.copy(self.frcmod, morph_frcmod)

    def _check_hybrid_frcmod(self, ambertools_tleap, ambertools_script_dir, protein_ff, ligand_ff):
        """
        Check that the output library can be used to create a valid amber topology.
        Add missing terms with no force to pass the topology creation.
        Returns the corrected .frcmod content, otherwise throws an exception.
        """
        # prepare the working directory
        cwd = self.config.pair_morphfrmocs_tests_dir / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        if protein_ff is None:
            protein_ff = '# no protein ff needed'
        else:
            protein_ff = 'source ' + protein_ff

        # prepare the superimposed .mol2 file if needed
        if not hasattr(self.suptop, 'mol2'):
            self.suptop.write_mol2()

        # prepare tleap input
        leap_in_test = 'leap_test_morph.in'
        leap_in_conf = open(ambertools_script_dir / leap_in_test).read()
        open(cwd / leap_in_test, 'w').write(leap_in_conf.format(
                                mol2=os.path.relpath(self.suptop.mol2, cwd),
                                frcmod=os.path.relpath(self.frcmod, cwd),
                                protein_ff=protein_ff, ligand_ff=ligand_ff))

        # attempt generating the .top
        logger.debug('Create amber7 topology .top')
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
            logger.debug('Test hybrid .frcmod: OK, no dummy angle/dihedrals inserted.')
            return False

        # extract the missing angles/dihedrals
        missing_bonds = set()
        missing_angles = []
        missing_dihedrals = []
        for line in tleap_process.stdout.splitlines():
            if "Could not find bond parameter for:" in line:
                bond = line.split(':')[-1].strip()
                missing_bonds.add(bond)
            elif "Could not find angle parameter:" in line or \
                    "Could not find angle parameter for atom types:" in line:
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
            logger.debug('Adding dummy bonds+angles+dihedrals to frcmod to generate .top')
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
                            logger.debug(f'Added dummy bond: "{dummy_bond}"')
                    if 'ANGLE' in line:
                        for angle in missing_angles:
                            dummy_angle = f'{angle:<14}0  120.010  \t\t# Dummy angle\n'
                            NEW_FRCMOD.write(dummy_angle)
                            logger.debug(f'Added dummy angle: "{dummy_angle}"')
                    if 'DIHE' in line:
                        for dihedral in missing_dihedrals:
                            dummy_dihedral = f'{dihedral:<14}1  0.00  180.000  2.000   \t\t# Dummy dihedrals\n'
                            NEW_FRCMOD.write(dummy_dihedral)
                            logger.debug(f'Added dummy dihedral: "{dummy_dihedral}"')

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


        logger.debug('Morph .frcmod after the insertion of dummy angle/dihedrals: OK')
        # set this .frcmod as the correct one now,
        self.frcmod_before_correction = self.frcmod
        self.frcmod = modified_hybrid_frcmod
        return True

    def overlap_fractions(self):
        """
        Calculate the size of the common area.

        :return: Four decimals capturing: 1) the fraction of the common size with respect to the ligA topology,
            2) the fraction of the common size with respect to the ligZ topology,
            3) the percentage of the disappearing atoms in the disappearing molecule
            4) the percentage of the appearing atoms  in the appearing molecule
        :rtype: [float, float, float, float]
        """

        if self.suptop is None:
            return 0, 0, float('inf'), float('inf')
        else:
            mcs_size = len(self.suptop.matched_pairs)

        matched_fraction_left = mcs_size / float(len(self.suptop.top1))
        matched_fraction_right = mcs_size / float(len(self.suptop.top2))
        disappearing_atoms_fraction = (len(self.suptop.top1) - mcs_size) \
                                   / float(len(self.suptop.top1)) * 100
        appearing_atoms_fraction = (len(self.suptop.top2) - mcs_size) \
                                   / float(len(self.suptop.top2)) * 100

        return matched_fraction_left, matched_fraction_right, disappearing_atoms_fraction, appearing_atoms_fraction
