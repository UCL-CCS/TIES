import os
import json
import copy
import subprocess
from pathlib import Path

from ties.helpers import load_MDAnalysis_atom_group, get_atom_names_counter, rename_ligand, parse_frcmod_sections
from ties.generator import get_atoms_bonds_from_mol2
from ties.topology_superimposer import superimpose_topologies


class Morph():
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

    def __init__(self, ligA, ligZ, config):
        self.ligA = ligA
        self.ligZ = ligZ
        self.config = config
        # fixme rm?
        self.workplace_root = config.workdir
        self.UNIQUE_ATOM_NAMES = self.workplace_root / Morph.UNIQUE_ATOM_NAMES_name
        self.FRCMOD_DIR = self.workplace_root / Morph.FRCMOD_DIR_name
        self.FRCMOD_TEST_DIR = self.workplace_root / Morph.FRCMOD_TEST_DIR_name

        self.internal_name = f'{ligA.internal_name}_{ligZ.internal_name}'
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

    def get_suptop(self, morph, manual_match=None, force_mismatch=None,
                   disjoint_components=False, net_charge_filter=True,
                   ignore_charges_completely=False,
                   use_general_type=True,
                   ignore_bond_types=True,
                   ignore_coords=False,
                   left_coords_are_ref=True,
                   align_molecules=True,
                   use_only_gentype=False,
                   check_atom_names_unique=True,
                   pair_charge_atol=0.1,
                   net_charge_threshold=0.1,
                   redistribute_charges_over_unmatched=True,
                   starting_pairs_heuristics=True):
        # use MDAnalysis to load the files
        # fixme - move this to the Morph class instead of this place,
        # fixme - should not squash all messsages. For example, wrong type file should not be squashed
        leftlig_atoms, leftlig_bonds, rightlig_atoms, rightlig_bonds, mda_l1, mda_l2 = \
            get_atoms_bonds_from_mol2(morph.renamed_ligA, morph.renamed_ligZ, use_general_type=use_general_type)
        # fixme - manual match should be improved here and allow for a sensible format.

        # empty lists count as None
        if not force_mismatch:
            force_mismatch = None

        # assign
        # fixme - Ideally I would reuse the mdanalysis data for this,
        # MDAnalysis can use bonds if they are present - fixme
        # map atom IDs to their objects
        ligand1_nodes = {}
        for atomNode in leftlig_atoms:
            ligand1_nodes[atomNode.get_id()] = atomNode
        # link them together
        for nfrom, nto, btype in leftlig_bonds:
            ligand1_nodes[nfrom].bind_to(ligand1_nodes[nto], btype)

        ligand2_nodes = {}
        for atomNode in rightlig_atoms:
            ligand2_nodes[atomNode.get_id()] = atomNode
        for nfrom, nto, btype in rightlig_bonds:
            ligand2_nodes[nfrom].bind_to(ligand2_nodes[nto], btype)

        # fixme - this should be moved out of here,
        #  ideally there would be a function in the main interface for this
        manual_match = [] if manual_match is None else manual_match
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
                                         starting_node_pairs=starting_node_pairs,
                                         force_mismatch=force_mismatch,
                                         disjoint_components=disjoint_components,
                                         net_charge_filter=net_charge_filter,
                                         pair_charge_atol=pair_charge_atol,
                                         ligand_l_mda=mda_l1, ligand_r_mda=mda_l2,
                                         ignore_charges_completely=ignore_charges_completely,
                                         ignore_bond_types=ignore_bond_types,
                                         ignore_coords=ignore_coords,
                                         left_coords_are_ref=left_coords_are_ref,
                                         align_molecules=align_molecules,
                                         use_general_type=use_general_type,
                                         use_only_element=use_only_gentype,
                                         check_atom_names_unique=check_atom_names_unique,
                                         net_charge_threshold=net_charge_threshold,
                                         redistribute_charges_over_unmatched=redistribute_charges_over_unmatched,
                                         starting_pairs_heuristics=starting_pairs_heuristics)

        assert len(suptops) == 1
        morph.set_suptop(suptops[0], mda_l1, mda_l2)

    def set_suptop(self, suptop, mda_l1, mda_l2):
        self.suptop = suptop
        self.mda_l1 = mda_l1
        self.mda_l2 = mda_l2

    def unique_atomres_names(self):
        """
        Ensure that each has a name that is unique to both ligands.

        Rename the ligand to ensure that no atom has the same name
        name atoms using the first letter (C, N, ..) and count them
        keep the names if possible (ie if they are already different)

        Resnames are set to "INI" and "FIN", this is useful for the hybrid dual topology

        fixme - rewrite this to be LigA.unique_names(LigZ) so that they can do it internally themselves?
        """
        # load both ligands
        left = load_MDAnalysis_atom_group(self.ligA.current)
        right = load_MDAnalysis_atom_group(self.ligZ.current)

        ligands_atom_names_overlap = len(set(right.atoms.names).intersection(set(left.atoms.names))) > 0

        if ligands_atom_names_overlap:
            print(f'Renaming right molecule ({self.ligZ.internal_name}) atom names because they are not unique')
            name_counter_L_nodes = get_atom_names_counter(left.atoms)
            rename_ligand(right.atoms, name_counter=name_counter_L_nodes)

        # rename the residue names to INI and FIN
        left.residues.resnames = ['INI']
        right.residues.resnames = ['FIN']

        # prepare the destination directory
        cwd = self.workplace_root / self.UNIQUE_ATOM_NAMES / f'{self.ligA.internal_name}_{self.ligZ.internal_name}'
        cwd.mkdir(parents=True, exist_ok=True)

        # save the updated atom names
        self.renamed_ligA = cwd / (self.ligA.internal_name + '.mol2')
        left.atoms.write(self.renamed_ligA)
        self.renamed_ligZ = cwd / (self.ligZ.internal_name + '.mol2')
        right.atoms.write(self.renamed_ligZ)

        # update the Ligands to use the renamed version

    def check_json_file(self):
        # An optimisation for testing
        # Check if the json was produced before, use it to
        # find the starting point to speed up the testing
        matching_json = self.workplace_root / f'fep_{self.ligA.internal_name}_{self.ligZ.internal_name}.json'
        if not matching_json.is_file():
            return None

        return [list(json.load(matching_json.open())['matched'].items())[0]]

    def write_summary_json(self):
        """
        Writes a simple .json file with a summary of which atoms are classified as appearing, disappearing.
        """

        # store at the root for now
        matching_json = self.workplace_root / f'fep_{self.ligA.internal_name}_{self.ligZ.internal_name}.json'

        with open(matching_json, 'w') as FOUT:
            # use json format, only use atomNames
            app, dis = self.suptop.get_single_topology_app()
            summary = {
                # the dual topology information
                'matched': {str(n1): str(n2) for n1, n2 in self.suptop.matched_pairs},
                'appearing': list(map(str, self.suptop.get_appearing_atoms())),
                'disappearing': list(map(str, self.suptop.get_disappearing_atoms())),
                # single topology information
                'single_top_matched': {str(n1): str(n2) for n1, n2 in self.suptop.get_single_topology_region()},
                # NAMD hybrid single-dual topology info
                'single_top_appearing': list(map(str, app)),
                'single_top_disappearing': list(map(str, dis)),
            }
            FOUT.write(json.dumps(summary, indent=4))

        self.summary = summary

    def write_pdb(self, hybrid_single_dual_top):
        morph_pdb_path = self.workplace_root / f'{self.ligA.internal_name}_{self.ligZ.internal_name}_morph.pdb'

        # def write_morph_top_pdb(filepath, mda_l1, mda_l2, suptop, hybrid_single_dual_top=False):
        if hybrid_single_dual_top:
            # the NAMD hybrid single dual topology
            # rename the ligand on the left to INI
            # and the ligand on the right to END

            # make a copy of the suptop here to ensure that the modifications won't affect it
            st = copy.copy(self.suptop)

            # first, set all the matched pairs to -2 and 2 (single topology)
            # regardless of how they were mismatched
            raise NotImplementedError('Cannot yet write hybrid single dual topology .pdb file')

            # then, set the different atoms to -1 and 1 (dual topology)

            # save in a single PDB file
            # Note that the atoms from left to right
            # in the single topology region have to
            # be separated by the same number
            # fixme - make a check for that
            return
        # fixme - find another library that can handle writing to a PDB file, MDAnalysis
        # save the ligand with all the appropriate atomic positions, write it using the pdb format
        # pdb file format: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        # write a dual .pdb file
        with open(morph_pdb_path, 'w') as FOUT:
            for atom in self.mda_l1.atoms:
                """
                There is only one forcefield which is shared across the two topologies. 
                Basically, we need to check whether the atom is in both topologies. 
                If that is the case, then the atom should have the same name, and therefore appear only once. 
                However, if there is a new atom, it should be specfically be outlined 
                that it is 1) new and 2) the right type
                """
                # write all the atoms if they are matched, that's the common part
                REMAINS = 0
                if self.suptop.contains_left_atom_name(atom.name.upper()):
                    line = f"ATOM  {atom.id:>5d} {atom.name:>4s} {atom.resname:>3s}  " \
                           f"{atom.resid:>4d}    " \
                           f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}" \
                           f"{1.0:>6.2f}{REMAINS:>6.2f}" + (' ' * 11) + \
                           '  ' + '  ' + '\n'
                    FOUT.write(line)
                else:
                    # this atom was not found, this means it disappears, so we should update the
                    DISAPPEARING_ATOM = -1.0
                    line = f"ATOM  {atom.id:>5d} {atom.name:>4s} {atom.resname:>3s}  " \
                           f"{atom.resid:>4d}    " \
                           f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}" \
                           f"{1.0:>6.2f}{DISAPPEARING_ATOM:>6.2f}" + \
                           (' ' * 11) + \
                           '  ' + '  ' + '\n'
                    FOUT.write(line)
            # add atoms from the right topology,
            # which are going to be created
            for atom in self.mda_l2.atoms:
                if not self.suptop.contains_right_atom_name(atom.name.upper()):
                    APPEARING_ATOM = 1.0
                    line = f"ATOM  {atom.id:>5d} {atom.name:>4s} {atom.resname:>3s}  " \
                           f"{atom.resid:>4d}    " \
                           f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}" \
                           f"{1.0:>6.2f}{APPEARING_ATOM:>6.2f}" + \
                           (' ' * 11) + \
                           '  ' + '  ' + '\n'
                    FOUT.write(line)

    def write_hybrid_mol2(self, use_left_charges=True, use_left_coords=True):
        hybrid_mol2 = self.workplace_root / f'{self.ligA.internal_name}_{self.ligZ.internal_name}_morph.mol2'

        # fixme - make this as a method of suptop as well
        # recreate the mol2 file that is merged and contains the correct atoms from both
        # mol2 format: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
        # fixme - build this molecule using the MDAnalysis builder instead of the current approach
        # however, MDAnalysis currently cannot convert pdb into mol2? ...
        # where the formatting is done manually
        with open(hybrid_mol2, 'w') as FOUT:
            bonds = self.suptop.get_dual_topology_bonds()

            FOUT.write('@<TRIPOS>MOLECULE ' + os.linesep)
            # name of the molecule
            FOUT.write('merged ' + os.linesep)
            # num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
            # fixme this is tricky
            FOUT.write(f'{self.suptop.get_unique_atom_count():d} '
                       f'{len(bonds):d}' + os.linesep)
            # mole type
            FOUT.write('SMALL ' + os.linesep)
            # charge_type
            FOUT.write('NO_CHARGES ' + os.linesep)
            FOUT.write(os.linesep)

            # write the atoms
            FOUT.write('@<TRIPOS>ATOM ' + os.linesep)
            # atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
            # e.g.
            #       1 O4           3.6010   -50.1310     7.2170 o          1 L39      -0.815300

            # so from the two topologies all the atoms are needed and they need to have a different atom_id
            # so we might need to name the atom_id for them, other details are however pretty much the same
            # the importance of atom_name is difficult to estimate

            # we are going to assign IDs in the superimposed topology in order to track which atoms have IDs
            # and which don't

            # fixme - for writing, modify things to achieve the desired output
            # note - we are modifying in place our atoms
            for left, right in self.suptop.matched_pairs:
                print(
                    f'Aligned {left.originalAtomName} id {left.atomId} with {right.originalAtomName} id {right.atomId}')
                if not use_left_charges:
                    left.charge = right.charge
                if not use_left_coords:
                    left.position = right.position

            subst_id = 1  # resid basically
            # write all the atoms that were matched first with their IDs
            # prepare all the atoms, note that we use primarily the left ligand naming
            all_atoms = [left for left, right in self.suptop.matched_pairs] + self.suptop.get_unmatched_atoms()
            unmatched_atoms = self.suptop.get_unmatched_atoms()
            # reorder the list according to the ID
            all_atoms.sort(key=lambda atom: self.suptop.get_generated_atom_id(atom))

            for atom in all_atoms:
                FOUT.write(f'{self.suptop.get_generated_atom_id(atom)} {atom.name} '
                           f'{atom.position[0]:.4f} {atom.position[1]:.4f} {atom.position[2]:.4f} '
                           f'{atom.type.lower()} {subst_id} {atom.resname} {atom.charge:.6f} {os.linesep}')

            # for left_atom, _ in suptop.matched_pairs:
            #     # note that the atom id is the most important
            #     FOUT.write(f'{suptop.get_generated_atom_ID(left_atom)} {left_atom.atomName} '
            #                f'{left_atom.position[0]:.4f} {left_atom.position[1]:.4f} {left_atom.position[2]:.4f} '
            #                f'{left_atom.type} {subst_id} {left_atom.resname} {left_atom.charge} {os.linesep}')

            # write the IDs for the atoms which are appearing/disappearing
            # for unmatched in suptop.get_unmatched_atoms():
            #     FOUT.write(f'{suptop.get_generated_atom_ID(unmatched)} {unmatched.atomName} '
            #                f'{unmatched.position[0]:.4f} {unmatched.position[1]:.4f} {unmatched.position[2]:.4f} '
            #                f'{unmatched.type} {subst_id} {unmatched.resname} {unmatched.charge} {os.linesep}')

            FOUT.write(os.linesep)

            # write bonds
            FOUT.write('@<TRIPOS>BOND ' + os.linesep)

            # we have to list every bond:
            # 1) all the bonds between the paired atoms, so that part is easy
            # 2) bonds which link the disappearing atoms, and their connection to the paired atoms
            # 3) bonds which link the appearing atoms, and their connections to the paired atoms

            bond_counter = 1
            list(bonds)
            for bond_from_id, bond_to_id, bond_type in sorted(list(bonds), key=lambda b: b[:2]):
                # Bond Line Format:
                # bond_id origin_atom_id target_atom_id bond_type [status_bits]
                FOUT.write(f'{bond_counter} {bond_from_id} {bond_to_id} {bond_type}' + os.linesep)
                bond_counter += 1

        self.mol2 = hybrid_mol2

    def join_frcmod_files(self, ambertools_tleap, ambertools_script_dir, protein_ff, ligand_ff):
        """
        Copied from the previous TIES. It's simpler and this appears to work fine.
        I tested the duplication and that seemed to have no effect on the final results.

        Note that we are also testing the .frcmod here with the specific user FF.
        """
        frcmod_info1 = parse_frcmod_sections(self.ligA.frcmod)
        frcmod_info2 = parse_frcmod_sections(self.ligZ.frcmod)

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
        cwd = self.workplace_root / self.FRCMOD_TEST_DIR / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        # prepare tleap input
        leap_in_test = 'leap_test_morph.in'
        leap_in_conf = open(ambertools_script_dir / leap_in_test).read()
        open(cwd / leap_in_test, 'w').write(leap_in_conf.format(
                                mol2=os.path.relpath(self.mol2, cwd),
                                frcmod=os.path.relpath(self.frcmod, cwd),
                                protein_ff=protein_ff, ligand_ff=ligand_ff))

        # attempt generating the .top
        print('Create amber7 topology .top')
        tleap_process = subprocess.run([ambertools_tleap, '-s', '-f', leap_in_test],
                                       cwd=cwd, text=True, timeout=20,
                                       capture_output=True, check=True)

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
                                mol2=os.path.relpath(self.mol2, cwd),
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
