import os
import json
import copy
import shutil
import subprocess
import math
from collections import OrderedDict
from pathlib import Path
import warnings

import MDAnalysis as mda
import numpy as np

from ties.topology_superimposer import get_atoms_bonds_from_mol2, \
    superimpose_topologies, element_from_type, get_atoms_bonds_from_ac


def prepare_inputs(morph,
                   dir_prefix,
                   protein=None,
                   namd_script_loc=None,
                   submit_script=None,
                   scripts_loc=None,
                   tleap_in=None,
                   protein_ff=None,
                   ligand_ff=None,
                   net_charge=None,
                   ambertools_script_dir=None,
                   ambertools_tleap=None,
                   hybrid_topology=False,
                   vmd_vis_script=None,
                   md_engine=False,
                   lambda_rep_dir_tree=False):

    cwd = dir_prefix / f'{morph.internal_name}'
    if not cwd.is_dir():
        cwd.mkdir(parents=True, exist_ok=True)

    # Agastya: simple format for appearing, disappearing atoms
    open(cwd / 'disappearing_atoms.txt', 'w').write(' '.join([a.name for a in morph.suptop.get_disappearing_atoms()]))
    open(cwd / 'appearing_atoms.txt', 'w').write(' '.join([a.name for a in morph.suptop.get_appearing_atoms()]))

    # copy the protein complex .pdb
    if protein is not None:
        shutil.copy(protein, cwd / 'protein.pdb')

    # copy the hybrid ligand (topology and .frcmod)
    shutil.copy(morph.mol2, cwd / 'morph.mol2')
    shutil.copy(morph.frcmod, cwd / 'morph.frcmod')

    # determine the number of ions to neutralise the ligand charge
    if net_charge < 0:
        Na_num = math.fabs(net_charge)
        Cl_num = 0
    elif net_charge > 0:
        Cl_num = net_charge
        Na_num = 0
    else:
        Cl_num = Na_num = 0

    # copy the protein tleap input file (ambertools)
    # set the number of ions manually
    assert Na_num == 0 or Cl_num == 0, 'At this point the number of ions should have be resolved'
    if Na_num == 0:
        tleap_Na_ions = '# no Na+ added'
    elif Na_num > 0:
        tleap_Na_ions = 'addIons sys Na+ %d' % Na_num
    if Cl_num == 0:
        tleap_Cl_ions = '# no Cl- added'
    elif Cl_num > 0:
        tleap_Cl_ions = 'addIons sys Cl- %d' % Cl_num

    leap_in_conf = open(ambertools_script_dir / tleap_in).read()
    open(cwd / 'leap.in', 'w').write(leap_in_conf.format(protein_ff=protein_ff,
                                                              ligand_ff=ligand_ff,
                                                              NaIons=tleap_Na_ions,
                                                              ClIons=tleap_Cl_ions))

    # ambertools tleap: combine ligand+complex, solvate, generate amberparm
    log_filename = cwd / 'generate_sys_top.log'
    with open(log_filename, 'w') as LOG:
        try:
            subprocess.run([ambertools_tleap, '-s', '-f', 'leap.in'],
                           stdout=LOG, stderr=LOG,
                           cwd = cwd,
                           check=True, text=True, timeout=30)
            hybrid_solv = cwd / 'sys_solv.pdb' # generated
            # check if the solvation is correct
        except subprocess.CalledProcessError as E:
            print('ERROR: occurred when trying to parse the protein.pdb with tleap. ')
            print(f'ERROR: The output was saved in the directory: {cwd}')
            print(f'ERROR: can be found in the file: {log_filename}')
            raise E

    # for hybrid single-dual topology approach, we have to use a different approach

    # generate the merged .fep file
    complex_solvated_fep = cwd / 'sys_solv_fep.pdb'
    correct_fep_tempfactor(morph.summary, hybrid_solv, complex_solvated_fep,
                           hybrid_topology)

    # fixme - check that the protein does not have the same resname?

    # calculate PBC for an octahedron
    solv_oct_boc = extract_PBC_oct_from_tleap_log(cwd / "leap.log")

    #these are the names for the source of the input, output names are fixed
    if md_engine in ('NAMD', 'namd'):
        min_script = "min.namd"
        eq_script = "eq.namd"
        prod_script = "prod.namd"

    elif md_engine in ('NAMD3', 'namd3'):
        min_script = "min3.namd"
        eq_script = "eq3.namd"
        prod_script = "prod3.namd"

    elif md_engine in ('OpenMM', 'openmm'):
        min_script = None
        eq_script = None
        prod_script = None

    else:
        raise ValueError('Unknown engine {}'.format(md_engine))

    #Make build and replica_conf dirs
    # replica_conf contains NAMD scripts
    if 'namd' in md_engine.lower():
        if not os.path.exists(cwd / 'replica_conf'):
            os.makedirs(cwd / 'replica_conf')
        else:
            shutil.rmtree(cwd / 'replica_conf')
            os.makedirs(cwd / 'replica_conf')

        # populate the replica_conf dir with scripts
        # minimization scripts
        init_namd_file_min(namd_script_loc, cwd / 'replica_conf', min_script,
                           structure_name='sys_solv', pbc_box=solv_oct_boc)
        # equilibriation scripts, note eq_namd_filenames unused
        generate_namd_eq(namd_script_loc / eq_script, cwd / 'replica_conf', structure_name='sys_solv')
        # production script
        generate_namd_prod(namd_script_loc / prod_script, cwd / 'replica_conf/sim1.conf', structure_name='sys_solv')

    elif 'openmm' in md_engine.lower():
        ties_script = open(scripts_loc / 'openmm' / 'TIES.cfg').read().format(structure_name='sys_solv',
                                                                              cons_file='cons.pdb', **solv_oct_boc)
        open(os.path.join(cwd, 'TIES.cfg'), 'w').write(ties_script)

    # build contains all input ie. positions, parameters and constraints
    if not os.path.exists(cwd / 'build'):
        os.makedirs(cwd / 'build')
    else:
        shutil.rmtree(cwd / 'build')
        os.makedirs(cwd / 'build')

    # populate the build dir with positions, parameters and constraints
    # generate 4 different constraint .pdb files (it uses B column), note constraint_files unused
    create_constraint_files(cwd / 'build' / hybrid_solv, os.path.join(cwd, 'build', 'cons.pdb'))
    #pdb, positions
    shutil.move(hybrid_solv, cwd / 'build')
    #pdb.fep, alchemical atoms
    shutil.move(complex_solvated_fep, cwd / 'build')
    #prmtop, topology
    shutil.move(cwd / "sys_solv.top", cwd / 'build')

    #copy replic-conf scripts
    rep_conf_scripts = ['eq0-replicas.conf', 'eq1-replicas.conf', 'eq2-replicas.conf', 'sim1-replicas.conf']
    for f in rep_conf_scripts:
        shutil.copy(namd_script_loc / f, cwd / 'replica_conf')

    # Generate the directory structure for all the lambdas, and copy the files
    if 'namd' in md_engine.lower():
        lambdas = [0, 0.05] + list(np.linspace(0.1, 0.9, 9)) + [0.95, 1]
    else:
        lambdas = range(13)
    if lambda_rep_dir_tree:
        for lambda_step in lambdas:
            if 'namd' in md_engine.lower():
                lambda_path = cwd / f'LAMBDA_{lambda_step:.2f}'
            else:
                lambda_path = cwd / 'LAMBDA_{}'.format(int(lambda_step))
            if not os.path.exists(lambda_path):
                os.makedirs(lambda_path)

            # for each lambda create 5 replicas
            for replica_no in range(1, 5 + 1):
                replica_dir = lambda_path / f'rep{replica_no}'
                if not os.path.exists(replica_dir):
                    os.makedirs(replica_dir)

                # set the lambda value for the directory,
                # this file is used by NAMD tcl scripts
                open(replica_dir / 'lambda', 'w').write(f'{lambda_step:.2f}')

                #create output dirs equilibriation and simulation
                if not os.path.exists(replica_dir / 'equilibration'):
                    os.makedirs(replica_dir / 'equilibration')

                if not os.path.exists(replica_dir / 'simulation'):
                    os.makedirs(replica_dir / 'simulation')

                # copy a submit script
                if submit_script is not None:
                    shutil.copy(namd_script_loc / submit_script, replica_dir / 'submit.sh')

    # copy the visualisation script as hidden
    shutil.copy(vmd_vis_script, cwd / 'vis.vmd')
    # simplify the vis.vmd use
    vis_sh = Path(cwd / 'vis.sh')
    vis_sh.write_text('#!/bin/sh \nvmd -e vis.vmd')
    vis_sh.chmod(0o755)


def _merge_frcmod_section(ref_lines, other_lines):
    """
    A helper function for merging lines in .frcmod files.
    Note that the order has to be kept. This is because some lines need to follow other lines.
    In this case, we exclude lines that are present in ref_lines.
    fixme - since there are duplicate lines, we need to check the duplicates and their presence,
    """
    merged_section = copy.copy(ref_lines)
    for line in other_lines:
        if line not in ref_lines:
            merged_section.append(line)

    return merged_section


def join_frcmod_files(f1, f2, output_filepath):
    """
    This implementation should be used. Switch to join_frcmod_files2.
    This version might be removed if the simple approach is fine.
    """
    # fixme - load f1 and f2

    def get_section(name, rlines):
        """
        Chips away from the lines until the section is ready

        fixme is there a .frcmod reader in ambertools?
        http://ambermd.org/FileFormats.php#frcmod
        """
        section_names = ['MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON']
        assert name in rlines.pop().strip()

        section = []
        while not (len(rlines) == 0 or any(rlines[-1].startswith(sname) for sname in section_names)):
            nextl = rlines.pop().strip()
            if nextl == '':
                continue
            # depending on the column name, parse differently
            if name == 'ANGLE':
                # e.g.
                # c -cc-na   86.700     123.270   same as c2-cc-na, penalty score=  2.6
                atom_types = nextl[:8]
                other = nextl[9:].split()[::-1]
                # The harmonic force constants for the angle "ITT"-"JTT"-
                #                     "KTT" in units of kcal/mol/(rad**2) (radians are the
                #                     traditional unit for angle parameters in force fields).
                harmonicForceConstant = float(other.pop())
                # TEQ        The equilibrium bond angle for the above angle in degrees.
                eq_bond_angle = float(other.pop())
                # the overall angle
                section.append([atom_types, harmonicForceConstant, eq_bond_angle])
            elif name == 'DIHE':
                # e.g.
                # ca-ca-cd-cc   1    0.505       180.000           2.000      same as c2-ce-ca-ca, penalty score=229.0
                atom_types = nextl[:11]
                other = nextl[11:].split()[::-1]
                """
                IDIVF      The factor by which the torsional barrier is divided.
                    Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
                    details. Basically, the actual torsional potential is

                           (PK/IDIVF) * (1 + cos(PN*phi - PHASE))

                 PK         The barrier height divided by a factor of 2.

                 PHASE      The phase shift angle in the torsional function.

                            The unit is degrees.

                 PN         The periodicity of the torsional barrier.
                            NOTE: If PN .lt. 0.0 then the torsional potential
                                  is assumed to have more than one term, and the
                                  values of the rest of the terms are read from the
                                  next cards until a positive PN is encountered.  The
                                  negative value of pn is used only for identifying
                                  the existence of the next term and only the
                                  absolute value of PN is kept.
                """
                IDIVF = float(other.pop())
                PK = float(other.pop())
                PHASE = float(other.pop())
                PN = float(other.pop())
                section.append([atom_types, IDIVF, PK, PHASE, PN])
            elif name == 'IMPROPER':
                # e.g.
                # cc-o -c -o          1.1          180.0         2.0          Using general improper torsional angle  X- o- c- o, penalty score=  3.0)
                # ...  IDIVF , PK , PHASE , PN
                atom_types = nextl[:11]
                other = nextl[11:].split()[::-1]
                # fixme - what is going on here? why is not generated this number?
                # IDIVF = float(other.pop())
                PK = float(other.pop())
                PHASE = float(other.pop())
                PN = float(other.pop())
                if PN < 0:
                    raise Exception('Unimplemented - ordering using with negative 0')
                section.append([atom_types, PK, PHASE, PN])
            else:
                section.append(nextl.split())
        return {name: section}

    def load_frcmod(filepath):
        # remark line
        rlines = open(filepath).readlines()[::-1]
        assert 'Remark' in rlines.pop()

        parsed = OrderedDict()
        for section_name in ['MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON']:
            parsed.update(get_section(section_name, rlines))

        return parsed

    def join_frcmod(left_frc, right_frc):
        joined = OrderedDict()
        for left, right in zip(left_frc.items(), right_frc.items()):
            lname, litems = left
            rname, ritems = right
            assert lname == rname

            joined[lname] = copy.copy(litems)

            if lname == 'MASS':
                if len(litems) > 0 or len(ritems) > 0:
                    raise Exception('Unimplemented')
            elif lname == 'BOND':
                for ritem in ritems:
                    if len(litems) > 0 or len(ritems) > 0:
                        if ritem not in joined[lname]:
                            raise Exception('Unimplemented')
            # ANGLE, e.g.
            # c -cc-na   86.700     123.270   same as c2-cc-na, penalty score=  2.6
            elif lname == 'ANGLE':
                for ritem in ritems:
                    # if the item is not in the litems, add it there
                    # extra the first three terms to determine if it is present
                    # fixme - note we are ignoring the "same as" note
                    if ritem not in joined[lname]:
                        joined[lname].append(ritem)
            elif lname == 'DIHE':
                for ritem in ritems:
                    if ritem not in joined[lname]:
                        joined[lname].append(ritem)
            elif lname == 'IMPROPER':
                for ritem in ritems:
                    if ritem not in joined[lname]:
                        joined[lname].append(ritem)
            elif lname == 'NONBON':
                # if they're empty
                if not litems and not ritems:
                    continue

                raise Exception('Unimplemented')
            else:
                raise Exception('Unimplemented')
        return joined

    def write_frcmod(frcmod, filename):
        with open(filename, 'w') as FOUT:
            FOUT.write('GENERATED .frcmod by joining two .frcmod files' + os.linesep)
            for sname, items in frcmod.items():
                FOUT.write(f'{sname}' + os.linesep)
                for item in items:
                    atom_types = item[0]
                    FOUT.write(atom_types)
                    numbers = ' \t'.join([str(n) for n in item[1:]])
                    FOUT.write(' \t' + numbers)
                    FOUT.write(os.linesep)
                # the ending line
                FOUT.write(os.linesep)

    left_frc = load_frcmod(f1)
    right_frc = load_frcmod(f2)
    joined_frc = join_frcmod(left_frc, right_frc)
    write_frcmod(joined_frc, output_filepath)


def _correct_fep_tempfactor_single_top(fep_summary, source_pdb_filename, new_pdb_filename):
    """
    Single topology version of function correct_fep_tempfactor.

    The left ligand has to be called INI
    And right FIN
    """
    u = mda.Universe(source_pdb_filename)
    if 'INI' not in u.atoms.resnames or 'FIN' not in u.atoms.resnames:
        raise Exception('Missing the resname "mer" in the pdb file prepared for fep')

    # dual-topology info
    # matched atoms are denoted -2 and 2 (morphing into each other)
    matched_disappearing = list(fep_summary['single_top_matched'].keys())
    matched_appearing = list( fep_summary['single_top_matched'].values())
    # disappearing is denoted by -1
    disappearing_atoms = fep_summary['single_top_disappearing']
    # appearing is denoted by 1
    appearing_atoms = fep_summary['single_top_appearing']

    # update the Temp column
    for atom in u.atoms:
        # ignore water and ions and non-ligand resname
        # we only modify the protein, so ignore the ligand resnames
        # fixme .. why is it called mer, is it tleap?
        if atom.resname != 'INI' and atom.resname != 'FIN':
            continue

        # if the atom was "matched", meaning present in both ligands (left and right)
        # then ignore
        # note: we only use the left ligand
        if atom.name.upper() in matched_disappearing:
            atom.tempfactor = -2
        elif atom.name.upper() in matched_appearing:
            atom.tempfactor = 2
        elif atom.name.upper() in disappearing_atoms:
            atom.tempfactor = -1
        elif atom.name.upper() in appearing_atoms:
            # appearing atoms should
            atom.tempfactor = 1
        else:
            raise Exception('This should never happen. It has to be one of the cases')

    u.atoms.write(new_pdb_filename)  # , file_format='PDB') - fixme?


def load_mda_u(filename):
    # Load a MDAnalysis file, but turn off warnings about the dual topology setting
    # fixme - fuse with the .pdb loading, no distinction needed
    warnings.filterwarnings(action='ignore', category=UserWarning,
                            message='Element information is absent or missing for a few atoms. '
                                    'Elements attributes will not be populated.'    # warning to ignore
                            )
    return mda.Universe(filename)


def correct_fep_tempfactor(fep_summary, source_pdb_filename, new_pdb_filename, hybrid_topology=False):
    """
    fixme - this function does not need to use the file?
    we have the json information available here.

    Sets the temperature column in the PDB file
    So that the number reflects the alchemical information
    Requires by NAMD in order to know which atoms
    appear (1) and which disappear (-1).
    """
    if hybrid_topology:
        # delegate correcting fep column in the pdb file
        return _correct_fep_tempfactor_single_top(fep_summary, source_pdb_filename, new_pdb_filename)

    u = load_mda_u(source_pdb_filename)
    if 'mer' not in u.atoms.resnames:
        raise Exception('Missing the resname "mer" in the pdb file prepared for fep')

    # dual-topology info
    matched = list(fep_summary['matched'].keys())
    appearing_atoms = fep_summary['appearing']
    disappearing_atoms = fep_summary['disappearing']

    # update the Temp column
    for atom in u.atoms:
        # ignore water and ions and non-ligand resname
        # we only modify the protein, so ignore the ligand resnames
        # fixme .. why is it called mer, is it tleap?
        if atom.resname != 'mer':
            continue

        # if the atom was "matched", meaning present in both ligands (left and right)
        # then ignore
        # note: we only use the left ligand
        if atom.name in matched:
            continue
        elif atom.name in appearing_atoms:
            # appearing atoms should
            atom.tempfactor = 1
        elif atom.name in disappearing_atoms:
            atom.tempfactor = -1
        else:
            raise Exception('This should never happen. It has to be one of the cases')

    u.atoms.write(new_pdb_filename)  # , file_format='PDB') - fixme?


def get_ligand_resname(filename):
    lig_resnames = mda.Universe(filename).residues.resnames
    assert len(lig_resnames) == 1
    return lig_resnames


def get_morphed_ligand_resnames(filename):
    lig_resnames = mda.Universe(filename).residues.resnames
    # assert len(lig_resnames) == 2
    return lig_resnames


def get_PBC_coords(pdb_file):
    """
    Return [x, y, z]
    """
    u = mda.Universe(pdb_file)
    x = np.abs(max(u.atoms.positions[:, 0]) - min(u.atoms.positions[:, 0]))
    y = np.abs(max(u.atoms.positions[:, 1]) - min(u.atoms.positions[:, 1]))
    z = np.abs(max(u.atoms.positions[:, 2]) - min(u.atoms.positions[:, 2]))
    return (x, y, z)


def extract_PBC_oct_from_tleap_log(leap_log):
    """
    http://ambermd.org/namd/namd_amber.html
    Return the 9 numbers for the truncated octahedron unit cell in namd
    cellBasisVector1  d         0.0            0.0
    cellBasisVector2  (-1/3)*d (2/3)sqrt(2)*d  0.0
    cellBasisVector3  (-1/3)*d (-1/3)sqrt(2)*d (-1/3)sqrt(6)*d
    """
    leapl_log_lines = open(leap_log).readlines()
    line_to_extract = "Total bounding box for atom centers:"
    line_of_interest = list(filter(lambda l: line_to_extract in l, leapl_log_lines))
    d1, d2, d3 = line_of_interest[-1].split(line_to_extract)[1].split()
    d1, d2, d3 = float(d1), float(d2), float(d3)
    assert d1 == d2 == d3
    # scale the d since after minimisation the system turns out to be much smaller
    d = d1 * 0.8
    return {
        'cbv1': d, 'cbv2': 0, 'cbv3': 0,
        'cbv4': (-1/3.0)*d, 'cbv5': (2/3.0)*np.sqrt(2)*d, 'cbv6': 0,
        'cbv7': (-1/3.0)*d, 'cbv8': (-1/3.0)*np.sqrt(2)*d, 'cbv9': (-1/3)*np.sqrt(6)*d,
    }


def prepare_antechamber_parmchk2(source_script, target_script, net_charge):
    """
    Prepare the ambertools scripts.
    Particularly, create the scritp so that it has the net charge
    # fixme - run antechamber directly with the right settings from here?
    # fixme - check if antechamber has a python interface?
    """
    net_charge_set = open(source_script).read().format(net_charge=net_charge)
    open(target_script, 'w').write(net_charge_set)


def set_charges_from_ac(mol2_filename, ac_ref_filename):
    # ! the mol file will be overwritten
    print('Overwriting the mol2 file with charges from the ac file')
    # load the charges from the .ac file
    ac_atoms, _ = get_atoms_bonds_from_ac(ac_ref_filename)
    # load the .mol2 files with MDAnalysis and correct the charges
    mol2 = load_mol2_wrapper(mol2_filename)

    for mol2_atom in mol2.atoms:
        found_match = False
        for ac_atom in ac_atoms:
            if mol2_atom.name.upper() == ac_atom.name.upper():
                found_match = True
                mol2_atom.charge = ac_atom.charge
                break
        assert found_match, "Could not find the following atom in the AC: " + mol2_atom.name

    # update the mol2 file
    mol2.atoms.write(mol2_filename)


def set_charges_from_mol2(mol2_filename, mol2_ref_filename, by_atom_name=False, by_index=False, by_general_atom_type=False):
    # mol2_filename will be overwritten!
    print(f'Overwriting {mol2_filename} mol2 file with charges from {mol2_ref_filename} file')
    # load the ref charges
    ref_mol2 = load_mol2_wrapper(mol2_ref_filename)
    # load the .mol2 files with MDAnalysis and correct the charges
    mol2 = load_mol2_wrapper(mol2_filename)

    if by_atom_name and by_index:
        raise ValueError('Cannot have both. They are exclusive')
    elif not by_atom_name and not by_index:
        raise ValueError('Either option has to be selected.')

    # save the sum of charges before
    ref_sum_q =  sum(a.charge for a in ref_mol2.atoms)

    if by_atom_name:
        for mol2_atom in mol2.atoms:
            if mol2_atom.type == 'DU':
                continue

            found_match = False
            for ref_atom in ref_mol2.atoms:
                if ref_atom.type == 'DU':
                    continue

                if mol2_atom.name.upper() == ref_atom.name.upper():
                    if found_match == True:
                        raise Exception('AtomNames are not unique or do not match')
                    found_match = True
                    mol2_atom.charge = ref_atom.charge
            assert found_match, "Could not find the following atom: " + mol2_atom.name
    elif by_general_atom_type:
        for mol2_atom in mol2.atoms:
            found_match = False
            for ref_atom in ref_mol2.atoms:
                if element_from_type[mol2_atom.type.upper()] == element_from_type[ref_atom.type.upper()]:
                    if found_match == True:
                        raise Exception('AtomNames are not unique or do not match')
                    found_match = True
                    mol2_atom.charge = ref_atom.charge
            assert found_match, "Could not find the following atom in the AC: " + mol2_atom.name
    elif by_index:
        for mol2_atom, ref_atom in zip(mol2.atoms, ref_mol2.atoms):
                atype = element_from_type[mol2_atom.type.upper()]
                reftype = element_from_type[ref_atom.type.upper()]
                if atype != reftype:
                    raise Exception(f"The found general type {atype} does not equal to the reference type {reftype} ")

                mol2_atom.charge = ref_atom.charge

    assert ref_sum_q == sum(a.charge for a in mol2.atoms)
    # update the mol2 file
    mol2.atoms.write(mol2_filename)


def get_protein_net_charge(working_dir, protein_file, ambertools_tleap, leap_input_file, prot_ff):
    """
    Use automatic ambertools solvation of a single component to determine what is the next charge of the system.
    This should be replaced with pka/propka or something akin.
    Note that this is unsuitable for the hybrid ligand: ambertools does not understand a hybrid ligand
    and might assign the wront net charge.
    """
    cwd = working_dir / 'prep' / 'prep_protein_to_find_net_charge'
    if not cwd.is_dir():
        cwd.mkdir()

    # copy the protein
    shutil.copy(working_dir / protein_file, cwd)

    # use ambertools to solvate the protein: set ion numbers to 0 so that they are determined automatically
    # fixme - consider moving out of the complex
    leap_in_conf = open(leap_input_file).read()
    ligand_ff = 'leaprc.gaff' # ignored but must be provided
    open(cwd / 'solv_prot.in', 'w').write(leap_in_conf.format(protein_ff=prot_ff, ligand_ff=ligand_ff,
                                                                          protein_file=protein_file))

    log_filename = cwd / "ties_tleap.log"
    with open(log_filename, 'w') as LOG:
        try:
            subprocess.run([ambertools_tleap, '-s', '-f', 'solv_prot.in'],
                           cwd = cwd,
                           stdout=LOG, stderr=LOG,
                           check=True, text=True,
                           timeout=60 * 2  # 2 minutes
                        )
        except subprocess.CalledProcessError as E:
            print('ERROR: tleap could generate a simple topology for the protein to check the number of ions. ')
            print(f'ERROR: The output was saved in the directory: {cwd}')
            print(f'ERROR: can be found in the file: {log_filename}')
            raise E


    # read the file to see how many ions were added
    u = mda.Universe(cwd / 'prot_solv.pdb')
    cl = len(u.select_atoms('name Cl-'))
    na = len(u.select_atoms('name Na+'))
    if cl > na:
        return cl-na
    elif cl < na:
        return -(na-cl)

    return 0


def prepareFile(src, dst, symbolic=True):
    """
    Either copies or sets up a relative link between the files.
    This allows for a quick switch in how the directory structure is organised.
    Using relative links means that the entire TIES ligand or TIES complex
    has to be moved together.
    However, one might want to be able to send a single replica anywhere and
    execute it independantly (suitable for BOINC).

    @type: 'sym' or 'copy'
    """
    if symbolic:
        # note that deleting all the files is intrusive, todo
        if os.path.isfile(dst):
            os.remove(dst)
        os.symlink(src, dst)
    else:
        if os.path.isfile(dst):
            os.remove(dst)
        shutil.copy(src, dst)



def set_coor_from_ref_by_named_pairs(mol2_filename, coor_ref_filename, output_filename, left_right_pairs_filename):
    """
    Set coordinates but use atom names provided by the user.

    Example of the left_right_pairs_filename content:
    # flip the first ring
    # move the first c and its h
    C32 C18
    H34 C19
    # second pair
    C33 C17
    # the actual matching pair
    C31 C16
    H28 H11
    # the second matching pair
    C30 C15
    H29 H12
    #
    C35 C14
    # flip the other ring with itself
    C39 C36
    C36 C39
    H33 H30
    H30 H33
    C37 C38
    C38 C37
    H31 H32
    H32 H31
    """
    # fixme - check if the names are unique

    # parse "left_right_pairs_filename
    # format per line: leftatom_name right_atom_name
    lines = open(left_right_pairs_filename).read().split(os.linesep)
    left_right_pairs = (l.split() for l in lines if not l.strip().startswith('#'))

    # load the ref coordinates
    ref_mol2 = load_mol2_wrapper(coor_ref_filename)
    # load the .mol2 files with MDAnalysis and correct the charges
    static_mol2 = load_mol2_wrapper(mol2_filename)
    # this is being modified
    mod_mol2 = load_mol2_wrapper(mol2_filename)


    for pair in left_right_pairs:
        print('find pair', pair)
        new_pos = False
        for mol2_atom in static_mol2.atoms:
            # check if we are assigning from another molecule
            for ref_atom in ref_mol2.atoms:
                if mol2_atom.name.upper() == pair[0] and ref_atom.name.upper() == pair[1]:
                    new_pos = ref_atom.position
            # check if we are trying to assing coords from the same molecule
            for another_atom in static_mol2.atoms:
                if mol2_atom.name.upper() == pair[0] and another_atom.name.upper() == pair[1]:
                    new_pos = another_atom.position

        if new_pos is False:
            raise Exception("Could not find this pair: " + str(pair))

        # assign the position to the right atom
        # find pair[0]
        found = False
        for atom in mod_mol2.atoms:
            if atom.name.upper() == pair[0]:
                atom.position = new_pos
                found = True
                break
        assert found


    # update the mol2 file
    mod_mol2.atoms.write(output_filename)


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


def update_PBC_in_namd_input(namd_filename, new_pbc_box, structure_filename, constraint_lines=''):
    """
    fixme - rename this file since it generates the .eq files
    These are the lines we modify:
    cellBasisVector1	{cell_x}  0.000  0.000
    cellBasisVector2	 0.000  {cell_y}  0.000
    cellBasisVector3	 0.000  0.000 {cell_z}

    With x/y/z replacing the 3 values
    """
    assert len(new_pbc_box) == 3

    reformatted_namd_in = open(namd_filename).read().format(
        cell_x=new_pbc_box[0], cell_y=new_pbc_box[1], cell_z=new_pbc_box[2],
        constraints=constraint_lines, output='test_output', structure=structure_filename)

    # write to the file
    open(namd_filename, 'w').write(reformatted_namd_in)

def create_constraint_files(original_pdb, output):
    '''

    :param original_pdb:
    :param output:
    :return:
    '''
    u = mda.Universe(original_pdb)
    # for each atom, give the B column the right value
    for atom in u.atoms:
        # ignore water
        if atom.resname in ['WAT', 'Na+', 'TIP3W', 'TIP3']:
            continue

        # set each atom depending on whether it is a H or not
        if atom.name.upper().startswith('H'):
            atom.tempfactor = 0
        else:
            # restrain the heavy atom
            atom.tempfactor = 4

    u.atoms.write(output)


def init_namd_file_min(from_dir, to_dir, filename, structure_name, pbc_box):
    '''

    :param from_dir:
    :param to_dir:
    :param filename:
    :param structure_name:
    :param pbc_box:
    :return:
    '''
    min_namd_initialised = open(os.path.join(from_dir, filename)).read() \
        .format(structure_name=structure_name, **pbc_box)
    out_name = 'eq0.conf'
    open(os.path.join(to_dir, out_name), 'w').write(min_namd_initialised)

def generate_namd_prod(namd_prod, dst_dir, structure_name):
    '''

    :param namd_prod:
    :param dst_dir:
    :param structure_name:
    :return:
    '''
    input_data = open(namd_prod).read()
    reformatted_namd_in = input_data.format(output='sim1', structure_name=structure_name)
    open(dst_dir, 'w').write(reformatted_namd_in)


def generate_namd_eq(namd_eq, dst_dir, structure_name):
    '''

    :param namd_eq:
    :param dst_dir:
    :param structure_name:
    :return:
    '''
    input_data = open(namd_eq).read()
    for i in range(1,3):

        if i == 1:
            run = """
constraintScaling 1
run 10000
            """
        else:
            run = """
while {$n <= $nall} {
   constraintScaling $factor
   run 40000
   set n [expr $n + 1]
   set factor [expr $factor * 0.5]
}

constraintScaling 0
run 600000
            """

        constraints = f"""
constraints  on
consexp  2
# use the same file for the position reference and the B column
consref  ../build/{structure_name}.pdb ;#need all positions
conskfile  ../build/cons.pdb
conskcol  B
        """

        prev_output = 'eq{}'.format(i-1)

        reformatted_namd_in = input_data.format(
            constraints=constraints, output='eq%d' % (i),
            prev_output=prev_output, structure_name=structure_name, run=run)

        next_eq_step_filename = dst_dir / ("eq%d.conf" % (i))
        open(next_eq_step_filename, 'w').write(reformatted_namd_in)


def redistribute_charges(mda):
    """
    Calculate the original charges in the matched component.
    """


    return
