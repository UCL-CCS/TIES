"""
Load two ligands, run the topology superimposer, and then using the results, generate the NAMD input files.
"""
from generator import *
import topology_superimposer

import os
import shutil
import sys
import subprocess
from pathlib import Path, PurePosixPath


# fixme - turn into a function and give it the hpc submit
hpc_submit = None #  "hpc_hartree_hsp.sh"
left_ligand = 'left_coor.pdb'
right_ligand = 'right_coor.pdb'
protein_filename = 'protein.pdb'
net_charge = -1
# force_mismatch_list = [('O2', 'O4'), ('N3', 'N6')] # None
# rather than using the empirical antechamber -c bcc, copy agastya's values
use_agastyas_charges = False
left_charges = 'left_q.mol2'
right_charges = 'right_q.mol2'
# the coordinates change slightly after antechamber, reassign the coordinates to the .mol2
use_original_coor = True



# set the working directory to the one where the script is called
workplace_root = Path(os.getcwd())
print('Working Directory: ', workplace_root)
# set the path to the scripts
code_root = Path(os.path.dirname(__file__))

# check if the input files are in the directory
if not (workplace_root / left_ligand).is_file():
    print(f'File {left_ligand} not found in {workplace_root}')
    sys.exit(1)
elif not (workplace_root / right_ligand).is_file():
    print(f'File {right_ligand} not found in {workplace_root}')
    sys.exit(1)
elif not (workplace_root / protein_filename).is_file():
    print(f'File {protein_filename} not found in {workplace_root}')
    sys.exit(1)

# copy ambertools scripts to create .mol2 with antechamber (and sqm)
script_dir = code_root / PurePosixPath('scripts')
namd_script_dir = script_dir / 'namd'
ambertools_script_dir = script_dir / 'ambertools'

# fixme - replace with the python API
assign_charges_antechamber_script = 'assign_charge_parameters_antechamber.sh'
check_charges_parmchk2_script = 'call_parmchk2.sh'

# antechamber creates the left.mol2 and right.mol2 files and assigns BCC charges
prepare_antechamber_parmchk2(ambertools_script_dir / assign_charges_antechamber_script,
                             workplace_root / assign_charges_antechamber_script, net_charge=net_charge)
if not (workplace_root / 'left.mol2').is_file() or not (workplace_root / 'right.mol2').is_file():
    try:
        print('Generating BCC charges with ambertool')
        output = subprocess.check_output(['sh', workplace_root / assign_charges_antechamber_script,
                                          left_ligand, right_ligand])
    except Exception as E:
        print('Ambertools antechamber could not generate at least one of the .mol2 files')
        print(E)
        raise E

if use_agastyas_charges:
    print('Copying REST charges from left_q.mol2 and right_q.mol2')
    # make a copy of the files before modifying them
    shutil.copy(workplace_root / 'left.mol2', workplace_root / 'left_before_RESP.mol2')
    shutil.copy(workplace_root / 'right.mol2', workplace_root / 'right_before_RESP.mol2')
    # take the .mol2 file and correct the charges to reflect Agastya's RESP
    set_charges_from_mol2(workplace_root / 'left.mol2', workplace_root / left_charges)
    set_charges_from_mol2(workplace_root / 'right.mol2', workplace_root / right_charges)

if use_original_coor:
    print(f'Copying coordinates from {left_ligand} and {right_ligand} since antechamber changes them slightly')
    # copy the files before applying the coordinates
    shutil.copy(workplace_root / 'left.mol2', workplace_root / 'left_before_COOR.mol2')
    shutil.copy(workplace_root / 'right.mol2', workplace_root / 'right_before_COOR.mol2')
    set_coor_from_ref(workplace_root / 'left.mol2', workplace_root / left_ligand)
    set_coor_from_ref(workplace_root / 'right.mol2', workplace_root / right_ligand)

# rename the atom names to ensure they are unique across the two molecules
ensureUniqueAtomNames(workplace_root / 'left.mol2', workplace_root / 'right.mol2')

# superimpose the two topologies
suptop, mda_l1, mda_l2 = getSuptop(workplace_root / 'left.mol2', workplace_root / 'right.mol2')

# save the superimposition results
left_right_matching_json = workplace_root / 'joint_meta_fep.json'
save_superimposition_results(left_right_matching_json, suptop)
# hybrid .pdb
write_dual_top_pdb(workplace_root / 'morph.pdb', mda_l1, mda_l2, suptop)
# hybrid .mol2
hybrid_mol2 = workplace_root / 'morph.mol2'
write_merged(suptop, hybrid_mol2)

# use parmchk2 to generate the .frcmod
shutil.copy(ambertools_script_dir / check_charges_parmchk2_script, workplace_root)
subprocess.check_output(['sh', workplace_root / check_charges_parmchk2_script, 'left', 'right'])
if not (workplace_root / 'left.frcmod').is_file() or \
        not (workplace_root / 'right.frcmod').is_file():
    raise Exception('Ambertools sqm did not generate at least one of the .frcmod files')

# are .frcmod files generated?
left_frcmod = workplace_root / 'left.frcmod'
right_frcmod = workplace_root / 'right.frcmod'
if not left_frcmod.is_file() or not right_frcmod.is_file():
    raise Exception(f'ERROR: Ambertools sqm generation of at least one of the .frcmod files failed. File not found. '
          f'Please check antechamber log for more details. ')

# join the .frcmod files
hybrid_frcmod = workplace_root / 'morph.frcmod'
join_frcmod_files2(left_frcmod, right_frcmod, hybrid_frcmod)

# the hybrid .frcmod might not contains "new terms" between the appearing/disappearing atoms.
# In that case, insert fake terms
# fixme - clean up this part: remove the direct call to tleap?
updated_frcmod = check_hybrid_frcmod(hybrid_mol2, hybrid_frcmod, '/home/dresio/software/amber18install/bin/tleap', 'gaff')
with open(hybrid_frcmod, 'w') as FOUT:
    FOUT.write(updated_frcmod)

# ---------------------
def prepare_inputs(workplace_root, directory='complex', protein=None,
                   hybrid_mol2=None, hybrid_frc=None,
                   left_right_mapping=None,
                   namd_script_loc=None,
                   submit_script=None,
                   scripts_loc=None,
                   tleap_in=None):
    # fixme - use tleap to merge+solvate, decide on the charges?

    # fixme rename
    dest_dir = workplace_root / directory
    if not dest_dir.is_dir():
        dest_dir.mkdir()

    # copy the protein complex .pdb
    if protein is not None:
        shutil.copy(workplace_root / protein, dest_dir)

    # copy the hybrid ligand (topology and .frcmod)
    shutil.copy(hybrid_mol2, dest_dir)
    shutil.copy(hybrid_frc, dest_dir)

    # copy the protein tleap input file (ambertools)
    shutil.copy(ambertools_script_dir / tleap_in, dest_dir / 'leap.in')
    shutil.copy(ambertools_script_dir / 'run_tleap.sh', dest_dir)

    try:
        # tleap: combine ligand+complex, solvate, generate amberparm
        output = subprocess.check_output(['sh', "run_tleap.sh"], cwd=dest_dir)
        assert "Errors = 0;" in str(output), "Errors when running tleap: " + str(output)
        # this file should have been generated by tleap
        morph_solv = dest_dir / "morph_solv.pdb"
        if not morph_solv.is_file():
            raise Exception(f'During the solvation in {dest_dir}, tleap did not generate the {morph_solv} file')
    except Exception as E:
        print('Error occured when running tleap script to solvate ligand: ', E)
        raise E
    assert 'Errors = 0' in str(output)

    # tleap generated:
    complex_solvated = dest_dir / 'morph_solv.pdb'

    # generate the merged .fep file
    complex_solvated_fep = dest_dir / 'morph_solv_fep.pdb'
    correct_fep_tempfactor(left_right_mapping, complex_solvated, complex_solvated_fep)

    # fixme - check that the protein does not have the same resname?

    # calculate PBC for an octahedron
    solv_oct_boc = extract_PBC_oct_from_tleap_log(dest_dir / "leap.log")

    # prepare NAMD input files for min+eq+prod
    init_namd_file_min(namd_script_loc, dest_dir, "min.namd",
                       structure_name='morph_solv', pbc_box=solv_oct_boc)
    eq_namd_filenames = generate_namd_eq(namd_script_loc / "eq.namd", dest_dir)
    shutil.copy(namd_script_loc / "prod.namd", dest_dir)

    # generate 4 different constraint .pdb files (it uses B column)
    constraint_files = create_4_constraint_files(complex_solvated, dest_dir)

    # Generate the directory structure for all the lambdas, and copy the files
    for lambda_step in [0, 0.05] + list(np.linspace(0.1, 0.9, 9)) + [0.95, 1]:
        lambda_path = dest_dir / f'lambda_{lambda_step:.2f}'
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

            # copy the necessary files
            shutil.copy(complex_solvated, replica_dir)
            shutil.copy(complex_solvated_fep, replica_dir)
            # copy ambertools-generated topology
            shutil.copy(dest_dir / "morph_solv.top", replica_dir)
            # copy the .pdb files with constraints in the B column
            [shutil.copy(constraint_file, replica_dir) for constraint_file in constraint_files]

            # copy the NAMD protocol files
            shutil.copy(dest_dir / "min.namd", replica_dir)
            [shutil.copy(eq, replica_dir) for eq in eq_namd_filenames]
            shutil.copy(dest_dir / "prod.namd", replica_dir)

            # copy the surfsara submit script - fixme - make this general
            if submit_script is not None:
                shutil.copy(namd_script_loc / submit_script, replica_dir / 'submit.sh')

    # copy handy scripts to the main directory
    shutil.copy(scripts_loc / "schedule_separately.py", dest_dir)
    shutil.copy(namd_script_loc / "check_namd_outputs.py", dest_dir)
    shutil.copy(namd_script_loc / "extract_energies.py", dest_dir)

##########################################################
# ------------------   LIGAND ----------------------

prepare_inputs(workplace_root, directory='lig',
               protein=None,
               hybrid_mol2=hybrid_mol2,
               hybrid_frc=hybrid_frcmod,
               left_right_mapping=left_right_matching_json,
               namd_script_loc=namd_script_dir,
               scripts_loc=script_dir,
               tleap_in='leap_ligand.in')

##########################################################
# ------------------ complex-complex --------------

prepare_inputs(workplace_root, directory='complex',
               protein=protein_filename,
               hybrid_mol2=hybrid_mol2,
               hybrid_frc=hybrid_frcmod,
               left_right_mapping=left_right_matching_json,
               namd_script_loc=namd_script_dir,
               scripts_loc=script_dir,
               tleap_in='leap_complex.in')

