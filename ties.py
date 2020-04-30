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
hpc_submit = 'hpc_surfsara_single_core_ligands.sh'
left_ligand = 'left_coor.pdb'
right_ligand = 'right_coor.pdb'
protein_filename = 'protein.pdb'
net_charge = 0
# force_mismatch_list = [('O2', 'O4'), ('N3', 'N6')]
# force_mismatch_list = None
# rather than using the empirical antechamber -c bcc, copy agastya's values
use_agastyas_charges = True
left_charges = 'left_q.mol2'
right_charges = 'right_q.mol2'

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
# copy the ambertools.sh for 1) creating .mol2 - antechamber, 2) optimising the structure with sqm
script_dir = code_root / PurePosixPath('scripts')
namd_script_dir = script_dir / 'namd'
ambertools_script_dir = script_dir / 'ambertools'

# fixme - replace with the python API
antechamber_sqm_script_name = 'assign_charge_parameters.sh'
assign_charges_antechamber_script_filename = 'assign_charge_parameters_antechamber.sh'
assign_charges_parmchk2_script_filename = 'assign_charge_parameters_parmchk2.sh'

if not use_agastyas_charges:
    prepare_antechamber_parmchk2(ambertools_script_dir / antechamber_sqm_script_name,
                            workplace_root / antechamber_sqm_script_name, net_charge=net_charge)
    # execute the script (the script has to source amber.sh)
    # do not do this if there is a .frcmod files
    if not (workplace_root / 'left.frcmod').is_file():
        # fixme - add checks for other files
        output = subprocess.check_output(['sh', workplace_root / antechamber_sqm_script_name, left_ligand, right_ligand])
        # todo - CHECK IF THE RESULTS ARE CORRECT
elif use_agastyas_charges:
    # fixme - create a file that explains that you copied the charges

    #
    prepare_antechamber_parmchk2(ambertools_script_dir/assign_charges_antechamber_script_filename,
                                 workplace_root/assign_charges_antechamber_script_filename, net_charge=net_charge)

    # create the .mol2 files
    if not (workplace_root / 'left.frcmod').is_file():
        try:
            output = subprocess.check_output(['sh', workplace_root / assign_charges_antechamber_script_filename,
                                              left_ligand, right_ligand])
        except Exception as E:
            print(E.output)
            raise E

    # take the .mol2 file and correct the charges to reflect Agastya's
    set_charges_from_mol2(workplace_root / 'left.mol2', workplace_root / left_charges)
    set_charges_from_mol2(workplace_root / 'right.mol2', workplace_root / right_charges)

    # copy and execute parmchk2 to generate the .frcmod
    shutil.copy(ambertools_script_dir/assign_charges_parmchk2_script_filename, workplace_root)
    if not (workplace_root / 'left.frcmod').is_file():
        subprocess.check_output(['sh', workplace_root/assign_charges_parmchk2_script_filename, 'left', 'right'])


# load the files (.mol2) and superimpose the two topologies
# fixme - superimpose the molecules or stop relaying on RMSD info
# fixme - call any of the tools you have (antechamber, parmchk2)
suptop, mda_l1, mda_l2 = getSuptop(workplace_root / 'left.mol2', workplace_root / 'right.mol2')

# save the results of the topology superimposition as a json
top_sup_joint_meta = workplace_root / 'joint_meta_fep.json'
save_superimposition_results(top_sup_joint_meta, suptop)
write_dual_top_pdb(workplace_root / 'morph.pdb', mda_l1, mda_l2, suptop)
# save the merged topologies as a .mol2 file
top_merged_filename = workplace_root / 'morph.mol2'
write_merged(suptop, top_merged_filename)

# write a copy of the morph.pdb - fixme
# correct_fep_tempfactor(top_sup_joint_meta, morph_solv, morph_solv_fep)

# this generates the "merged_solvated.pdb" which does not have .fep information in the .pdb tempfactor columns
morph_solv = workplace_root / "morph_solv.pdb"
morph_solv_fep = workplace_root / "morph_solv_fep.pdb"

# check if the .frcmod were generated
left_frcmod = workplace_root / 'left.frcmod'
right_frcmod = workplace_root / 'right.frcmod'
if not left_frcmod.is_file():
    sys.exit(5)
elif not right_frcmod.is_file():
    sys.exit(5)

# generate the joint .frcmod file
merged_frc_filename = workplace_root / 'morph.frcmod'
join_frcmod_files2(left_frcmod, right_frcmod, merged_frc_filename)
# fixme
updated_frcmod = check_hybrid_frcmod(top_merged_filename, merged_frc_filename, '/home/dresio/software/amber18install/bin/tleap', 'gaff')
with open(merged_frc_filename, 'w') as FOUT:
    FOUT.write(updated_frcmod)

# copy the solvate script for tleap
shutil.copy(ambertools_script_dir / "run_tleap.sh", workplace_root)
shutil.copy(ambertools_script_dir / "leap.in", workplace_root)
# solvate using AmberTools, copy leap.in and use tleap
try:
    output = subprocess.check_output(['sh', workplace_root / "run_tleap.sh"])
    assert "Errors = 0;" in str(output), "Errors when running tleap: " + str(output)
except subprocess.CalledProcessError as E:
    print('Error occured when running tleap on the ligand: ', E.output)
    sys.exit(2)
# make a copy of the tleap generated topology file with a more useful extension
if (workplace_root / "morph_solv.prmtop").is_file():
    (workplace_root / "morph_solv.prmtop").unlink() # remove the file
shutil.copy(workplace_root / "morph_solv.top", workplace_root / "morph_solv.prmtop")

# this is applied to the solvated version
correct_fep_tempfactor(top_sup_joint_meta, morph_solv, morph_solv_fep)

# take care of the ligand-ligand without the protein
liglig_workplace = workplace_root / 'lig'
if not liglig_workplace.is_dir():
    liglig_workplace.mkdir()

# generate the constraint files in liglig
constraint_files = create_4_constraint_files(morph_solv, liglig_workplace)
pbc_dimensions = get_PBC_coords(morph_solv)
init_namd_file_min(namd_script_dir, liglig_workplace, "min.namd",
                           structure_name='morph_solv', pbc_box=pbc_dimensions)
# prepare the namd eq
eq_namd_files = generate_namd_eq(namd_script_dir / "eq.namd", liglig_workplace)
# namd production protocol
prod_filename = init_namd_file_prod(namd_script_dir, liglig_workplace, "prod.namd", structure_name='morph_solv')

# Generate the directory structure for all the lambdas, and copy the files
for lambda_step in [0, 0.05] + list(np.linspace(0.1, 0.9, 9)) + [0.95, 1]:
    lambda_path = liglig_workplace / f'lambda_{lambda_step:.2f}'
    if not os.path.exists(lambda_path):
        os.makedirs(lambda_path)

    # for each lambda create 5 replicas
    for replica_no in range(1, 5 + 1):
        replica_dir = lambda_path / f'rep{replica_no}'
        if not os.path.exists(replica_dir):
            os.makedirs(replica_dir)
            # set the lambda value for the directory

        open(replica_dir / 'lambda', 'w').write(f'{lambda_step:.2f}')

        # copy the files for each simulation
        # "coordinates" with the fep metadata
        shutil.copy(morph_solv_fep, replica_dir)
        # copy the ambertools generated topology file
        shutil.copy(workplace_root / "morph_solv.top", replica_dir)
        # the normal coordinates
        shutil.copy(workplace_root / "morph_solv.pdb", replica_dir)

        # copy the restraint files, which use b column
        [shutil.copy(constraint_file, replica_dir) for constraint_file in constraint_files]

        # copy NAMD protocol
        shutil.copy(liglig_workplace / 'min.namd', replica_dir)
        [shutil.copy(eq, replica_dir) for eq in eq_namd_files]
        shutil.copy(prod_filename, replica_dir)

        # todo - create/copy the 4 different EQ files
        # copy the submit script
        shutil.copy(namd_script_dir / hpc_submit, replica_dir / 'submit.sh')

# copy the scheduler to the main directory
shutil.copy(script_dir / "schedule_separately.py", liglig_workplace)
shutil.copy(script_dir /"check_namd_outputs.py", liglig_workplace)

##########################################################
# ------------------ complex-complex --------------

# fixme - use tleap to merge+solvate, decide on the charges?

complex_workplace = workplace_root / 'complex'
if not complex_workplace.is_dir():
    complex_workplace.mkdir()

# prepare the simulation files
# copy the complex .pdb
shutil.copy(workplace_root / protein_filename, complex_workplace)
# todo - ensure the protein protonation is correct
# todo - call antechamber?

# copy the ligand (the morphed ligand), and its .frcmod
shutil.copy(top_merged_filename, complex_workplace)
shutil.copy(merged_frc_filename, complex_workplace)

# todo - dock with the ligand to create a complex
# copy the protein tleap input file (ambertools)
shutil.copy(ambertools_script_dir / 'leap_complex.in', complex_workplace)
shutil.copy(ambertools_script_dir / 'run_tleap_complex.sh', complex_workplace)

# solvate the complex (tleap, ambertools)
#      rn tleap also combines complex+ligand, and generates amberparm
try:
    output = subprocess.check_output(['sh', "run_tleap_complex.sh"], cwd=complex_workplace)
except Exception as ex:
    print(ex.output)
    raise ex
assert 'Errors = 0' in str(output)

# tleap generates these
complex_solvated = complex_workplace / 'morph_solv.pdb'

# update the complex to create complex.fep file
# generate the merged .fep file
complex_solvated_fep = complex_workplace / 'morph_solv_fep.pdb'
correct_fep_tempfactor(top_sup_joint_meta, complex_solvated, complex_solvated_fep)

# fixme - ensure that the _fep is only applied to the ligand, not the protein,
# fixme - check that the protein does not have the same resname

# get the PBC data from MDAnalysis
solv_box_complex_pbc = get_PBC_coords(complex_solvated)

# copy the NAMD input files to the main directory first
# prepare the namd minmisation
init_namd_file_min(namd_script_dir, complex_workplace, "min.namd",
                           structure_name='morph_solv', pbc_box=solv_box_complex_pbc)
eq_namd_filenames = generate_namd_eq(namd_script_dir / "eq.namd", complex_workplace)
init_namd_file_prod(namd_script_dir, complex_workplace, "prod.namd", structure_name='morph_solv')

# generate the 4 different constrain .pdb files files, which use b column
constraint_files = create_4_constraint_files(complex_solvated, complex_workplace)

# Generate the directory structure for all the lambdas, and copy the files
for lambda_step in [0, 0.05] + list(np.linspace(0.1, 0.9, 9)) + [0.95, 1]:
    lambda_path = complex_workplace / f'lambda_{lambda_step:.2f}'
    if not os.path.exists(lambda_path):
        os.makedirs(lambda_path)

    # for each lambda create 5 replicas
    for replica_no in range(1, 5 + 1):
        replica_dir = lambda_path / f'rep{replica_no}'
        if not os.path.exists(replica_dir):
            os.makedirs(replica_dir)

        # set the lambda value for the directory
        open(replica_dir / 'lambda', 'w').write(f'{lambda_step:.2f}')

        # copy all the necessary files
        shutil.copy(complex_solvated, replica_dir)
        shutil.copy(complex_solvated_fep, replica_dir)
        # copy the ambertools generated topology
        shutil.copy(complex_workplace / "morph_solv.top", replica_dir)
        # copy the .pdb files with constraints in the B column
        [shutil.copy(constraint_file, replica_dir) for constraint_file in constraint_files]

        # copy the NAMD protocol files
        shutil.copy(complex_workplace / "min.namd", replica_dir)
        [shutil.copy(eq, replica_dir) for eq in eq_namd_filenames]
        shutil.copy(complex_workplace / "prod.namd", replica_dir)

        # copy the surfsara submit script - fixme - make this general
        shutil.copy(namd_script_dir / hpc_submit, replica_dir / 'submit.sh')

# copy the scheduler to the main directory
shutil.copy(script_dir / "schedule_separately.py", complex_workplace)
shutil.copy(script_dir / "check_namd_outputs.py", complex_workplace)

# fixme States show the progress of the simulation.

# use the preprepared pdb complex with the ligand
# solvate the preprepared pdb complex with the ligand
# generate all the merged files


