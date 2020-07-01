"""
Load two ligands, run the topology superimposer, and then using the results, generate the NAMD input files.
"""
import os
import shutil
import sys
import tempfile
import subprocess
import math
import argparse
from pathlib import Path, PurePosixPath

from generator import *


# fixme - turn into a function and give it the hpc submit
hpc_submit = None #  "hpc_hartree_hsp.sh"
left_ligand = 'left_coor.pdb'
right_ligand = 'right_coor.pdb'
protein_filename = 'protein.pdb'
net_charge = -2
# pair_q_atol = 0.101 # in electrons
pair_q_atol = 0.131 # in electrons
manually_matched = None # ['C2', 'C3']
force_mismatch = None # [('C9', 'C60')]

use_agastyas_charges = True
# the coordinates change slightly after antechamber, reassign the coordinates to the .mol2
use_original_coor = True

# ambertools forcefiled
# amber_forcefield = "leaprc.protein.ff14SB"
amber_forcefield = "leaprc.ff99SBildn" # used by Agastya before, latest is: "leaprc.protein.ff14SB"
atom_type = 'gaff' # fixme - check?
if use_agastyas_charges:
    left_charges = 'left_q.mol2'
    right_charges = 'right_q.mol2'
    # ignore the charges,
    charge_type = 'dc'  # 'Delete Charge'
    # use berendsen
else:
    charge_type = 'bcc'  # 'AM1-BCC'

# amber_ligand_ff = "leaprc.gaff" # fixme - check
namd_prod = "prod_2017.namd"    # only different because uses Berendsen
# namd_prod = "prod.namd"
# namd_prod = "prod_2017_manual.namd"    # only different because uses Berendsen
align_molecules = False
ignore_protein = False

# First interface
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TIES 20')
    parser.add_argument('action', metavar='command', type=str,
                        help='Action to be performed. E.g. "ties rename .." ')
    parser.add_argument('-l', '--left-ligand', metavar='Left_Ligand_File', dest='left_ligand',
                        type=str, required=False,
                        help='The left ligand filename')
    parser.add_argument('-r', '--right_ligand', metavar='Right_Ligand_File', dest='right_ligand',
                        type=str, required=False,
                        help='The right ligand filename')
    # protein location
    # net charge (for bcc)
    # cwd

    args = parser.parse_args()

    # ------------------ Software Configuration
    # set the working directory to the one where the script is called
    workplace_root = Path(os.getcwd())
    print('Working Directory: ', workplace_root)
    # set the path to the scripts
    code_root = Path(os.path.dirname(__file__))

    # conf ambertools
    # todo - check if $AMBERHOME is set
    ambertools_bin = PurePosixPath("/home/dresio/software/amber18install/bin/")
    # directory with scripts/input files
    script_dir = code_root / PurePosixPath('scripts')
    namd_script_dir = script_dir / 'namd'
    ambertools_script_dir = script_dir / 'ambertools'

    # check if the input files are in the directory
    if not (workplace_root / left_ligand).is_file():
        print(f'File {left_ligand} not found in {workplace_root}')
        sys.exit(1)
    elif not (workplace_root / right_ligand).is_file():
        print(f'File {right_ligand} not found in {workplace_root}')
        sys.exit(1)

    if args.action == 'rename':
        # check if ligand files are fine
        if args.left_ligand is None or args.right_ligand is None:
            print('Please supply files for renaming with -l and -r')
            sys.exit()
        usr_ll = args.left_ligand
        usr_rl = args.right_ligand
        print('Atom names will be renamed to ensure that the atom names are unique across the two molecules.')
        renameAtomNamesUnique(usr_ll, usr_rl)
        sys.exit()

    if args.action == 'create':
        print('Main protocol will be used')
    else:
        sys.exit()


# ------------------------------------- Start TIES
if not ignore_protein and not (workplace_root / protein_filename).is_file():
    print(f'File {protein_filename} not found in {workplace_root}')
    sys.exit(1)

# subprocess options for calling ambertools
subprocess_kwargs = {
    "check" : True, "text" : True,
    "cwd" : workplace_root,
    "stdout" : sys.stdout, #subprocess.PIPE,
    "stderr" : sys.stderr, #subprocess.PIPE,
    "timeout" : 60 * 10 # 10 minute timeout
}

# prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
if not (workplace_root / 'left.mol2').is_file() or not (workplace_root / 'right.mol2').is_file():
    print('Ambertools antechamber stage: converting to .mol2 and generating charges')
    subprocess.run([ambertools_bin / 'antechamber', '-i', left_ligand, '-fi', 'pdb',
                                        '-o', 'left.mol2', '-fo', 'mol2',
                                        '-c', charge_type, '-at', atom_type,
                                        '-nc', str(net_charge)],
                               **subprocess_kwargs)

# prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
if not (workplace_root / 'right.mol2').is_file():
    print('Ambertools antechamber stage: converting to .mol2 and generating charges')
    subprocess.run([ambertools_bin / 'antechamber', '-i', right_ligand, '-fi', 'pdb',
                                '-o', 'right.mol2', '-fo', 'mol2',
                                '-c', charge_type, '-at', atom_type,
                                '-nc', str(net_charge)],
                               **subprocess_kwargs)


if use_agastyas_charges:
    print('Copying REST charges from left_q.mol2 and right_q.mol2')
    # make a copy of the files before modifying them
    shutil.copy(workplace_root / 'left.mol2', workplace_root / 'left_before_RESP.mol2')
    shutil.copy(workplace_root / 'right.mol2', workplace_root / 'right_before_RESP.mol2')
    # take the .mol2 file and correct the charges to reflect Agastya's RESP
    set_charges_from_mol2(workplace_root / 'left.mol2', workplace_root / left_charges, by_atom_name=True)
    set_charges_from_mol2(workplace_root / 'right.mol2', workplace_root / right_charges, by_atom_name=True)

if use_original_coor:
    print(f'Copying coordinates from {left_ligand} and {right_ligand} since antechamber changes them slightly')
    # copy the files before applying the coordinates
    shutil.copy(workplace_root / 'left.mol2', workplace_root / 'left_before_COOR.mol2')
    shutil.copy(workplace_root / 'right.mol2', workplace_root / 'right_before_COOR.mol2')
    set_coor_from_ref(workplace_root / 'left.mol2', workplace_root / left_ligand, by_atom_name=True)
    set_coor_from_ref(workplace_root / 'right.mol2', workplace_root / right_ligand, by_atom_name=True)
    # set_coor_from_ref(workplace_root / 'left.mol2', workplace_root / left_ligand, by_index=True)
    # set_coor_from_ref(workplace_root / 'right.mol2', workplace_root / right_ligand, by_index=True)

# rename the atom names to ensure they are unique across the two molecules
renameAtomNamesUnique(workplace_root / 'left.mol2', workplace_root / 'right.mol2')

# superimpose the two topologies
suptop, mda_l1, mda_l2 = getSuptop(workplace_root / 'left.mol2', workplace_root / 'right.mol2',
                                   align_molecules=align_molecules,
                                   pair_charge_atol=pair_q_atol,
                                   manual_match=manually_matched, force_mismatch=force_mismatch)

# save the superimposition results
left_right_matching_json = workplace_root / 'joint_meta_fep.json'
save_superimposition_results(left_right_matching_json, suptop)
# hybrid .pdb
write_dual_top_pdb(workplace_root / 'morph.pdb', mda_l1, mda_l2, suptop)
# hybrid .mol2
hybrid_mol2 = workplace_root / 'morph.mol2'
write_merged(suptop, hybrid_mol2)

# generate the functional forms
print('Ambertools parmchk2 generating .frcmod')
left_chk2 = subprocess.run([ambertools_bin / 'parmchk2', '-i', 'left.mol2', '-o', 'left.frcmod',
                             '-f', 'mol2', '-s', atom_type], **subprocess_kwargs)
right_chk2 = subprocess.run([ambertools_bin / 'parmchk2', '-i', 'right.mol2', '-o', 'right.frcmod',
                             '-f', 'mol2', '-s', atom_type], **subprocess_kwargs)

# join the .frcmod files
left_frcmod = workplace_root / 'left.frcmod'
right_frcmod = workplace_root / 'right.frcmod'
hybrid_frcmod = workplace_root / 'morph.frcmod'
join_frcmod_files2(left_frcmod, right_frcmod, hybrid_frcmod)

# if the hybrid .frcmod needs new terms between the appearing/disappearing atoms, insert dummy ones
updated_frcmod_content = check_hybrid_frcmod(hybrid_mol2, hybrid_frcmod, amber_forcefield,
                                             ambertools_bin, ambertools_script_dir, cwd=workplace_root)
with open(hybrid_frcmod, 'w') as FOUT:
    FOUT.write(updated_frcmod_content)

shutil.copy(namd_script_dir / "check_namd_outputs.py", workplace_root)
shutil.copy(namd_script_dir / "ddg.py", workplace_root)

# ---------------------

def prepare_inputs(workplace_root, directory='complex', protein=None,
                   hybrid_mol2=None, hybrid_frc=None,
                   left_right_mapping=None,
                   namd_script_loc=None,
                   submit_script=None,
                   scripts_loc=None,
                   tleap_in=None,
                   protein_ff=None,
                   net_charge=None):

    dest_dir = workplace_root / directory
    if not dest_dir.is_dir():
        dest_dir.mkdir()

    # copy the protein complex .pdb
    if protein is not None:
        shutil.copy(workplace_root / protein, dest_dir)

    # copy the hybrid ligand (topology and .frcmod)
    shutil.copy(hybrid_mol2, dest_dir)
    shutil.copy(hybrid_frc, dest_dir)

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
    leap_in_conf = open(ambertools_script_dir / tleap_in).read()
    if Na_num == 0:
        tleap_Na_ions = ''
    elif Na_num > 0:
        tleap_Na_ions = 'addIons sys Na+ %d' % Na_num
    if Cl_num == 0:
        tleap_Cl_ions = ''
    elif Cl_num > 0:
        tleap_Cl_ions = 'addIons sys Cl- %d' % Cl_num
    open(dest_dir / 'leap.in', 'w').write(leap_in_conf.format(protein_ff=protein_ff,
                                                              NaIons=tleap_Na_ions,
                                                              ClIons=tleap_Cl_ions))

    # ambertools tleap: combine ligand+complex, solvate, generate amberparm
    subprocess_kwargs['cwd'] = dest_dir
    subprocess.run([ambertools_bin / 'tleap', '-s', '-f', 'leap.in'], **subprocess_kwargs)
    hybrid_solv = dest_dir / 'sys_solv.pdb' # generated
    # check if the solvation is correct

    # generate the merged .fep file
    complex_solvated_fep = dest_dir / 'sys_solv_fep.pdb'
    correct_fep_tempfactor(left_right_mapping, hybrid_solv, complex_solvated_fep)

    # fixme - check that the protein does not have the same resname?

    # calculate PBC for an octahedron
    solv_oct_boc = extract_PBC_oct_from_tleap_log(dest_dir / "leap.log")

    # prepare NAMD input files for min+eq+prod
    init_namd_file_min(namd_script_loc, dest_dir, "min.namd",
                       structure_name='sys_solv', pbc_box=solv_oct_boc)
    eq_namd_filenames = generate_namd_eq(namd_script_loc / "eq.namd", dest_dir)
    shutil.copy(namd_script_loc / namd_prod, dest_dir / 'prod.namd')

    # generate 4 different constraint .pdb files (it uses B column)
    constraint_files = create_4_constraint_files(hybrid_solv, dest_dir)

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
            prepareFile(os.path.relpath(hybrid_solv, replica_dir), replica_dir / 'sys_solv.pdb')
            prepareFile(os.path.relpath(complex_solvated_fep, replica_dir), replica_dir / 'sys_solv_fep.pdb')
            # copy ambertools-generated topology
            prepareFile(os.path.relpath(dest_dir / "sys_solv.top", replica_dir), replica_dir / "sys_solv.top")
            # copy the .pdb files with constraints in the B column
            for constraint_file in constraint_files:
                prepareFile(os.path.relpath(constraint_file, replica_dir),
                            replica_dir / os.path.basename(constraint_file))

            # copy the NAMD protocol files
            prepareFile(os.path.relpath(dest_dir / 'min.namd', replica_dir), replica_dir / 'min.namd')
            [prepareFile(os.path.relpath(eq, replica_dir), replica_dir / os.path.basename(eq))
                        for eq in eq_namd_filenames]
            prepareFile(os.path.relpath(dest_dir / 'prod.namd', replica_dir), replica_dir / 'prod.namd')

            # copy a submit script
            if submit_script is not None:
                shutil.copy(namd_script_loc / submit_script, replica_dir / 'submit.sh')

    # copy handy scripts to the main directory
    shutil.copy(scripts_loc / "schedule_separately.py", dest_dir)


##########################################################
# ------------------   Ligand ----------------------

prepare_inputs(workplace_root, directory='lig',
               protein=None,
               hybrid_mol2=hybrid_mol2,
               hybrid_frc=hybrid_frcmod,
               left_right_mapping=left_right_matching_json,
               namd_script_loc=namd_script_dir,
               scripts_loc=script_dir,
               tleap_in='leap_ligand.in',
               protein_ff=amber_forcefield,
               net_charge=net_charge)

##########################################################
# ------------------ Complex  --------------

# calculate the charges of the protein (using ambertools)
protein_net_charge = get_protein_net_charge(workplace_root, protein_filename,
                       ambertools_bin, ambertools_script_dir / 'solv_prot.in',
                       subprocess_kwargs, amber_forcefield)

prepare_inputs(workplace_root, directory='complex',
               protein=protein_filename,
               hybrid_mol2=hybrid_mol2,
               hybrid_frc=hybrid_frcmod,
               left_right_mapping=left_right_matching_json,
               namd_script_loc=namd_script_dir,
               scripts_loc=script_dir,
               tleap_in='leap_complex.in',
               protein_ff=amber_forcefield,
               net_charge=net_charge + protein_net_charge)

