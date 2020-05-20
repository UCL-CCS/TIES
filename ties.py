"""
Load two ligands, run the topology superimposer, and then using the results, generate the NAMD input files.
"""
from generator import *
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
net_charge = 0
# force_mismatch_list = [('O2', 'O4'), ('N3', 'N6')] # None

use_agastyas_charges = True
left_charges = 'left_q.mol2'
right_charges = 'right_q.mol2'
# the coordinates change slightly after antechamber, reassign the coordinates to the .mol2
use_original_coor = True

# ambertools forcefiled
amber_forcefield = "leaprc.ff99SBildn" # used by Agastya before, latest is: "leaprc.protein.ff14SB"
atom_type = 'gaff' # fixme - check?
if use_agastyas_charges:
    # ignore the charges,
    charge_type = 'dc'  # 'Delete Charge'
    # use berendsen
    namd_prod = "prod_2017.namd"
else:
    charge_type = 'bcc'  # 'AM1-BCC'
    namd_prod = "prod.namd"
# amber_ligand_ff = "leaprc.gaff" # fixme - check


# ------------------ Software Configuration
# set the working directory to the one where the script is called
workplace_root = Path(os.getcwd())
print('Working Directory: ', workplace_root)
# set the path to the scripts
code_root = Path(os.path.dirname(__file__))

# conf ambertools
ambertools_bin = PurePosixPath("/home/dresio/software/amber18install/bin/")
antechamber_path = ambertools_bin / "antechamber"
parmchk2_path = ambertools_bin / "parmchk2"
tleap_path = ambertools_bin / "tleap"

# directory with scripts/input files
script_dir = code_root / PurePosixPath('scripts')
namd_script_dir = script_dir / 'namd'
ambertools_script_dir = script_dir / 'ambertools'


# ------------------------------------- Start TIES
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
    subprocess.run([antechamber_path, '-i', left_ligand, '-fi', 'pdb',
                                        '-o', 'left.mol2', '-fo', 'mol2',
                                        '-c', charge_type, '-at', atom_type,
                                        '-nc', str(net_charge)],
                               **subprocess_kwargs)

# prepare the .mol2 files with antechamber (ambertools), assign BCC charges if necessary
if not (workplace_root / 'right.mol2').is_file():
    print('Ambertools antechamber stage: converting to .mol2 and generating charges')
    subprocess.run([antechamber_path, '-i', right_ligand, '-fi', 'pdb',
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

# generate the functional forms
print('Ambertools parmchk2 generating .frcmod')
left_chk2 = subprocess.run([str(parmchk2_path), '-i', 'left.mol2', '-o', 'left.frcmod',
                             '-f', 'mol2', '-s', atom_type], **subprocess_kwargs)
right_chk2 = subprocess.run([parmchk2_path, '-i', 'right.mol2', '-o', 'right.frcmod',
                             '-f', 'mol2', '-s', atom_type], **subprocess_kwargs)

# join the .frcmod files
left_frcmod = workplace_root / 'left.frcmod'
right_frcmod = workplace_root / 'right.frcmod'
hybrid_frcmod = workplace_root / 'morph.frcmod'
join_frcmod_files2(left_frcmod, right_frcmod, hybrid_frcmod)

def check_hybrid_frcmod(mol2_file, hybrid_frcmod, tleap_path, atomff_type):
    """
    Previous code: https://github.com/UCL-CCS/BacScratch/blob/master/agastya/ties_hybrid_topology_creator/output.py
    Check that the output library can be used to create a valid amber topology.
    Add missing terms with no force to pass the topology creation.
    Returns the corrected .frcmod content, otherwise throws an exception.
    """
    # prepare files
    tmp_dir = tempfile.mkdtemp()
    shutil.copy(mol2_file, tmp_dir)
    shutil.copy(hybrid_frcmod, tmp_dir)
    copy_hybrid_frcmod = os.path.join(tmp_dir, os.path.basename(hybrid_frcmod))
    test_system_build_script = '''
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff

frcmod = loadamberparams morph.frcmod
hybrid = loadMol2 morph.mol2
saveamberprep hybrid test.prepc
saveamberparm hybrid test.top test.crd
savepdb hybrid test.pdb
quit
'''
    # prepare the input file
    leap_in_filename = 'test.leapin'
    with open(os.path.join(tmp_dir, leap_in_filename), 'w') as FOUT:
        FOUT.write(test_system_build_script)

    # fixme - temporary solution
    output = os.popen(f'cd {tmp_dir} && {tleap_path} -s -f {leap_in_filename}').read()

    missing_angles = []
    missing_dihedrals = []
    for line in output.splitlines():
        if "Could not find angle parameter:" in line:
            cols = line.split(':')
            angle = cols[-1].strip()
            if angle not in missing_angles:
                missing_angles.append(angle)
        elif "No torsion terms for" in line:
            cols = line.split()
            torsion = cols[-1].strip()
            if torsion not in missing_dihedrals:
                missing_dihedrals.append(torsion)

    if missing_angles or missing_dihedrals:
        # raise Exception('hybrid frcmod missing dihedrals or angles')
        print('WARNING: Adding default values for missing dihedral to frcmod')
        with open(copy_hybrid_frcmod) as FRC:
            frcmod_lines = FRC.readlines()

        new_frcmod = open(copy_hybrid_frcmod, 'w')
        for line in frcmod_lines:
            new_frcmod.write(line)
            if 'ANGLE' in line:
                for angle in missing_angles:
                    new_frcmod.write(
                        '{:<14}0     120.010   same as ca-ca-ha\n'.format(angle))
            if 'DIHE' in line:
                for angle in missing_dihedrals:
                    new_frcmod.write(
                        '{:<14}1    0.00       180.000           2.000      same as X -c2-ca-X\n'.format(angle))

        new_frcmod.close()
        # returncode = subprocess.call(['tleap', '-s', '-f', leap_in_filename])
        returncode = os.popen(f'cd {tmp_dir} && {tleap_path} -s -f {leap_in_filename}').read()

        if not "Errors = 0" in returncode:
            print('ERROR: Test of the hybrid topology failed')
            sys.exit(1)

    print('\nHybrid topology created correctly')
    with open(copy_hybrid_frcmod) as FRC:
        frcmod_with_null_terms = FRC.read()
    return frcmod_with_null_terms


# the hybrid .frcmod might not contains "new terms" between the appearing/disappearing atoms.
# In that case, insert fake terms
# fixme - clean up this part: remove the direct call to tleap?
updated_frcmod = check_hybrid_frcmod(hybrid_mol2, hybrid_frcmod, tleap_path, atom_type)
with open(hybrid_frcmod, 'w') as FOUT:
    FOUT.write(updated_frcmod)

# ---------------------
def prepare_inputs(workplace_root, directory='complex', protein=None,
                   hybrid_mol2=None, hybrid_frc=None,
                   left_right_mapping=None,
                   namd_script_loc=None,
                   submit_script=None,
                   scripts_loc=None,
                   tleap_in=None,
                   protein_ff=None):

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
    leap_in_conf = open(ambertools_script_dir / tleap_in).read()
    open(dest_dir / 'leap.in', 'w').write(leap_in_conf.format(protein_ff=protein_ff))

    # ambertools tleap: combine ligand+complex, solvate, generate amberparm
    subprocess_kwargs['cwd'] = dest_dir
    subprocess.run([tleap_path, '-s', '-f', 'leap.in'], **subprocess_kwargs)
    hybrid_solv = dest_dir / 'morph_solv.pdb' # generated

    # generate the merged .fep file
    complex_solvated_fep = dest_dir / 'morph_solv_fep.pdb'
    correct_fep_tempfactor(left_right_mapping, hybrid_solv, complex_solvated_fep)

    # fixme - check that the protein does not have the same resname?

    # calculate PBC for an octahedron
    solv_oct_boc = extract_PBC_oct_from_tleap_log(dest_dir / "leap.log")

    # prepare NAMD input files for min+eq+prod
    init_namd_file_min(namd_script_loc, dest_dir, "min.namd",
                       structure_name='morph_solv', pbc_box=solv_oct_boc)
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
            shutil.copy(hybrid_solv, replica_dir)
            shutil.copy(complex_solvated_fep, replica_dir)
            # copy ambertools-generated topology
            shutil.copy(dest_dir / "morph_solv.top", replica_dir)
            # copy the .pdb files with constraints in the B column
            [shutil.copy(constraint_file, replica_dir) for constraint_file in constraint_files]

            # copy the NAMD protocol files
            shutil.copy(dest_dir / "min.namd", replica_dir)
            [shutil.copy(eq, replica_dir) for eq in eq_namd_filenames]
            shutil.copy(dest_dir / 'prod.namd', replica_dir / 'prod.namd')

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
               tleap_in='leap_ligand.in',
               protein_ff=amber_forcefield)

##########################################################
# ------------------ complex-complex --------------

prepare_inputs(workplace_root, directory='complex',
               protein=protein_filename,
               hybrid_mol2=hybrid_mol2,
               hybrid_frc=hybrid_frcmod,
               left_right_mapping=left_right_matching_json,
               namd_script_loc=namd_script_dir,
               scripts_loc=script_dir,
               tleap_in='leap_complex.in',
               protein_ff=amber_forcefield)

