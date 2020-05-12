"""
Load two ligands, run the topology superimposer, and then
using the results, generate the NAMD input files.

todo - use ambertools python interface rather than ambertools directly,
       or find some other API bbblocks building blocks? or some other API,
       or create minimal your own
todo - adjust the input for surfsara for nwo
todo - estimate the size and legnth of the simulation, to adjust the number of cores/nodes
todo -

frcmod file format: `http://ambermd.org/FileFormats.php#frcmod`

Improvements Long Term:
 - adapt pathlib to work with the directories

fixme - create a function that takes care of each thing
fixme - liglig and protprot follow pretty much the same protocol, so it should be one function?

"""
import os
import numpy as np
import shutil
import sys
import subprocess
from pathlib import Path, PurePosixPath
from generator import *

# fixme - turn into a function and give it the hpc submit
hpc_submit = 'hpc_surfsara.sh'
workplace_root = Path('/home/dresio/code/BAC2020/namd_study/validation/ethane')

def prepare_antechamber_parmchk2(source_script, target_script, net_charge):
    """
    Prepare the ambertools scripts.
    Particularly, create the scritp so that it has the net charge
    # fixme - run antechamber directly with the right settings from here?
    # fixme - check if antechamber has a python interface?
    """
    net_charge_set = open(source_script).read().format(net_charge=net_charge)
    open(target_script, 'w').write(net_charge_set)

# take care of the ligand-ligand without the protein
liglig_workplace = workplace_root / 'lig'
if not liglig_workplace.is_dir():
    liglig_workplace.mkdir()

morph_solv = workplace_root / "morph_solv.pdb"
morph_solv_fep = workplace_root / "morph_solv_fep.pdb"

# generate the constraint files in liglig
constraint_files = create_4_constraint_files(morph_solv, liglig_workplace)
pbc_dimensions = get_PBC_coords(morph_solv)
init_namd_file_min(workplace_root, liglig_workplace, "min.namd",
                           structure_name='morph_solv', pbc_box=pbc_dimensions)
# prepare the namd eq
eq_namd_files = generate_namd_eq(workplace_root / "eq.namd", liglig_workplace)
# namd production protocol
prod_filename = init_namd_file_prod(workplace_root, liglig_workplace, "prod.namd", structure_name='morph_solv')

shutil.copy(workplace_root / "node_schedule_separately.py", liglig_workplace)

# Generate the directory structure for all the lambdas, and copy the files
for lambda_step in [0, 0.05] + list(np.linspace(0.1, 0.9, 9)) + [0.95, 1]:
    lambda_path = liglig_workplace / f'lambda_{lambda_step:.2f}'
    if not os.path.exists(lambda_path):
        os.makedirs(lambda_path)

    shutil.copy(workplace_root / "node_hpc_surfsara.sh", lambda_path / 'submit.sh')

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
        shutil.copy(workplace_root / "morph_solv.psf", replica_dir)
        # the normal coordinates
        shutil.copy(workplace_root / "morph_solv.pdb", replica_dir)

        # copy the restraint files, which use b column
        [shutil.copy(constraint_file, replica_dir) for constraint_file in constraint_files]

        # copy NAMD protocol
        shutil.copy(liglig_workplace / 'min.namd', replica_dir)
        [shutil.copy(eq, replica_dir) for eq in eq_namd_files]
        shutil.copy(prod_filename, replica_dir)

        shutil.copy(workplace_root / "par_all22_prot.inp", replica_dir)

        # todo - create/copy the 4 different EQ files
        # copy the submit script
        shutil.copy(workplace_root / hpc_submit, replica_dir / 'submit.sh')

# copy the scheduler to the main directory
shutil.copy(workplace_root / "schedule_separately.py", liglig_workplace)
shutil.copy(workplace_root /"check_namd_outputs.py", liglig_workplace)

print('THIS IS A VALIDATION STUDY WITHOUT A COMPLEX, TERMINATING!!!!!')
sys.exit(0)
