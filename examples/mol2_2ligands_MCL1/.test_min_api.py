# verify that the min_api.py worked fine
import glob
import pathlib

# the default directory was created
tiesdir = pathlib.Path('ties')
assert tiesdir.exists()

# the files have been saved
assert (tiesdir / 'fep_MCL1_lig_02_MCL1_lig_03.json').exists()
assert (tiesdir / 'MCL1_lig_02_MCL1_lig_03_morph.mol2').exists()
assert (tiesdir / 'MCL1_lig_02_MCL1_lig_03_morph.pdb').exists()
