# verify that the min_api.py worked fine
import glob
import pathlib

# the default directory was created
tiesdir = pathlib.Path('ties')
assert tiesdir.exists()

# the files have been saved
assert (tiesdir / 'fep_l02_l03.json').exists()
assert (tiesdir / 'l02_l03_morph.mol2').exists()
assert (tiesdir / 'l02_l03_morph.pdb').exists()
