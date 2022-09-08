# verify that the min_api.py worked fine
import glob
import pathlib

# the default directory was created
tiesdir = pathlib.Path('cmdties_protein')
assert tiesdir.exists()

# the files have been saved
assert (tiesdir / 'meta_l02_l03.json').exists()
assert (tiesdir / 'l02_l03_morph.mol2').exists()
assert (tiesdir / 'l02_l03_morph.pdb').exists()
assert (tiesdir / 'ties-l02-l03' / 'com' / 'build' / 'sys_solv.top' ).exists()
