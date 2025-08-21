# verify that the min_api.py worked fine
from pathlib import Path

# the files have been saved
assert Path("meta_l02_l03.json").exists()
assert Path("l02_l03_morph.pdb").exists()
assert Path("l02_l03_morph.mol2").exists()
