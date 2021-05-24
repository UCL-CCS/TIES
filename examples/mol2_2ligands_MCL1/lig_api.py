import tempfile

from ties import Ligand
from ties import Config

# we can avoid saving any of the files by using the temporary directory
td = tempfile.TemporaryDirectory()

c = Config()
# specify where the work should take place
c.workdir = td.name

lig = Ligand('l02_same_atom_name.mol2', c)

lig.make_atom_names_correct()
assert lig.atom_names_correct()

# what about?
lig.antechamber_prepare_mol2()

# extract the final results
print(lig.current)
# print the corrections introduced
print(lig.renaming_map)