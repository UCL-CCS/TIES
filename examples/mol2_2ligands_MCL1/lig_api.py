from ties import Ligand


lig = Ligand('l02_same_atom_name.mol2')

lig.make_atom_names_correct()
assert lig.atom_names_correct()

# prepare the .mol2 input
lig.antechamber_prepare_mol2()

# the final .mol2 file
assert lig.current.exists()
#
print(lig.renaming_map)