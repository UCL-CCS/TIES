from ties import Pair

pair = Pair("l02.mol2", "l03.mol2")
hybrid = pair.superimpose()

# save the results
hybrid.write_metadata("meta_l02_l03.json")
hybrid.write_pdb("l02_l03_morph.pdb")
hybrid.write_mol2("l02_l03_morph.mol2")
