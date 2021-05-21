from ties import Pair

pair = Pair('l02.mol2', 'l03.mol2', ligand_net_charge=-1)
pair.make_atom_names_correct()

hybrid = pair.superimpose()

# save meta data
hybrid.write_metadata()
hybrid.write_pdb()
hybrid.write_mol2()

hybrid.prepare_inputs()
