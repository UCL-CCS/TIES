from ties import Pair

pair = Pair('MCL1_lig_02.mol2', 'MCL1_lig_03.mol2')
pair.make_atom_names_unique()
hybrid = pair.superimpose()

# save meta data
hybrid.write_summary_json()
hybrid.write_pdb()
hybrid.write_mol2()