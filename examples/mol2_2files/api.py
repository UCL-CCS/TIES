from ties import Morph

morph = Morph('MCL1_lig_02.mol2', 'MCL1_lig_03.mol2')
morph.make_atom_names_unique()
hybrid = morph.compute_suptop()

# save meta data
hybrid.write_summary_json()
hybrid.write_pdb()
hybrid.write_mol2()