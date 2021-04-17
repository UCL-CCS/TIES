from ties import Morph

morph = Morph('MCL1_lig_02.mol2', 'MCL1_lig_03.mol2')
suptop = morph.compute_suptop()

# save meta data
morph.write_summary_json()
morph.write_pdb()
morph.write_hybrid_mol2()