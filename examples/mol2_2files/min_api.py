from ties import Pair

pair = Pair('MCL1_lig_02.mol2', 'MCL1_lig_03.mol2')
hybrid = pair.superimpose()

# save the results
hybrid.write_summary_json()
hybrid.write_pdb()
hybrid.write_mol2()