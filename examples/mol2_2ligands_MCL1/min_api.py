from ties import Pair

pair = Pair('l02.mol2', 'l03.mol2')
hybrid = pair.superimpose()

# save the results
hybrid.write_summary_json()
hybrid.write_pdb()
hybrid.write_mol2()