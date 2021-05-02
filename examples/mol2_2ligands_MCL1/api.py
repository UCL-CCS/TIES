from ties import Pair
from ties import Config

config = Config()
config.ligand_net_charge = -1

pair = Pair('l02.mol2', 'l03.mol2', config=config)
pair.make_atom_names_unique()
hybrid = pair.superimpose()

# save meta data
hybrid.write_summary_json()
hybrid.write_pdb()
hybrid.write_mol2()

hybrid.prepare_inputs()
