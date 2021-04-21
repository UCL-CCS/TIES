from ties import Pair
from ties import Config

config = Config()
config.ambertools_home = '/home/dresio/software/amber18install'
config.ligand_net_charge = -1

pair = Pair('MCL1_lig_02.mol2', 'MCL1_lig_03.mol2', config=config)
pair.make_atom_names_unique()
hybrid = pair.superimpose()

# save meta data
hybrid.write_summary_json()
hybrid.write_pdb()
hybrid.write_mol2()

hybrid.prepare_inputs()

print('done')