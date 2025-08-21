import os
from ties import Pair
from ties import Config

# explicitly create config (which will be used by all classes underneath)
config = Config()
config.ligand_net_charge = -1

pair = Pair("l02.mol2", "l03.mol2", config=config)
pair.make_atom_names_unique()

# overwrite the previous config settings with relevant parameters
hybrid = pair.superimpose(
    use_element_in_superimposition=True, redistribute_q_over_unmatched=True
)

# prep for the output
os.mkdir("explicit") if not os.path.exists else None

# save meta data to specific locations
hybrid.write_metadata("explicit/result.json")
hybrid.write_pdb("explicit/result.pdb")
hybrid.write_mol2("explicit/result.mol2")

hybrid.prepare_inputs()
