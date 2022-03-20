from ties import Pair
from ties import Config
from ties import Protein

config = Config()
config.workdir = 'ties20'
config.md_engine = 'openmm'
pair = Pair('l02.mol2', 'l03.mol2', ligand_net_charge=-1, config=config)
pair.make_atom_names_unique()

hybrid = pair.superimpose()

# save meta data
hybrid.write_metadata()
hybrid.write_pdb()
hybrid.write_mol2()

hybrid.prepare_inputs()

config.protein = 'protein.pdb'
config.protein_ff = 'leaprc.protein.ff14SB'
protein = Protein(config.protein, config)
hybrid.prepare_inputs(protein=protein)
