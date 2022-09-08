from ties import Pair, Config, Protein, MD

#Settings for simulation
config = Config()
config.workdir = 'ties20'
config.md_engine = 'openmm'
config.protein = 'protein.pdb'
config.protein_ff = 'leaprc.protein.ff14SB'

#load the two ligands and create a hybrid
pair = Pair('l02.mol2', 'l03.mol2', ligand_net_charge=-1, config=config)
pair.make_atom_names_unique()
hybrid = pair.superimpose()

#setup ligand simulation
hybrid.prepare_inputs()

#add protien and setup complex simulation
protein = Protein(config.protein, config)
hybrid.prepare_inputs(protein=protein)

#run ligand and complex simulations
legs = ['lig', 'com']
for leg in legs:
   md = MD('./ties20/ties-l02-l03/{}'.format(leg), fast=True)
   md.setup()
   #md.run()

#run the analysis of these simulations
exp_data = {'ties20': {'ties-l02-l03': [0.0, 0.0]}}
md.analysis(exp_data, legs)
