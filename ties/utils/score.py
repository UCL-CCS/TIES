"""
Score conformers in the binding pocket using MM potential energy.
"""

import os
import time
from ast import literal_eval

from rdkit import Chem
import prody
import fegrow


def extract_best_conformer(
    mol: Chem.Mol, prody_protein, protein_pdb="fe_rec_final.pdb"
):
    """

    :param mol_mol2: Should contain all the conformers for evaluation
    :param prody_protein:
    :param protein_pdb:
    :return:
    """
    start = time.time()

    # convert to FEgrow RMol
    rmol = fegrow.RMol(mol)
    rmol.remove_clashing_confs(prody_protein, min_dst_allowed=0.5)

    if rmol.GetNumConformers() == 0:
        print("Warning: no conformers for energy minimisation: ", rmol.GetProp("_Name"))
        return

    # extract the mcs mapping
    # we will freeze the MCS area in the ligand before optimisation

    mcs = literal_eval(rmol.GetProp("MCS(ref,lig)"))
    # extract
    ligand_atoms_to_freeze = [idx for _, idx in mcs]

    energies = rmol.optimise_in_receptor(
        receptor_file=protein_pdb,
        ligand_force_field="openff",
        use_ani=False,
        sigma_scale_factor=0.8,
        relative_permittivity=4,
        water_model=None,
        platform_name="CPU",  # or e.g. 'CUDA'
        ligand_freeze=ligand_atoms_to_freeze,
    )

    # save the lowest energy conformer
    best = energies[energies.Energy == energies.Energy.min()]
    best_conf_id = int(best.index.values[0])
    # print(f"Optimised, {time.time() - start:.2f}s")

    # keep the lowest energy conformer
    for conf in list(rmol.GetConformers()):
        if conf.GetId() != best_conf_id:
            rmol.RemoveConformer(conf.GetId())

    return rmol


def prepare_protein(name, parsed_pdb="fe_rec_final.pdb"):
    # prepare the protein
    # load the complex with the ligand
    if not os.path.exists(parsed_pdb):
        print("Preparing the protein")

        sys = prody.parsePDB(name)
        # remove any unwanted molecules
        rec = sys.select("not (nucleic or hetatm or water)")

        # save the processed protein
        prody.writePDB("fe_rectmp.pdb", rec)

        # fix the receptor file (missing residues, protonation, etc)
        fegrow.fix_receptor("fe_rectmp.pdb", parsed_pdb)

    protein = prody.parsePDB(parsed_pdb)
    return protein, parsed_pdb


if __name__ == "__main__":
    # minimise the energies for each conformer
    # and save the most energetically favourable
    prody_protein, parsed_pdb = prepare_protein("protein.pdb")
    extract_best_conformer("mols/74.mol2", prody_protein, parsed_pdb)
