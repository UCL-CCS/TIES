"""
Score conformers in the binding pocket using MM potential energy.
"""

import argparse
import os
from ast import literal_eval
from pathlib import Path

from rdkit import Chem
import prody
import fegrow
from openff.toolkit import Molecule

from ties.modules.utils.utils import paths_from_glob, load_conformers


def extract_best_conformer(
    mol: Chem.Mol, prody_protein, protein_pdb="fe_rec_final.pdb"
):
    """

    :param mol_mol2: Should contain all the conformers for evaluation
    :param prody_protein:
    :param protein_pdb:
    :return:
    """
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
        ligand_indices_to_freeze=ligand_atoms_to_freeze,
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

    # if __name__ == "__main__":
    # minimise the energies for each conformer
    # and save the most energetically favourable


parser = argparse.ArgumentParser(
    description="Extract the best conformer for the FEP calculations. "
    "Minimise the energy of all conformers, and extract the "
    "lowest energy conformer. ",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-sdfs",
    metavar="str",
    dest="filenames",
    type=paths_from_glob,
    required=False,
    default="fitted/*.sdf",
    help="An SDF file with molecules",
)
parser.add_argument(
    "-prot",
    metavar="str",
    dest="protein",
    type=Path,
    required=True,
    help="The protein for scoring with FEgrow (see FEgrow on how to prepare the PDB)",
)
parser.add_argument(
    "-dir",
    metavar="str",
    dest="out_dir",
    type=Path,
    required=False,
    default=Path("results"),
    help="Directory to which the output should be saved",
)

if __name__ == "__main__":
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(exist_ok=True)

    prody_protein, parsed_pdb = prepare_protein(args.protein)

    for filename in args.filenames:
        print("Next", filename)

        fitconfs = load_conformers(filename)

        filename_result = out_dir / f"{filename.stem}.sdf"
        if filename_result.exists():
            continue

        try:
            bestconf = extract_best_conformer(fitconfs, prody_protein, parsed_pdb)
            Molecule.from_rdkit(bestconf).to_file(filename_result, file_format="sdf")
        except Exception as E:
            print("failed", filename, E)
