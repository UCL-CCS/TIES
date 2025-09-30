import argparse
import warnings
from pathlib import Path
from openff.toolkit import Molecule

from ties.modules import fit, score
from ties.modules.utils.utils import write_mol, paths_from_glob


def docker(
    pairs: tuple[str, str, float],
    mols: dict[str, Molecule],
    crystal_structure_label="0",
    protein_filename: str = "fe_rec_final.pdb",
    confs_dir: Path = Path("confs"),
    results_dir: Path = Path("results"),
):
    # set the reference structure
    fitted = {crystal_structure_label: mols[crystal_structure_label]}

    prody_protein, parsed_pdb = score.prepare_protein(protein_filename)

    for node_ref, node, dst in pairs:
        if node in fitted:
            warnings.warn("node already in fitted")
            continue

        try:
            node_to_ref = fit.mcs_fit(
                fitted[node_ref].to_rdkit(), mols[node].to_rdkit(), conformers=1000
            )

            # save these
            write_mol(node_to_ref, confs_dir / f"{node_ref}-{node}.sdf")

            # score
            bestconf = score.extract_best_conformer(
                node_to_ref, prody_protein, parsed_pdb
            )

            off_bestconf = Molecule.from_rdkit(bestconf, allow_undefined_stereo=True)

            # update fitted
            fitted[node] = off_bestconf

            off_bestconf.to_file(
                results_dir / f"{node_ref}-{node}.sdf", file_format="sdf"
            )
        except Exception as E:
            print("failed", node_ref, node, E)
            # raise Exception from E


parser = argparse.ArgumentParser()
parser.add_argument(
    "-pairs",
    metavar="TXT_FILE_PAIRS",
    dest="pairs_filename",
    type=Path,
    required=True,
    help="A file with pairs (3 columns: node1, node2, dst) ",
)
parser.add_argument(
    "-sdfs",
    metavar="FILENAMES_GLOB",
    dest="mols_filenames",
    type=paths_from_glob,
    required=False,
    default="mols/*.sdf",
    help="SDF filenames with molecules to work on (after prep)",
)
parser.add_argument(
    "-reference",
    metavar="LIGAND_NAME",
    dest="reference_ligand_label",
    type=Path,
    required=False,
    default="0",
    help="The label of the reference ligand for docking",
)
parser.add_argument(
    "-prot",
    metavar="PDB_FILENAME",
    dest="protein_filename",
    type=Path,
    required=False,
    default="fe_rec_final.pdb",
    help="The protein filename",
)
parser.add_argument(
    "-confs_dir",
    metavar="OUT_DIR",
    dest="confs_dir",
    type=Path,
    required=False,
    default="confs",
    help="A directory for the generated initial guesses. ",
)
parser.add_argument(
    "-results_dir",
    metavar="OUT_DIR",
    dest="results_dir",
    type=Path,
    required=False,
    default="results",
    help="A directory for the final docked structures. ",
)

if __name__ == "__main__":
    args = parser.parse_args()

    # get the pairs and the distances
    pairs = [line.split() for line in open(args.pairs_filename).readlines()]

    mols = [
        Molecule.from_file(f, allow_undefined_stereo=True) for f in args.mols_filenames
    ]
    mols = {m.name: m for m in mols}

    docker(
        pairs,
        mols,
        crystal_structure_label=args.reference_ligand_label,
        protein_filename=args.protein_filename,
        confs_dir=args.confs_dir,
        results_dir=args.results_dir,
    )
