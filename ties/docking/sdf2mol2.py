"""
Convert SDF to .MOL2 file.

Extract the charges and the atom types from the properties.
"""

from pathlib import Path
import warnings
import argparse

from ties import Ligand


def sdf_to_mol2(filename: Path, resname="MOL"):
    warnings.warn("Reading only 1 frame from the SDF")
    ligand = Ligand(filename)
    ligand.pmd_structure.save(str(filename.parent / f"{filename.stem}.mol2"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-sdf",
        metavar="str",
        dest="filename",
        type=Path,
        required=True,
        help="An SDF file",
    )
    args = parser.parse_args()

    sdf_to_mol2(args.filename)
