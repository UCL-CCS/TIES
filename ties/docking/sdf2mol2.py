"""
Convert SDF to .MOL2 file.

Extract the charges and the atom types from the properties.
"""

from pathlib import Path
import sys
import warnings

from rdkit import Chem
import parmed as pmd


def sdf_to_mol2(filename: Path, resname="MOL"):
    # load the SDF file with parmed
    warnings.warn("Reading only 1 frame from the SDF")
    rd_mol = Chem.SDMolSupplier(str(filename), removeHs=False)[0]

    # extract the props
    bcc_types = rd_mol.GetProp("atom.dprop.GAFFAtomType").split()
    partial_charges = map(float, rd_mol.GetProp("atom.dprop.PartialCharge").split())

    pmd_mol = pmd.load_rdkit(rd_mol)
    pmd_mol.residues[0].name = resname

    for atom, bcc_type, partial_charge in zip(
        pmd_mol.atoms, bcc_types, partial_charges
    ):
        atom.charge = partial_charge
        atom.type = bcc_type

    pmd_mol.save(str(filename.parent / f"{filename.stem}.mol2"))


if __name__ == "__main__":
    if len(sys.argv) == 2:
        sdf = Path(sys.argv[1])
        assert sdf.exists()

        sdf_to_mol2(sdf)

    else:
        raise Exception("forgot filename?")
