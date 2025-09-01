"""
Convert MOL2 to SDF file.

Extract the charges and the atom types and save them in the properties.
"""

from pathlib import Path
import sys
import re

from rdkit import Chem
import parmed as pmd


def mol2_to_sdf(filename: Path, resname="MOL"):
    # load the .mol2 file with parmed
    mol2 = pmd.load_file(str(filename), structure=True)

    # parmed appears to occasionally add a bond(?),
    # so verify that parmed is not adding any bonds or atoms
    re_sults = re.search(
        r"^\s*(?P<atom_n>\d+)\s*(?P<bond_n>\d+)\s*\d+\s*\d+\s*\d+\s*$",
        open(filename).read(),
        re.MULTILINE,
    )
    assert len(mol2.atoms) == int(re_sults.group("atom_n"))
    assert len(mol2.bonds) == int(re_sults.group("bond_n"))

    # parmed reads some atoms wrong, e.g. Cl1 as C,
    # correct it
    for atom in mol2.atoms:
        if atom.type == "cl":
            atom.atomic_number = 17
        elif atom.type == "br":
            atom.atomic_number = 35

    # convert
    rdmol = mol2.rdkit_mol

    # validate: check if the atoms are in the same order
    for rdatom, pmdatom in zip(rdmol.GetAtoms(), mol2.atoms):
        assert rdatom.GetAtomicNum() == pmdatom.atomic_number

    # copy pmd bond order to rdmol
    for rd_bond, pmd_bond in zip(rdmol.GetBonds(), mol2.bonds):
        # verify the bonds are the same
        assert rd_bond.GetBeginAtomIdx() == pmd_bond.atom1.idx
        assert rd_bond.GetEndAtomIdx() == pmd_bond.atom2.idx

        # see https://parmed.github.io/ParmEd/html/topobj/parmed.topologyobjects.Bond.html
        if pmd_bond.order == 1:
            rd_bond.SetBondType(Chem.BondType.SINGLE)
        elif pmd_bond.order == 2:
            rd_bond.SetBondType(Chem.BondType.DOUBLE)
        elif pmd_bond.order == 3:
            rd_bond.SetBondType(Chem.BondType.TRIPLE)
        elif pmd_bond.order == 1.5:
            rd_bond.SetBondType(Chem.BondType.AROMATIC)
        else:
            raise NotImplementedError("Missing bonds?")

    # extract the props
    rdmol.SetProp("atom.dprop.GAFFAtomType", str([a.type for a in mol2.atoms]))
    rdmol.SetProp(
        "atom.dprop.PartialCharge", " ".join(str(a.charge) for a in mol2.atoms)
    )

    rdmol.SetProp("_Name", filename.stem)

    with Chem.SDWriter(filename.parent / f"{filename.stem}.sdf") as SD:
        SD.write(rdmol)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        mol2 = Path(sys.argv[1])
        assert mol2.exists()

        mol2_to_sdf(mol2)

    else:
        raise Exception("forgot filename?")
