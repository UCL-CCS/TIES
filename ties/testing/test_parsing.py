from rdkit import Chem
from rdkit.Chem.AllChem import EmbedMolecule

from ties.parsing import pmd_structure_from_rdmol


def test_mol_from_rdmol():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    EmbedMolecule(mol)

    pmd_structure = pmd_structure_from_rdmol(mol)

    assert len(pmd_structure.bonds) == mol.GetNumBonds()
    assert len(pmd_structure.atoms) == mol.GetNumAtoms()
