"""
These tests focus on the Ligand
"""

import rdkit.Chem

from fegrow.testing.library_gen import smiles
from ties.ligand import Ligand


def test_correct_atom_names_HAY(data):
    # atom names do not follow the right format
    lig = Ligand(data / "l_HAY.mol2")
    assert lig.are_atom_names_correct()


def test_correct_atom_names_HAY_renamed(data):
    # atom names do not follow the right format
    lig = Ligand(data / "l_HAY.mol2")
    lig.correct_atom_names()
    assert lig.are_atom_names_correct()


def test_ligand_from_rdkit_mol():
    rdmol = rdkit.Chem.MolFromSmiles("COO")
    lig = Ligand(rdmol)
    # check if the file was saved into an SDF
    assert lig.current.exists()


def test_ligand_from_to_rdkit():
    rd_mol = rdkit.Chem.AddHs(rdkit.Chem.MolFromSmiles("COO"))
    ligand = Ligand(rd_mol)

    # convert it back
    rd_mol_back = ligand.to_rdkit()

    assert rd_mol.GetNumAtoms() == rd_mol_back.GetNumAtoms()
    assert rd_mol.GetNumBonds() == rd_mol_back.GetNumBonds()

    # check if the atoms are in the same order
    for ref_atom, recreated_atom in zip(rd_mol.GetAtoms(), rd_mol_back.GetAtoms()):
        assert ref_atom.GetAtomicNum() == recreated_atom.GetAtomicNum()


def test_rdkit_losing_chiral_data_parmed_conversion():
    """
    When using ParmEd to convert a molecule to RDKit, the chiral information is lost.

    We force ParmEd conversion by removing the original RDKit molecule.
    """
    smiles_with_chiral_centre = "F[C@](Cl)(Br)I"
    m1 = rdkit.Chem.MolFromSmiles(smiles_with_chiral_centre)  # one enantiomer
    l1 = Ligand(m1)
    del l1.original_rdmol
    m1_recovered = rdkit.Chem.MolToSmiles(l1.to_rdkit())
    assert smiles_with_chiral_centre != m1_recovered


def test_rdkit_keeps_chiral_data_via_original_rdkit_ligand():
    """
    The Ligand by default keeps the original RDKit molecule.
    """
    smiles_with_chiral_centre = "F[C@](Cl)(Br)I"
    m1 = rdkit.Chem.MolFromSmiles(smiles_with_chiral_centre)  # one enantiomer
    l1 = Ligand(m1)
    m1_recovered = rdkit.Chem.MolToSmiles(l1.to_rdkit())
    assert smiles_with_chiral_centre == m1_recovered
