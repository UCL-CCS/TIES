"""
These tests focus on the Ligand
"""

import rdkit.Chem

from ties import Pair


def test_rdkit_mols():
    cco = rdkit.Chem.MolFromSmiles("CCO")
    ccco = rdkit.Chem.MolFromSmiles("CCCO")
    pair = Pair(cco, ccco)
    suptop = pair.superimpose(use_rdkit_mcs=False)
    assert len(suptop) == 3


def test_rdkit_mols_rdkit_mcs():
    cco = rdkit.Chem.MolFromSmiles("CCO")
    ccco = rdkit.Chem.MolFromSmiles("CCCO")
    pair = Pair(cco, ccco)
    suptop = pair.superimpose(use_rdkit_mcs=True)
    assert len(suptop) == 3


def test_rdkit_mols_rdkit_mcs_enantiomers():
    m1 = rdkit.Chem.MolFromSmiles("F[C@](Cl)(Br)I")  # one enantiomer
    m2 = rdkit.Chem.MolFromSmiles("F[C@@](Cl)(Br)I")  # opposite
    pair = Pair(m1, m2)
    suptop = pair.superimpose(use_rdkit_mcs=True)
    assert len(suptop) != 5


def test_rdkit_mols_atom_types():
    """
    modify the atom types using the property "BCCAtomTypes"
    """
    cco = rdkit.Chem.MolFromSmiles("CCO")
    ccco = rdkit.Chem.MolFromSmiles("CCCO")

    # only the first two atom types match
    cco.SetProp("atom.dprop.GAFFAtomType", "c c2 o")
    ccco.SetProp("atom.dprop.GAFFAtomType", "c c2 c3 o")

    pair = Pair(cco, ccco)
    suptop = pair.superimpose(use_element_in_superimposition=False)
    assert len(suptop) == 2


def test_rdkit_mols_partial_charges():
    """
    set the atom's partial charges using the property "atom.dprop.PartialCharge"
    """
    cco = rdkit.Chem.MolFromSmiles("CCCO")
    ccco = rdkit.Chem.MolFromSmiles("CCCO")

    # only two atoms match using the default filter
    cco.SetProp("atom.dprop.PartialCharge", "0 0 0.1 -0.1")
    ccco.SetProp("atom.dprop.PartialCharge", "0 0 -0.1 0.1")

    pair = Pair(cco, ccco)
    suptop = pair.superimpose(use_element_in_superimposition=False)
    assert len(suptop) == 2
