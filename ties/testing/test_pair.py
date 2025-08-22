"""
These tests focus on the Ligand
"""

import parmed
import rdkit.Chem

from ties import Pair


def test_atom_names_uniqe(data):
    # the atoms names across these two files overlaps
    # test renaming them to avoid this issue
    pair = Pair(data / "ligq.mol2", data / "l_H18q.mol2")
    pair.make_atom_names_unique()

    # it's convoluted, but we have to open separately the newly created file
    ligA = parmed.load_file(str(pair.current_ligA))
    ligZ = parmed.load_file(str(pair.current_ligZ))

    common_atom_names = {a.name for a in ligA.atoms}.intersection(
        {a.name for a in ligZ.atoms}
    )
    assert len(common_atom_names) == 0


def test_rdkit_mols():
    cco = rdkit.Chem.MolFromSmiles("CCO")
    ccco = rdkit.Chem.MolFromSmiles("CCCO")
    pair = Pair(cco, ccco)
    suptop = pair.superimpose()
    assert len(suptop) == 3


def test_rdkit_mols_atom_types():
    """
    modify the atom types using the property "BCCAtomTypes"
    """
    cco = rdkit.Chem.MolFromSmiles("CCO")
    ccco = rdkit.Chem.MolFromSmiles("CCCO")

    # only the first two atom types match
    cco.SetProp("BCCAtomTypes", "['c', 'c2', 'o']")
    ccco.SetProp("BCCAtomTypes", "['c', 'c2', 'c3', 'o']")

    pair = Pair(cco, ccco)
    suptop = pair.superimpose(superimposition_starting_heuristic=1.0)
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
    suptop = pair.superimpose()
    assert len(suptop) == 2
