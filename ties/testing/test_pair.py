"""
These tests focus on the Ligand
"""
import parmed

from ties import Ligand, Pair


def test_atom_names_uniqe():
    # the atoms names across these two files overlaps
    # test renaming them to avoid this issue
    pair = Pair('data/ligq.mol2', 'data/l_H18q.mol2')
    pair.make_atom_names_unique()

    # it's convoluted, but we have to open separately the newly created file
    ligA = parmed.load_file(str(pair.current_ligA))
    ligZ = parmed.load_file(str(pair.current_ligZ))

    common_atom_names = {a.name for a in ligA.atoms}.intersection({a.name for a in ligZ.atoms})
    assert len(common_atom_names) == 0