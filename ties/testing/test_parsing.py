import tempfile

import parmed
from rdkit import Chem
from rdkit.Chem.AllChem import EmbedMolecule

from ties.bb.gaff_atom_types import gaff_element_map as gaff_atom_type_map
from ties.parsing import pmd_structure_from_rdmol, correct_mol2_gaff_atoms


def test_mol_from_rdmol():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    EmbedMolecule(mol)

    pmd_structure = pmd_structure_from_rdmol(mol)

    assert len(pmd_structure.bonds) == mol.GetNumBonds()
    assert len(pmd_structure.atoms) == mol.GetNumAtoms()


def test_mol2_fileformat_GAFF_atom_types_read_with_parmed():
    """
    Parmed reads the MOL2 GAFF atom types as the wrong elements.
    Check all GAFF elements
    """
    mol_mol2_minimal_template = (
        "@<TRIPOS>MOLECULE\nUNL\n 1 1 1 0 0\nSMALL\nbcc\n\n\n"
        "@<TRIPOS>ATOM\n 1 C1 0 0 0 {atom_type} 1 UNL 0\n"
    )

    for atom_type, element in gaff_atom_type_map.items():
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".mol2") as TF:
            TF.write(mol_mol2_minimal_template.format(atom_type=atom_type))
            TF.flush()
            pmd_mol = parmed.load_file(TF.name)

        if element != pmd_mol.atoms[0].element_name:
            correct_mol2_gaff_atoms(pmd_mol)

            # correction worked
            assert element == pmd_mol.atoms[0].element_name
