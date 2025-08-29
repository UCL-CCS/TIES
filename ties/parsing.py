import itertools
import logging
import warnings
import ast

import parmed
import rdkit
import rdkit.Chem

from ties.bb.atom import Atom


logger = logging.getLogger(__name__)
periodic_table = rdkit.Chem.GetPeriodicTable()


def pmd_structure_from_rdmol(mol: rdkit.Chem.Mol):
    """
    Generate a parmed structure from an RDKit Mol.

    The atom types and charges are extracted from the properties.

    :param mol:
    :return:
    """
    parmed_structure = parmed.load_rdkit(mol)

    # verify that they are the same structures
    for rd_a, pmd_a in zip(mol.GetAtoms(), parmed_structure.atoms):
        assert rd_a.GetAtomicNum() == pmd_a.atomic_number

    # extract the charges and the atom types
    pq_prop_openff = "atom.dprop.PartialCharge"
    if mol.HasProp(pq_prop_openff):
        pqs = list(map(float, mol.GetProp(pq_prop_openff).split()))
        assert len(pqs) == mol.GetNumAtoms()
        for atom, pq in zip(parmed_structure.atoms, pqs):
            atom.charge = pq
    else:
        warnings.warn(
            f"Missing partial charges property ({pq_prop_openff}) from the RDKit Mol"
        )

    at_prop = "BCCAtomTypes"
    if mol.HasProp("%s" % at_prop):
        ats = ast.literal_eval(mol.GetProp("%s" % at_prop))
        assert len(ats) == mol.GetNumAtoms()
        for atom, pq in zip(parmed_structure.atoms, ats):
            atom.type = pq
    else:
        warnings.warn(f"Missing atom types property ({at_prop}) in the RDKit molecule")

    # remove extra bonds in parmed (conversion issues)
    for rd_bond, pmd_bond in list(
        itertools.zip_longest(mol.GetBonds(), parmed_structure.bonds, fillvalue=None)
    )[::-1]:
        if rd_bond is None:
            parmed_structure.bonds.remove(pmd_bond)
            continue

        # sanity
        assert rd_bond.GetBeginAtomIdx() == pmd_bond.atom1.idx
        assert rd_bond.GetEndAtomIdx() == pmd_bond.atom2.idx

    assert mol.GetNumBonds() == len(parmed_structure.bonds)

    return parmed_structure


def get_atoms_bonds_and_parmed_structure(filename, use_general_type=True):
    """

    :param filename:
    :param use_general_type:
    :return: (atoms, bonds, parmed_structure)
    """

    # load if .sdf
    if filename.suffix.lower() == ".sdf":
        warnings.warn(
            f"Reading .sdf ({filename}) via RDKit - using only the first conformer. "
        )
        mol = rdkit.Chem.SDMolSupplier(filename, removeHs=False)[0]
        parmed_structure = pmd_structure_from_rdmol(mol)
    else:
        parmed_structure = parmed.load_file(str(filename), structure=True)

    atoms, bonds = get_atoms_bonds_from_pmd_structure(parmed_structure)

    return atoms, bonds, parmed_structure


def get_atoms_bonds_from_pmd_structure(parmed_structure: parmed.Structure):
    """
    # convert the Parmed atoms into Atom objects.
    """
    atoms = []
    for parmed_atom in parmed_structure.atoms:
        atom_type = parmed_atom.type

        if not atom_type:
            # extract the element symbol from the atomic number
            atom_type = periodic_table.GetElementSymbol(parmed_atom.atomic_number)

        # atom type might be empty if
        if not atom_type:
            # use the atom name as the atom type, e.g. C7
            atom_type = parmed_atom.name

        try:
            atom = Atom(
                name=parmed_atom.name,
                atom_type=atom_type,
                charge=parmed_atom.charge,
            )
        except AttributeError:
            # most likely the charges were missing, manually set the charges to 0
            atom = Atom(
                name=parmed_atom.name,
                atom_type=atom_type,
                charge=0.0,
            )
            logger.warning(
                "One of the input files is missing charges. Setting the charge to 0"
            )
        atom.id = parmed_atom.idx
        atom.position = [parmed_atom.xx, parmed_atom.xy, parmed_atom.xz]
        atom.resname = parmed_atom.residue.name
        atoms.append(atom)

    bonds = [(b.atom1.idx, b.atom2.idx, b.order) for b in parmed_structure.bonds]

    for from_idx, to_idx, bond_type in bonds:
        print("bonds", from_idx, to_idx, bond_type)
        atoms[from_idx].bind_to(atoms[to_idx], bond_type)

    # verify: the ids and the order in the list is the same
    for i, atom in enumerate(atoms):
        assert i == atom.id

    return atoms, bonds


def do_atom_names_have_simple_format(names):
    """
    Check if the atom name is followed by a number, e.g. "C15"
    Note that the full atom name cannot be more than 4 characters.
    This is because the PDB format does not allow for more
    characters which can lead to inconsistencies.

    :param names: a list of atom names
    :type names: list[str]
    :return True if they all follow the correct format.
    """
    for name in names:
        # cannot exceed 4 characters
        if len(name) > 4:
            return False

        # count letters before any digit
        letter_count = 0
        for letter in name:
            if not letter.isalpha():
                break

            letter_count += 1

        # at least one character
        if letter_count == 0:
            return False

        # extrac the number suffix
        atom_number = name[letter_count:]
        try:
            int(atom_number)
        except Exception:
            return False

    return True
