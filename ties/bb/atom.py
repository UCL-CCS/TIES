import hashlib
import re

import numpy as np

from ties.bb.bond import Bond, Bonds
from ties.bb.gaff_atom_types import gaff_element_map


class Atom:
    counter = 1

    def __init__(self, name, atom_type, charge=0, use_general_type=False):
        self._original_name = None

        self._id = None
        self.name = name
        self._original_name = name.upper()
        self.type = atom_type

        self._resname = None
        self.charge = charge
        self._original_charge = charge

        self.resid = None
        self.bonds: Bonds = Bonds()
        self.use_general_type = use_general_type
        self.hash_value = None

        self._unique_counter = Atom.counter
        Atom.counter += 1

    @property
    def original_name(self):
        # This atom name remains the same. It is never modified.
        return self._original_name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name.upper()

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, id):
        self._id = id

    @property
    def resname(self):
        return self._resname

    @resname.setter
    def resname(self, resname):
        self._resname = resname

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, charge):
        self._charge = charge

    @property
    def original_charge(self):
        return self._original_charge

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, atom_type):
        self._type = atom_type.upper()

        # save the general type
        # fixme - ideally it would use the config class that would use the right mapping
        # strip the element from the associated digits/numbers
        no_trailing_digits = re.match("[A-Za-z]*", self.type)[0]
        self.element = gaff_element_map[no_trailing_digits]

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, pos):
        # switch the type to float32 for any coordinate work with MDAnalysis
        self._position = np.array([pos[0], pos[1], pos[2]], dtype="float32")

    def is_hydrogen(self):
        if self.element == "H":
            return True

        return False

    def bound_to(self, atom):
        for bond in self.bonds:
            if bond.atom is atom:
                return True

        return False

    @property
    def united_charge(self):
        """
        United atom charge: summed charges of this atom and the bonded hydrogens.
        """
        return self.charge + sum(
            bond.atom.charge for bond in self.bonds if bond.atom.is_hydrogen()
        )

    def __hash__(self):
        # Compute the hash key once
        if self.hash_value is not None:
            return self.hash_value

        m = hashlib.md5()
        # fixme - ensure that each node is characterised by its chemistry,
        # fixme - id might not be unique, so check before the input data
        m.update(str(self.charge).encode("utf-8"))
        # use the unique counter to distinguish between created atoms
        m.update(str(self._unique_counter).encode("utf-8"))
        # include the number of bonds
        m.update(str(len(self.bonds)).encode("utf-8"))
        self.hash_value = int(m.hexdigest(), 16)
        return self.hash_value

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def bind_to(self, other, bond_type):
        self.bonds.add(Bond(other, bond_type))
        other.bonds.add(Bond(self, bond_type))

    def eq(self, atom, atol=0):
        """
        Check if the atoms are of the same type and have a charge within the given absolute tolerance.
        """
        if self.type == atom.type and np.isclose(self.charge, atom.charge, atol=atol):
            return True

        return False

    def united_eq(self, atom, atol=0):
        """
        Like .eq, but treat the atoms as united atoms.
        Check if the atoms have the same atom type, and
        if if their charges are within the absolute tolerance.
        If the atoms have hydrogens, add up the attached hydrogens and use a united atom representation.
        """
        if self.type != atom.type:
            return False

        if not np.isclose(self.united_charge, atom.united_charge, atol=atol):
            return False

        return True

    def same_element(self, other):
        # check if the atoms are the same elements
        if self.element == other.element:
            return True
        return False

    def same_type(self, atom):
        if self.type == atom.type:
            return True

        return False
