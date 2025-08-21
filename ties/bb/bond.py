class Bond:
    def __init__(self, atom, type):
        self.atom = atom
        self.type = type

    def __repr__(self):
        return f"Bond to {self.atom}"


class Bonds(set):
    def without(self, atom):
        return {bond for bond in self if bond.atom is not atom}
