import MDAnalysis as mda

"""
Create four different files containing the complex with .pdb B column
having the value for the contraints. 

This file is for a namd
"""


def create_constraint(input_filename, output, constraint):
    u = mda.Universe(input_filename)
    # for each atom, give the B column the right value
    for atom in u.atoms:
        # ignore water
        if atom.resname == "WAT":
            continue

        # set each atom depending on whether it is a H or not
        if atom.name.upper().startswith("H"):
            atom.tempfactor = 0
        else:
            # restrain the heavy atom
            atom.tempfactor = constraint

    u.atoms.write(output)


def create_4_constraint_files():
    """
    Generate 4 constraint files and return the filenames
    """
    """
    coordinates  ../build/complex.pdb
    constraints  on
    consexp  2
    consref  ../build/complex.pdb ;#need all positions
    conskfile  ../constraint/f4.pdb
    conskcol  B
    """
    # create the 4 constraint files
    filenames = []
    for i in range(1, 4 + 1):
        next_constraint_filename = "constraint%d_complex_merged_solvated.pdb" % i
        create_constraint("complex_merged_solvated.pdb", next_constraint_filename, i)
        filenames.append(next_constraint_filename)

    return filenames
