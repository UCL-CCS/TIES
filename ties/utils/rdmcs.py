"""
Use RDKit to find the MCS score.
"""

from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem import rdFMCS
from rdkit import Chem
from itertools import combinations


Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)


def get_mcs(mol1: Chem.Mol, mol2: Chem.Mol):
    res = rdFMCS.FindMCS(
        [mol1, mol2],
        completeRingsOnly=True,
        ringMatchesRingOnly=True,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        # atomCompare=rdFMCS.AtomCompare.CompareAny,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
    )
    match = Chem.MolFromSmarts(res.smartsString)

    # get the matching atoms
    m1m = mol1.GetSubstructMatch(match)
    m2m = mol2.GetSubstructMatch(match)

    score = len(m1m) * 2 / (mol1.GetNumAtoms() + mol2.GetNumAtoms())

    # return res.smartsString, res.numAtoms
    return {
        "mcs": list(zip(m1m, m2m)),
        "score": score,
        "m1": mol1.GetProp("_Name"),
        "m2": mol2.GetProp("_Name"),
    }


def arg_parser(args):
    return get_mcs(*args)


def compute_mcs_parallel(files):
    # RDKit load molecules
    mols = [Chem.SDMolSupplier(mol)[0] for mol in files]

    # ensure the filename is also the ID
    assert all(f.stem == m.GetProp("_Name") for f, m in zip(files, mols))

    # compare all pairs
    with ProcessPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(arg_parser, combinations(mols, r=2)))

    return results
