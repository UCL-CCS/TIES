"""
Use RDKit to find the MCS score.
"""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit.Chem import rdFMCS
from rdkit import Chem
from itertools import combinations

from ties.docking.utils import paths_from_glob

Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)


def get_mcs(mol1: Chem.Mol, mol2: Chem.Mol):
    res = rdFMCS.FindMCS(
        [mol1, mol2],
        completeRingsOnly=True,
        ringMatchesRingOnly=True,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        # atomCompare=rdFMCS.AtomCompare.CompareAny,
        bondCompare=rdFMCS.BondCompare.CompareAny,
    )
    match = Chem.MolFromSmarts(res.smartsString)

    # get the matching atoms
    m1m = mol1.GetSubstructMatch(match)
    m2m = mol2.GetSubstructMatch(match)

    m1m_noh = [idx for idx in m1m if mol1.GetAtomWithIdx(idx).GetAtomicNum() != 1]
    m2m_noh = [idx for idx in m2m if mol2.GetAtomWithIdx(idx).GetAtomicNum() != 1]

    # remove any hydrogens
    score = (len(m1m_noh) * 2) / (mol1.GetNumHeavyAtoms() + mol2.GetNumHeavyAtoms())

    # return res.smartsString, res.numAtoms
    return {
        "mcs": list(zip(m1m_noh, m2m_noh)),
        "score": score,
        "m1": mol1.GetProp("_Name"),
        "m2": mol2.GetProp("_Name"),
    }


def arg_parser(args):
    return get_mcs(*args)


def compute_mcs_parallel(files):
    # RDKit load molecules
    mols = [Chem.SDMolSupplier(mol, removeHs=False)[0] for mol in files]

    # ensure the filename is also the ID
    assert all(f.stem == m.GetProp("_Name") for f, m in zip(files, mols))

    # compare all pairs
    with ProcessPoolExecutor(max_workers=12) as executor:
        results = list(executor.map(arg_parser, combinations(mols, r=2)))

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-sdfs",
        metavar="glob or file",
        dest="sdfs",
        type=paths_from_glob,
        required=False,
        default="mols/*sdf",
        help="An SDF file with molecules",
    )
    parser.add_argument(
        "-out",
        metavar="filename",
        dest="out_json",
        type=Path,
        required=False,
        default="mcs_data.json",
        help="Output to store all to all MCS data as computed with RDKit. ",
    )

    args = parser.parse_args()

    if args.out_json.exists():
        raise FileExistsError(f"File exists: {args.out_json}")

    # list files
    results = compute_mcs_parallel(args.sdfs)

    with open(args.out_json, "w") as F:
        json.dump(results, F, indent=4)
