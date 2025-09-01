import glob
from pathlib import Path

from rdkit import Chem


def paths_from_glob(path_pattern):
    """
    Evaluate the pattern and extract the paths
    """
    return [Path(p) for p in glob.glob(path_pattern)]


def write_mol(mol: Chem.Mol, filename):
    with Chem.SDWriter(filename) as SD:
        for conf in mol.GetConformers():
            SD.write(mol, confId=conf.GetId())


def load_conformers(sdf) -> Chem.Mol:
    """
    Load the SDF as one mol with many conformers
    """
    mol = None
    for conf in Chem.SDMolSupplier(sdf, removeHs=False):
        if mol is None:
            mol = conf

        mol.AddConformer(conf.GetConformer(), assignId=True)

    return mol
