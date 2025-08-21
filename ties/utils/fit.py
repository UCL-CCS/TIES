"""
Use FEgrow for MCS docking, with TIES for the MCS.
"""

import warnings
import time

from rdkit import Chem
import fegrow

from ties import Pair, Config


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


def mcs_fit(ref: Chem.Mol, lig: Chem.Mol, conformers=1000) -> Chem.Mol:
    start = time.time()

    rlig = fegrow.RMol(lig)
    rlig._save_template(ref)

    config = Config()
    pair = Pair(ref, lig, config=config)
    sup = pair.superimpose()

    # save the mapping to be used later in FEgrow [Optional]
    cc = sup._removed_because_disjointed_cc  #  noqa: F841

    matched_pairs = sup.matched_pairs
    matched_pairs += [
        a for a, _ in sup._removed_because_unmatched_rings
    ]  # partial rings
    matched_pairs += [a for a, _ in sup._removed_due_to_net_charge]  # charges
    matched_pairs += [
        a for a, _ in sup._removed_pairs_with_charge_difference
    ]  # charges

    mapping = [(a.id, b.id) for a, b in matched_pairs if a.element != "H"]

    top1_noh = [a for a in sup.top1 if a.element != "H"]
    top2_noh = [a for a in sup.top2 if a.element != "H"]
    score = len(mapping) * 2 / (len(top1_noh) + len(top2_noh))

    if score < 0.4:
        warnings.warn(
            f"Low MCS score {score}, ditching {lig.GetProp('_Name')} to {ref.GetProp('_Name')}"
        )
        return None

    print(f"Generated SupTop, {time.time() - start:.2f}s")
    start = time.time()

    # mcs_all = get_mcs(ref, lig)["mcs"]
    # # we have to strip the hydrogens from this mapping
    # mapping_heavy = [
    #     (a, b) for a, b in mcs_all if ref.GetAtomWithIdx(a).GetAtomicNum() != 1
    # ]
    # print("rdmap", mapping_heavy)

    rlig.RemoveAllConformers()

    # generate lots of conformers to start from
    try:
        rlig.generate_conformers(
            num_conf=conformers, minimum_conf_rms=0.4, mapping=mapping
        )
    except ValueError as E:
        print(
            "Could not generate conformers",
            ref.GetProp("_Name"),
            lig.GetProp("_Name"),
            E,
        )
        return None

    print(f"Generated Confs, {time.time() - start:.2f}s")

    # attach the MCS details
    rlig.SetProp("MCS(ref,lig)", str(mapping))

    return rlig
