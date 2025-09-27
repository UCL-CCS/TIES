"""
Use FEgrow for MCS docking, with TIES for the MCS.
"""

import argparse
import traceback
from pathlib import Path
import warnings
import time

from rdkit import Chem
import fegrow
from openff.toolkit import Molecule

from ties import Pair, Config
from ties.modules.utils.utils import paths_from_glob, write_mol


def mcs_fit(
    ref: Chem.Mol,
    lig: Chem.Mol,
    connected_component_mcs=True,
    conformers=1000,
    score_threshold_mcsheavy=0.4,
) -> Chem.Mol:
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

    if connected_component_mcs:
        matched_pairs += cc

    mapping = [(a.id, b.id) for a, b in matched_pairs if a.element != "H"]

    # calculate the score
    top1_noh = [a for a in sup.top1 if a.element != "H"]
    top2_noh = [a for a in sup.top2 if a.element != "H"]
    score = len(mapping) * 2 / (len(top1_noh) + len(top2_noh))

    if score < score_threshold_mcsheavy:
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


parser = argparse.ArgumentParser(
    description="Generate candidate conformers for relative perturbation calculations, "
    "where each conformer is fitted to the reference structure and"
    " to the binding pocket of the protein. ",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-sdfs",
    metavar="str",
    dest="filenames",
    type=paths_from_glob,
    required=False,
    default="mols/*.sdf",
    help="An SDF file with molecules",
)
parser.add_argument(
    "-ref",
    metavar="str",
    dest="reference",
    type=Path,
    required=True,
    help="A SDF reference molecule file for fitting",
)
parser.add_argument(
    "-dir",
    metavar="str",
    dest="out_dir",
    type=Path,
    required=False,
    default=Path("fitted"),
    help="Directory to which the output should be saved",
)
parser.add_argument(
    "-cc",
    metavar="bool",
    dest="connected_component_mcs",
    type=bool,
    required=False,
    default=True,
    help="Whether to include a connected component that was removed. ",
)
parser.add_argument(
    "-minmcs",
    metavar="threshold",
    dest="score_threshold_mcsheavy",
    type=float,
    required=False,
    default=0.4,
    help="MCS score threshold below which the fitting will be ignored.  ",
)

if __name__ == "__main__":
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(exist_ok=True)

    ref = Molecule.from_file(args.reference, allow_undefined_stereo=True)

    for filename in args.filenames:
        print(filename)

        try:
            mol = Molecule.from_file(filename, allow_undefined_stereo=True)
            confs = mcs_fit(
                ref.to_rdkit(),
                mol.to_rdkit(),
                connected_component_mcs=args.connected_component_mcs,
                score_threshold_mcsheavy=args.score_threshold_mcsheavy,
            )
            write_mol(confs, out_dir / f"{mol.name}.sdf")
        except Exception:
            print(f"Failed {filename} with the exception below: ")
            print(traceback.format_exc())
