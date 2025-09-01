"""
Prepare the structures from smiles:
 - assign charges
 - create a low energy conformer (minimised)

Use OpenFF stack.

Room for improvement:
 - protein binding pocket aware protonation
"""

from pathlib import Path

import openff
from openff.toolkit import ForceField
from openff.units import unit
from openff.interchange import Interchange
import openmm
import openmm.unit
from openff.units import Quantity
from openmmforcefields.generators import GAFFTemplateGenerator
from bs4 import BeautifulSoup

from ties.docking.utils import paths_from_glob

forcefield = ForceField("openff-2.2.1.offxml")


def param_general_conf(
    mol: openff.toolkit.Molecule,
    modify_conformer=True,
    max_min_iterations=10_000,
):
    # Generate multiple conformers and use the lowest energy one
    mol.assign_partial_charges(
        partial_charge_method="am1bccelf10",
    )

    ## add GAFF BCC atom types
    gaff = GAFFTemplateGenerator(molecules=mol)
    bcc_mol = gaff.generate_residue_template(mol)

    soup = BeautifulSoup(bcc_mol, "xml")
    atoms = list(soup.find("Residue").children)
    # establish that the order is the same
    assert all([a.attrs["name"] == off_a.name for a, off_a in zip(atoms, mol.atoms)])
    # get the types
    gaff_types = [a.attrs["type"] for a in atoms if a.name == "Atom"]

    mol.properties["atom.dprop.GAFFAtomType"] = " ".join(gaff_types)

    if not modify_conformer:
        # do not change the existing conformer
        return mol

    ## Generate initial conformers
    mol.clear_conformers()
    mol.generate_conformers(n_conformers=200, rms_cutoff=0.1 * unit.angstrom)

    interchange = Interchange.from_smirnoff(
        force_field=forcefield, topology=mol.to_topology()
    )

    tolerance = Quantity(
        10.0,
        unit.kilojoule_per_mole / unit.nanometer,
    )

    simulation = interchange.to_openmm_simulation(
        integrator=openmm.LangevinMiddleIntegrator(
            293.15 * openmm.unit.kelvin,
            1.0 / openmm.unit.picosecond,
            2.0 * openmm.unit.femtosecond,
        ),
        combine_nonbonded_forces=False,
    )
    simulation.context.computeVirtualSites()

    best_conf = None
    best_energy = None

    for conf in mol.conformers:
        # Set positions and minimize
        simulation.context.setPositions(conf.to_openmm())

        simulation.minimizeEnergy(
            tolerance=tolerance.to_openmm(),
            maxIterations=max_min_iterations,
        )

        state = simulation.context.getState(getPositions=True, getEnergy=True)

        pos = state.getPositions(asNumpy=True)
        pot_ene = state.getPotentialEnergy()

        if best_energy is None or pot_ene < best_energy:
            best_conf = pos
            best_energy = pot_ene

    # keep only the best conformer
    mol.clear_conformers()
    mol.add_conformer(best_conf)
    mol.properties["best_conf_ene"] = best_energy

    return mol


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-sdfs",
        metavar="str",
        dest="sdfs",
        type=paths_from_glob,
        required=False,
        help="An SDF file with molecules",
    )
    parser.add_argument(
        "-smiid",
        metavar="str",
        dest="smiid",
        type=Path,
        required=False,
        default=False,
        help="A files with smiles and ID/name in the second column. ",
    )
    parser.add_argument(
        "-dir",
        metavar="str",
        dest="out_dir",
        type=Path,
        required=False,
        default=Path("mols"),
        help="Directory to which the output should be saved",
    )
    parser.add_argument(
        "-confs",
        metavar="str",
        dest="generate_new_conformers",
        type=Path,
        required=False,
        default=True,
        help="Generate new conformers (ignore existing ones)",
    )
    args = parser.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.sdfs:
        # assume one molecule per SDF
        for sdf in args.sdfs:
            print(sdf)

            mol = openff.toolkit.Molecule.from_file(sdf, allow_undefined_stereo=True)

            param_mol = param_general_conf(mol)

            out_sdf = out_dir / Path(param_mol.name + ".sdf")
            if out_sdf.exists():
                print("file exists already: ", out_sdf)
                continue

            param_mol.to_file(out_sdf, file_format="sdf")

    if args.smiid:
        for line in open(args.smiid).readlines():
            print("Next Smiles", line)

            smi, mol_id = line.strip().rsplit("\t", maxsplit=1)
            mol = openff.toolkit.Molecule.from_smiles(smi, allow_undefined_stereo=True)
            mol.name = mol_id

            param_mol = param_general_conf(mol)

            out_sdf = out_dir / Path(param_mol.name + ".sdf")
            if out_sdf.exists():
                print("file exists already: ", out_sdf)
                continue

            param_mol.to_file(out_sdf, file_format="sdf")
