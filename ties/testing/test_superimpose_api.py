import json
from glob import glob
import numpy as np

from ties import Pair


def test_superimposition_api():
    """
    Use the Pair api to superimpose two molecules.
    """

    # load
    p = Pair(
        "data/superimpose_api/docked-rmsd_001.pdb", "data/superimpose_api/obj01.pdb"
    )

    # superimpose
    suptop = p.superimpose(
        use_element_in_superimposition=True,
        align_molecules_using_mcs=False,
        ignore_charges_completely=True,
        starting_pairs_heuristics=False,
    )

    # check if the full rings was superimposed
    assert len(suptop) == 5


def test_rmsd_suptop():
    """
    Check the rmsd of the superimposed area.
    :return:
    """
    p = Pair(
        "data/superimpose_api/docked-rmsd_001.pdb", "data/superimpose_api/obj01.pdb"
    )
    suptop = p.superimpose(
        use_element_in_superimposition=True,
        align_molecules_using_mcs=False,
        ignore_charges_completely=True,
        starting_pairs_heuristics=False,
    )

    # calculate the benzene/oxygen ring of the superimposed area
    # this is done without the superimposition
    rmsd = suptop.rmsd()
    np.testing.assert_allclose(rmsd, 3.3815635)
    aligned_rmsd = suptop.align_ligands_using_mcs()
    np.testing.assert_allclose(aligned_rmsd, 0.05504, rtol=1e05)
