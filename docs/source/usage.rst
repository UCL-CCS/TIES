Usage
=====

TIES offers both a Python API as well as a command line interface.
The API can be very minimal, for example:

.. code-block:: Python

    from ties import Pair


    # load the two ligands and use the default configuration
    pair = Pair('l02.mol2', 'l03.mol2')
    # superimpose the ligands passed above
    hybrid = pair.superimpose()

    # save the results
    hybrid.write_metadata('meta_l02_l03.json')
    hybrid.write_pdb('l02_l03_morph.pdb')
    hybrid.write_mol2('l02_l03_morph.mol2')

This minimal example is sufficient to generate the input for the TIES_MD package
for the simulations in either NAMD or OpenMM.

Note that in this example we do not set any explicit settings.
For that we need to employ the :class:`Config` which we
can then pass to the Pair.

 :class:`Config` contains the settings for all classes in the TIES package, and
 therefore can be used to define a **protocol**.

Whereas all settings can be done in :class:`Config`, for clarity
some can be passed separately here to the :class:`Pair`. This way,
it overwrites the settings in the `config` object:

.. code-block:: Python

    from ties import Pair
    from ties import Config
    from ties import Protein


    config = Config()
    # configure the two settings
    config.workdir = 'ties20'
    config.md_engine = 'openmm'
    # set ligand_net_charge as a parameter,
    # which is equivalent to using config.ligand_net_charge
    pair = Pair('l02.mol2', 'l03.mol2', ligand_net_charge=-1, config=config)
    # rename atoms to help with any issues
    pair.make_atom_names_unique()

    hybrid = pair.superimpose()

    # save meta data to files
    hybrid.write_metadata()
    hybrid.write_pdb()
    hybrid.write_mol2()

    # add the protein for the full RBFE protocol
    config.protein = 'protein.pdb'
    config.protein_ff = 'leaprc.protein.ff14SB'
    protein = Protein(config.protein, config)
    hybrid.prepare_inputs(protein=protein)

Below we show the variation in which we are using :class:`Config` to pass the
net charge of the molecule.

.. code-block:: Python

    import os
    from ties import Pair
    from ties import Config

    # explicitly create config (which will be used by all classes underneath)
    config = Config()
    config.ligand_net_charge = -1

    pair = Pair('l02.mol2', 'l03.mol2', config=config)
    pair.make_atom_names_unique()

    # overwrite the previous config settings with relevant parameters
    hybrid = pair.superimpose(use_element_in_superimposition=True, redistribute_q_over_unmatched=True)

    # prep for the output
    os.mkdir('explicit') if not os.path.exists else None

    # save meta data to specific locations
    hybrid.write_metadata('explicit/result.json')
    hybrid.write_pdb('explicit/result.pdb')
    hybrid.write_mol2('explicit/result.mol2')

    hybrid.prepare_inputs()

Note that there is also the :class:`Ligand` that supports additional operations,
and can be passed directly to :class:`Ligand`.

.. code-block:: Python

    from ties import Ligand


    lig = Ligand('l02_same_atom_name.mol2')

    lig.make_atom_names_correct()
    assert lig.atom_names_correct()

    # prepare the .mol2 input
    lig.antechamber_prepare_mol2()

    # the final .mol2 file
    assert lig.current.exists()

    # Atom naming {new_name: old_name}
    print(lig.renaming_map)
    assert sum('O1' == a for a in lig.renaming_map.values()) == 3
