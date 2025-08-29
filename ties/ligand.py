import os
import subprocess
import shutil
import logging
import uuid
from pathlib import Path

import parmed
import rdkit.Chem

from ties.config import Config
from ties.helpers import get_new_atom_names
from ties import parsing


logger = logging.getLogger(__name__)


class Ligand:
    """
    The ligand helper class. Helps to load and manage the different copies of the ligand file.
    Specifically, it tracks the different copies of the original input files as it is transformed (e.g. charge assignment).

    :param ligand: ligand filepath or RDKit molecule
    :type ligand: string
    :param config: Optional configuration from which the relevant ligand settings can be used
    :type config: :class:`Config`
    :param save: write a file with unique atom names for further inspection
    :type save: bool
    """

    def __init__(self, ligand, config=None, save=True, use_general_type=True):
        """Constructor method"""

        self.save = save
        # save workplace root
        self.config = Config() if config is None else config

        if isinstance(ligand, rdkit.Chem.Mol):
            pmd_structure = parsing.pmd_structure_from_rdmol(ligand)
            atoms, bonds = parsing.get_atoms_bonds_from_pmd_structure(pmd_structure)

            # at the moment we rely on paths as well
            # make this molecule available as a file
            short_uuid = str(uuid.uuid4())[:8]
            lig_path = self.config.workdir / f"{short_uuid}.sdf"
            with rdkit.Chem.SDWriter(lig_path) as SD:
                SD.write(ligand)

            ligand = lig_path
        else:
            # fixme - move use_general_type parameter to config for later
            atoms, bonds, pmd_structure = parsing.get_atoms_bonds_and_parmed_structure(
                ligand, use_general_type=use_general_type
            )

        self.pmd_structure = pmd_structure
        self.atoms = atoms
        self.bonds = bonds

        self.config.ligand_files = ligand

        self.original_input = Path(ligand).absolute()

        # internal name without an extension
        self.internal_name = self.original_input.stem

        # last used representative Path file
        self.current = self.original_input

        self._renaming_map = None
        self.ligand_with_uniq_atom_names = None

    def __repr__(self):
        # return self.original_input.stem
        return self.internal_name

    def use_element(self, value):
        for atom in self.atoms:
            atom.use_general_type = value

    # make this into a python property
    def suffix(self):
        return self.current.suffix.lower()

    def antechamber_prepare_mol2(self, **kwargs):
        """
        Converts the ligand into a .mol2 format.

        BCC charges are generated if missing or requested.
        It calls antechamber (the charge type -c is not used if user prefers to use their charges).
        Any DU atoms created in the antechamber call are removed.

        :param atom_type: Atom type bla bla
        :type atom_type:
        :param net_charge:
        :type net_charge: int
        """
        self.config.set_configs(**kwargs)

        if self.config.ligands_contain_q or not self.config.antechamber_charge_type:
            logger.info(
                f"Antechamber: User-provided atom charges will be reused ({self.current.name})"
            )

        mol2_cwd = self.config.lig_dir / self.internal_name

        # prepare the directory
        mol2_cwd.mkdir(parents=True, exist_ok=True)
        mol2_target = mol2_cwd / f"{self.internal_name}.mol2"

        # do not redo if the target file exists
        if not (mol2_target).is_file():
            log_filename = mol2_cwd / "antechamber.log"
            with open(log_filename, "w") as LOG:
                try:
                    cmd = [
                        self.config.ambertools_antechamber,
                        "-i",
                        self.current,
                        "-fi",
                        self.current.suffix[1:],
                        "-o",
                        mol2_target,
                        "-fo",
                        "mol2",
                        "-at",
                        self.config.ligand_ff_name,
                        "-nc",
                        str(self.config.ligand_net_charge),
                        "-dr",
                        str(self.config.antechamber_dr),
                    ] + self.config.antechamber_charge_type
                    subprocess.run(
                        cmd,
                        cwd=mol2_cwd,
                        stdout=LOG,
                        stderr=LOG,
                        check=True,
                        text=True,
                        timeout=60 * 30,  # 30 minutes
                    )
                except subprocess.CalledProcessError as ProcessError:
                    raise Exception(
                        f"Could not convert the ligand into .mol2 file with antechamber. "
                        f"See the log and its directory: {log_filename} . "
                        f"Command used: {' '.join(map(str, cmd))}"
                    ) from ProcessError
            logger.debug(
                f"Converted {self.original_input} into .mol2, Log: {log_filename}"
            )
        else:
            logger.info(f"File {mol2_target} already exists. Skipping. ")

        self.antechamber_mol2 = mol2_target
        self.current = mol2_target

        # remove any DUMMY DU atoms in the .mol2 atoms
        self._removeDU_atoms()

    def _removeDU_atoms(self):
        """
        Ambertools antechamber creates sometimes DU dummy atoms.
        These are not created when BCC charges are computed from scratch.
        They are only created if you reuse existing charges.
        They appear to be a side effect. We remove the dummy atoms therefore.
        """
        mol2 = parmed.load_file(str(self.current), structure=True)
        # check if there are any DU atoms
        has_DU = any(a.type == "DU" for a in mol2.atoms)
        if not has_DU:
            return

        # make a backup copy before (to simplify naming)
        shutil.move(
            self.current,
            self.current.parent / ("lig.beforeRemovingDU" + self.current.suffix),
        )

        # remove DU type atoms and save the file
        for atom in mol2.atoms:
            if atom.name != "DU":
                continue

            atom.residue.delete_atom(atom)
        # save the updated molecule
        mol2.save(str(self.current))
        logger.debug('Removed dummy atoms with type "DU"')

    def generate_frcmod(self, **kwargs):
        """
        params
         - parmchk2
         - atom_type
        """
        self.config.set_configs(**kwargs)

        logger.debug(f"INFO: frcmod for {self} was computed before. Not repeating.")
        if hasattr(self, "frcmod"):
            return

        # fixme - work on the file handles instaed of the constant stitching
        logger.debug(f"Parmchk2: generate the .frcmod for {self.internal_name}.mol2")

        # prepare cwd
        cwd = self.config.lig_frcmod_dir / self.internal_name
        if not cwd.is_dir():
            cwd.mkdir(parents=True, exist_ok=True)

        target_frcmod = f"{self.internal_name}.frcmod"
        log_filename = cwd / "parmchk2.log"
        with open(log_filename, "w") as LOG:
            try:
                subprocess.run(
                    [
                        self.config.ambertools_parmchk2,
                        "-i",
                        self.current,
                        "-o",
                        target_frcmod,
                        "-f",
                        "mol2",
                        "-s",
                        self.config.ligand_ff_name,
                    ],
                    stdout=LOG,
                    stderr=LOG,
                    check=True,
                    text=True,
                    cwd=cwd,
                    timeout=20,  # 20 seconds
                )
            except subprocess.CalledProcessError as E:
                raise Exception(
                    f"GAFF Error: Could not generate FRCMOD for file: {self.current} . "
                    f"See more here: {log_filename}"
                ) from E

        logger.debug(f"Parmchk2: created frcmod: {target_frcmod}")
        self.frcmod = cwd / target_frcmod

    def correct_atom_names(self):
        """
        Ensure that each atom name:
         - is unique
         - has letter followed by digits
         - has max 4 characters
        E.g. C17, NX23

        :param self.save: if the path is provided, the updated file
            will be saved with the unique names and a handle to the new file (ParmEd) will be returned.
        """
        if self.are_atom_names_correct():
            return

        logger.debug(f"Ligand {self.internal_name} will have its atom names renamed. ")

        ligand = parmed.load_file(str(self.current), structure=True)

        logger.debug(
            f"Atom names in the molecule ({self.original_input}/{self.internal_name}) are either not unique "
            f"or do not follow NameDigit format (e.g. C15). Renaming"
        )
        _, renaming_map = get_new_atom_names(ligand.atoms)
        self._renaming_map = renaming_map
        logger.debug(f"Rename map: {renaming_map}")

        # save the output here
        os.makedirs(self.config.lig_unique_atom_names_dir, exist_ok=True)

        ligand_with_uniq_atom_names = self.config.lig_unique_atom_names_dir / (
            self.internal_name + self.current.suffix
        )
        if self.save:
            ligand.save(str(ligand_with_uniq_atom_names))

        self.ligand_with_uniq_atom_names = ligand_with_uniq_atom_names
        self.parmed = ligand
        # this object is now represented by the updated ligand
        self.current = ligand_with_uniq_atom_names

    def are_atom_names_correct(self):
        """
        Checks if atom names:
         - are unique
         - have a correct format "LettersNumbers" e.g. C17
        """
        ligand = parmed.load_file(str(self.current), structure=True)
        atom_names = [a.name for a in ligand.atoms]

        are_uniqe = len(set(atom_names)) == len(atom_names)

        return are_uniqe
