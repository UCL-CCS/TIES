import logging
import os
import sys
import pathlib
import subprocess
import tempfile

import csv

import parmed
import rdkit.Chem


logger = logging.getLogger(__name__)


class Config:
    """
    The configuration with parameters that can be used to define the entire protocol.
    The settings can be overridden later in the actual classes.

    The settings are stored as properties in the object and can be overwritten.
    """

    def __init__(self, **kwargs):
        # set the path to the scripts
        self.code_root = pathlib.Path(os.path.dirname(__file__))

        # scripts/input files,
        # these are specific to the host
        self.script_dir = self.code_root / "scripts"
        self.namd_script_dir = self.script_dir / "namd"
        self.ambertools_script_dir = self.script_dir / "ambertools"
        self.tleap_check_protein = self.ambertools_script_dir / "check_prot.in"
        self.vmd_vis_script = self.script_dir / "vmd" / "vis_morph.vmd"
        self.vmd_vis_script_sh = self.script_dir / "vmd" / "vis_morph.sh"

        self.unique_atom_names = False

        self._workdir = None
        self._antechamber_dr = False
        self._ambertools_home = None

        self._protein = None

        self._ligand_net_charge = None
        self._atom_pair_q_atol = 0.1
        self._net_charge_threshold = 0.1
        self._redistribute_q_over_unmatched = True
        self._allow_disjoint_components = False
        # use only the element in the superimposition rather than the specific atom type
        self._use_element = False
        self._use_element_in_superimposition = True
        self._partial_ring_allowed = False
        self.starting_pairs_heuristics = True
        # weights in choosing the best MCS, the weighted sum of "(1 - MCS fraction) and RMSD".
        self.weights_ratio = [1, 0]

        # coordinates
        self._align_molecules_using_mcs = False
        self.align_add_removed_mcs = False
        self._use_original_coor = False
        self._coordinates_file = None

        self._ligand_files = set()
        self._manually_matched_atom_pairs = None
        self._manually_mismatched_pairs = None
        self._ligands_contain_q = None

        self._ligand_tleap_in = None
        self._complex_tleap_in = None

        self._superimposition_starting_pairs = None
        self._superimposition_starting_heuristic = 0.6

        self._protein_ff = None
        self._ligand_ff = "leaprc.gaff"
        self._ligand_ff_name = "gaff"

        # MD/NAMD production input file
        self._md_engine = "namd"
        # default to modern CPU version
        self.namd_version = "2.14"
        self._lambda_rep_dir_tree = False

        # experimental
        self._use_hybrid_single_dual_top = False
        self._ignore_charges_completely = False

        self.ligands = None

        # if True, do not allow ligands with the same ligand name
        self.uses_cmd = False

        # assign all the initial configuration values
        self.set_configs(**kwargs)

        # logging
        self.logging_breakdown = False
        self.logging_formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        self.logging_level = logging.INFO

    @property
    def workdir(self):
        """
        Working directory for antechamber calls.
        If None, a temporary directory in /tmp/ will be used.

        :return: Work dir
        :rtype: str
        """
        if self._workdir is None:
            # this is API, so avoid making a directory
            # keep the 'ties' in case the output directory structured needs to be generated/harvested
            temp_dir = tempfile.TemporaryDirectory(prefix="ties")
            self._workdir = pathlib.Path(temp_dir.name)

            # store the temporary directory to trigger its removal when the Config stops existing
            self._workdir_tempdir = temp_dir

        return self._workdir

    @workdir.setter
    def workdir(self, cwd):
        if cwd is not None:
            if type(cwd) is str:
                cwd = pathlib.Path(cwd)
            # user provided
            self._workdir = cwd.absolute()
        else:
            # current directory
            self._workdir = pathlib.Path(os.getcwd()) / "ties"
            self._workdir.mkdir(exist_ok=True)
        logger.debug(f"Working Directory: {self._workdir}")

    # --------------- general
    @property
    def protein(self):
        """
        Path to the protein

        :return: Protein filename
        :rtype: str
        """
        return self._protein

    @protein.setter
    def protein(self, path):
        if path is None:
            logger.info("No protein was select. Skipping protein dG.")
            return

        #  can be loaded by parmed?
        logger.info(f"Opening the protein file {path} with ParmEd..")
        parmed.load_file(str(path), structure=True)
        self._protein = pathlib.Path(path)

    @property
    def ligand_files(self):
        """
        A list of ligand filenames.
        :return:
        """
        return self._ligand_files

    @ligand_files.setter
    def ligand_files(self, files):
        """
        str or Path,
        a Ligand object will converted into a Path
        """
        if isinstance(files, str) or isinstance(files, pathlib.Path):
            files = [pathlib.Path(files)]

        if len(files) < 1:
            logger.error(
                "Please supply at least one ligand file with -l (--ligands). E.g. -l file1.pdb file2.pdb"
            )
            sys.exit(1)

        # check if the ligands have the same extension
        if len({lig.suffix for lig in files}) != 1:
            # fixme - does it actually matter? you parse it anyway?
            logger.error(
                "ligands (-l) have different extensions. "
                "Please ensure all ligands have the same extension"
            )
            sys.exit()

        # files must have unique names
        filenames = [lig.stem for lig in files]
        if len(filenames) != len(set(filenames)):
            logger.error(
                "Some ligand (-l) names are the same. "
                "ensure your ligands have unique names. "
            )
            sys.exit()

        # verify the files with ParmEd if possible
        for ligand in files:
            if ligand.suffix.lower() in (".ac", ".prep"):
                logger.warning(
                    f"Cannot verify .ac/.prep {ligand} with ParmEd. Skipping. "
                )
            elif ligand.suffix.lower() in (".sdf"):
                lig = rdkit.Chem.SDMolSupplier(ligand)
            else:
                logger.debug(f"Opening the ligand {ligand} with ParmEd..")
                lig = parmed.load_file(str(ligand), structure=True)
                # there should be one residue
                if len({a.residue.name for a in lig.atoms}) > 1:
                    logger.warning(
                        f"more than one residue name detected in the ligand {ligand}"
                    )

        # TODO - ensure that it is a connected component and there is no extra atoms

        self._ligand_files = self._ligand_files.union(set(files))

    # --------------- ambertools
    @property
    def ambertools_home(self):
        """
        Ambertools HOME path. If not configured, the env variable AMBERHOME as AMBER_PREFIX will be checked.

        :return: ambertools path
        """
        # ambertools path given?
        path = None
        if self._ambertools_home is None:
            # otherwise check env paths
            if os.getenv("AMBERHOME"):
                path = pathlib.Path(os.getenv("AMBERHOME"))
            elif os.getenv("AMBER_PREFIX"):
                path = pathlib.Path(os.getenv("AMBER_PREFIX"))
            else:
                # try to deduce the env from the location of antechamber binary
                logger.warning(
                    "$AMBERHOME not found. Guessing the location by looking up antechamber. "
                )
                proc = subprocess.run(["which", "antechamber"], capture_output=True)
                decode = proc.stdout.decode("utf-8").strip()
                if "not found" not in decode:
                    ant_path = pathlib.Path(decode)
                    if ant_path.is_file():
                        path = ant_path.parent.parent

            if path is None:
                logger.error(
                    "Error: Cannot find ambertools. $AMBERHOME and $AMBER_PREFIX are empty"
                )
                logger.error("Option 1: source your ambertools script amber.sh")
                logger.error(
                    "Option 2: specify manually the path to amberhome with -ambertools option"
                )
                raise Exception("No ambertools")

            assert path.exists()
            self._ambertools_home = path

        return self._ambertools_home

    @property
    def ambertools_antechamber(self):
        """
        Antechamber path based on the .ambertools_home

        :return:
        """
        return self.ambertools_home / "bin" / "antechamber"

    @property
    def ambertools_parmchk2(self):
        """
        Parmchk2 path based on the .ambertools_home
        :return:
        """
        return self.ambertools_home / "bin" / "parmchk2"

    @property
    def ambertools_tleap(self):
        """
        Tleap path based on the .ambertools_home
        :return:
        """
        return self.ambertools_home / "bin" / "tleap"

    @ambertools_home.setter
    def ambertools_home(self, path):
        if path is None:
            self._ambertools_home = None
            return

        path = pathlib.Path(path)
        if not path.exists():
            logger.error(
                "The provided ambertools home path does not point towards the directory."
                f"{path}"
            )

        self._ambertools_home = path

    @property
    def antechamber_dr(self):
        """
        Whether to use -dr setting when calling antechamber.

        :return:
        """
        return self._antechamber_dr

    @antechamber_dr.setter
    def antechamber_dr(self, bool):
        # using ambertools antechamber dr mode?
        if bool is True:
            self._antechamber_dr = "yes"
        else:
            self._antechamber_dr = "no"
        logger.debug(f"Antechamber dr: {self._antechamber_dr}")

    @property
    def ligand_net_charge(self):
        """
        The ligand charge. If not provided, neutral charge is assumed.
        The charge is necessary for calling antechamber (-nc).

        :return:
        """
        if self._ligand_net_charge is None:
            if self._ligands_contain_q:
                net_q = self._get_first_ligand_net_q()
                logger.info(
                    f"Using the first ligand's net q = {net_q} (from partial charges)"
                )
                self._ligand_net_charge = net_q
            else:
                logger.warning("Ligand net charge not provided (-nc). Assuming 0. ")
                self._ligand_net_charge = 0

        return self._ligand_net_charge

    @ligand_net_charge.setter
    def ligand_net_charge(self, net_charge):
        if net_charge is None:
            logger.warning(
                "Ligand net charge was not supplied (-nc). Neutral charge will be assumed. "
            )

        self._ligand_net_charge = net_charge
        logger.debug(f"Ligand net charge (-nc): {net_charge}")

    @property
    def coordinates_file(self):
        """
        A file from which coordinate can be taken.

        :return:
        """
        return self._coordinates_file

    @coordinates_file.setter
    def coordinates_file(self, file):
        if file is None:
            return

        logger.debug(f"ParmEd: verifying the coordinate file {file}")
        parmed.load_file(file, structure=True)
        # fixme - warn if the atom names are not uniq, warn if there is more than one residue, no water, etc
        self._coordinates_file = file

    @property
    def atom_pair_q_atol(self):
        """
        It defines the maximum difference in charge
        between any two superimposed atoms a1 and a2.
        If the two atoms differ in charge more than this value,
        they will be unmatched and added to the alchemical regions.

        :return: default (0.1e)
        :rtype: float
        """
        return self._atom_pair_q_atol

    @atom_pair_q_atol.setter
    def atom_pair_q_atol(self, atol):
        self._atom_pair_q_atol = atol
        logger.debug(
            f"The maximum acceptable difference in charge between two paired atoms: {atol:.2f}"
        )

    @property
    def net_charge_threshold(self):
        """
        Defines how much the superimposed regions can, in total, differ in charge.
        If the total exceeds the thresholds, atom pairs will be unmatched
        until the threshold is met.

        :return: default (0.1e)
        :rtype: float
        """
        return self._net_charge_threshold

    @net_charge_threshold.setter
    def net_charge_threshold(self, threshold):
        self._net_charge_threshold = threshold
        logger.debug(f"Using MCS net charge difference threshold of {threshold}")

    @property
    def ignore_charges_completely(self):
        """
        Ignore the charges during the superimposition. Useful for debugging.
        :return: default (False)
        :rtype: bool
        """
        return self._ignore_charges_completely

    @ignore_charges_completely.setter
    def ignore_charges_completely(self, bool):
        self._ignore_charges_completely = bool
        if bool:
            logger.debug("Ignoring the charges. ")
            self.redistribute_q_over_unmatched = False

    @property
    def allow_disjoint_components(self):
        """
        Defines whether there might be multiple superimposed areas that are
        separated by alchemical region.

        :return: default (False)
        :rtype: bool
        """
        return self._allow_disjoint_components

    @allow_disjoint_components.setter
    def allow_disjoint_components(self, boolean):
        self._allow_disjoint_components = boolean
        logger.debug(f"Allowing disjoint components: {self._allow_disjoint_components}")

    @property
    def use_element_in_superimposition(self):
        """
        Use element rather than the actual atom type for the superimposition
        during the joint-traversal of the two molecules.

        :return: default (False)
        :rtype: bool
        """
        return self._use_element_in_superimposition

    @use_element_in_superimposition.setter
    def use_element_in_superimposition(self, bool):
        self._use_element_in_superimposition = bool
        if bool:
            logger.warning(
                "Ignoring the specific atom types during superimposition. "
                "The results should not be used in production simulations."
            )

    @property
    def partial_ring_allowed(self):
        return self._partial_ring_allowed

    @property
    def align_molecules_using_mcs(self):
        """
        After determining the maximum common substructure (MCS),
        use it to align the coordinates of the second molecule to the first.

        :return: default (False)
        :rtype: bool
        """
        return self._align_molecules_using_mcs

    @align_molecules_using_mcs.setter
    def align_molecules_using_mcs(self, boolean):
        # align the coordinates in ligZ to the ligA using the MCS
        self._align_molecules_using_mcs = boolean
        # fixme - should be using the MCS before charges change
        logger.debug(f"Will align the coordinates using the final MCS: {boolean}")

    @property
    def use_original_coor(self):
        """
        Antechamber when assigning charges can modify the charges slightly.
        If that's the case, use the original charges in order to correct this slight
        divergence in coordinates.

        :return: default (?)
        :rtype: bool
        """
        return self._use_original_coor

    @use_original_coor.setter
    def use_original_coor(self, boolean):
        # use the original coords because antechamber can change them slightly
        self._use_original_coor = boolean

    def _guess_ligands_contain_q(self):
        """
        Checks if the first .mol2 file contains charges.
        :return:
        """
        # if all ligands are .mol2, then charges are provided
        if all(lig.suffix.lower() == ".mol2" for lig in self.ligand_files):
            # if all atoms have q = 0 that means they're a placeholder
            u = parmed.load_file(str(list(self.ligand_files)[0]), structure=True)
            all_q_0 = all(a.charge == 0 for a in u.atoms)
            if all_q_0:
                return False

            return True

        return False

    def _get_first_ligand_net_q(self):
        """
        :return: Returns the net charge from the parmed file with partial charges.
        """

        # if all atoms have q = 0 that means they're a placeholder
        u = parmed.load_file(str(list(self.ligand_files)[0]), structure=True)
        net_q = sum(a.charge == 0 for a in u.atoms)
        return net_q

    @property
    def ligands_contain_q(self):
        """
        If not provided, it tries to deduce whether charges are provided.
        If all charges are set to 0, then it assumes that charges are not provided.

        If set to False explicitly, charges are ignored and computed again.

        :return: default (None)
        :rtype: bool
        """

        # does the file type contain charges?
        ligand_ext = set(lig.suffix.lower() for lig in self.ligand_files).pop()
        if self._ligands_contain_q is True:
            if ligand_ext in {".mol2", ".ac"}:
                self._ligands_contain_q = True
            else:
                logger.error(
                    "If charges are provided with the ligands, "
                    "the filetypes .mol2 or .ac have to be used."
                )
                sys.exit(1)
        elif self._ligands_contain_q is False:
            self._ligands_contain_q = False
        elif self._ligands_contain_q is None:
            # determine whether charges are provided using the file extensions
            if ligand_ext in {".mol2", ".ac", ".prep"}:
                # if all charges are 0, then recompute
                self._ligands_contain_q = self._guess_ligands_contain_q()
                logger.info("Existing atom charges detected (filetype .ac/.mol2)")
            elif ligand_ext == ".pdb":
                self._ligands_contain_q = False
                logger.debug(
                    "Assuming that charges are not provided in the given .pdb ligand files. "
                )
            else:
                logger.error(f"Error: unrecognized file type {ligand_ext}. ")
                sys.exit(1)

        # fixme?
        if self._ligands_contain_q:
            # leave charges the way they are in the files
            # TODO - ensure that when antechamber_charge_type is accessed ,this function is used? implement in another
            self.antechamber_charge_type = []

        logger.debug(f"Ligand files already contain charges: {self._ligands_contain_q}")

        return self._ligands_contain_q

    @ligands_contain_q.setter
    def ligands_contain_q(self, boolean):
        self._ligands_contain_q = boolean

    @property
    def superimposition_starting_pairs(self):
        """
        Set a starting pair for the superimposition to narrow down the MCS search.
        E.g. "C2-C12"

        :rtype: str
        """
        return self._superimposition_starting_pairs

    @superimposition_starting_pairs.setter
    def superimposition_starting_pairs(self, value):
        if value is None:
            self._superimposition_starting_pairs = None
            return

        # split the different starting pairs
        pairs = value.split(";")
        self._superimposition_starting_pairs = [pair.split("-") for pair in pairs]

    @property
    def superimposition_starting_heuristic(self):
        return self._superimposition_starting_heuristic

    @superimposition_starting_heuristic.setter
    def superimposition_starting_heuristic(self, value):
        self._superimposition_starting_heuristic = value

    @property
    def manually_matched_atom_pairs(self):
        """
        Either a list of pairs or a file with a list of pairs of atoms
        that should be superimposed/matched.

        :return:
        """
        return self._manually_matched_atom_pairs

    @manually_matched_atom_pairs.setter
    def manually_matched_atom_pairs(self, file_or_pairs):
        # A user-provided list of pairs that should be matched
        # This functionality works only if there are two ligands
        if file_or_pairs is None:
            self._manually_matched_atom_pairs = []
            return
        if self._ligand_files is None:
            raise ValueError(
                "Wrong use of Config class. Please set the ._ligand_files attribute first. "
            )

        # TODO
        # pair_format: C1-C7 or C1-C7,C2-C8
        # if different format, it is likely a file, check file

        manually_matched = []
        if file_or_pairs is not None and pathlib.Path(file_or_pairs).is_file():
            with open(file_or_pairs) as IN:
                for left_atom, right_atom in csv.reader(IN, delimiter="-"):
                    manually_matched.append((left_atom.strip(), right_atom.strip()))

        if len(manually_matched) > 1:
            raise NotImplementedError(
                "Currently only one atom pair can be matched - others were not tested"
            )

        self._manually_matched_atom_pairs = manually_matched

    @property
    def manually_mismatched_pairs(self):
        """
        A path to a file with a list of a pairs that should be mismatched.
        """
        if self._manually_mismatched_pairs is None:
            return []

        return self._manually_mismatched_pairs

    @manually_mismatched_pairs.setter
    def manually_mismatched_pairs(self, value):
        mismatch = []

        # only allow the file for now
        if value is not None:
            path = pathlib.Path(value)
            if not path.is_file():
                raise Exception(
                    f"Input file containing mismatching pairs cannot be found: {value}"
                )

            with open(path) as IN:
                for left_atom, right_atom in csv.reader(IN, delimiter="-"):
                    mismatch.append((left_atom.strip(), right_atom.strip()))

        self._manually_mismatched_pairs = mismatch
        return self._manually_mismatched_pairs

    @property
    def protein_ff(self):
        """
        The protein forcefield to be used by ambertools for the protein parameterisation.

        :return: default (leaprc.ff19SB)
        :rtype: string
        """
        if self.protein is None:
            return None

        if self._protein_ff is None:
            logger.warning(
                "Protein FF is not configured in the config.protein_ff. "
                "Setting the default leaprc.ff19SB"
            )
            # fixme - update to a later ff
            self._protein_ff = "leaprc.protein.ff19SB"
        return self._protein_ff

    @protein_ff.setter
    def protein_ff(self, ff):
        self._protein_ff = ff
        logger.debug(f"Protein force field name: {self._protein_ff}")

    @property
    def md_engine(self):
        """
        The MD engine, with the supported values NAMD2.13, NAMD2.14, NAMD3 and OpenMM

        :return: NAMD2.13, NAMD2.14, NAMD3 and OpenMM
        :rtype: string
        """
        return self._md_engine

    @md_engine.setter
    def md_engine(self, value):
        supported = ["NAMD2.13", "NAMD2.14", "NAMD3", "OpenMM"]
        if isinstance(value, str) and value.lower() == "namd2.14":
            self._md_engine = "namd"
            self.namd_version = "2.14"
        elif isinstance(value, str) and value.lower() == "namd2.13":
            self._md_engine = "namd"
            self.namd_version = "2.13"
        elif isinstance(value, str) and value.lower() == "namd3":
            self._md_engine = "namd3"
            self.namd_version = "3"
        elif isinstance(value, str) and value.lower() == "openmm":
            self._md_engine = "openmm"
            self.namd_version = ""
        else:
            raise ValueError(
                "Unknown engine {}. Supported engines {}".format(value, supported)
            )

        logger.debug(f"MD Engine: {value}")

    @property
    def lambda_rep_dir_tree(self):
        return self._lambda_rep_dir_tree

    @lambda_rep_dir_tree.setter
    def lambda_rep_dir_tree(self, value):
        self._lambda_rep_dir_tree = value

    @property
    def ligand_ff(self):
        """
        The forcefield for the ligand.
        """
        return self._ligand_ff

    @property
    def ligand_ff_name(self):
        """
        Either GAFF or GAFF2

        :return:
        """
        return self._ligand_ff_name

    @ligand_ff_name.setter
    def ligand_ff_name(self, atom_type):
        # save also the atom type
        if atom_type == "gaff":
            # they both use the same ff
            self._ligand_ff = "leaprc.gaff"
        elif atom_type == "gaff2":
            self._ligand_ff = "leaprc.gaff2"
        else:
            raise ValueError(
                'Argument -lff cannot be anything else but "gaff" or "gaff2". '
            )

        self._ligand_ff_name = atom_type

    @property
    def redistribute_q_over_unmatched(self):
        """
        The superimposed and matched atoms have every slightly different charges.
        Taking an average charge between any two atoms introduces imbalances
        in the net charge of the alchemical regions, due to
        the different charge distribution.

        :return: default(True)
        """
        return self._redistribute_q_over_unmatched

    @redistribute_q_over_unmatched.setter
    def redistribute_q_over_unmatched(self, boolean):
        # Redistribute charge imbalances created due to atom-pair averaging
        self._redistribute_q_over_unmatched = boolean
        logger.debug(
            f"Distribute the introduced charge disparity in the alchemical region: "
            f"{self._redistribute_q_over_unmatched}"
        )

    @property
    def use_hybrid_single_dual_top(self):
        """
        Hybrid single dual topology (experimental). Currently not implemented.

        :return: default(False).
        """
        # hybrid single dual topology
        return self._use_hybrid_single_dual_top

    @use_hybrid_single_dual_top.setter
    def use_hybrid_single_dual_top(self, boolean):
        if boolean:
            self._use_hybrid_single_dual_top = True
            self._complex_tleap_in = "leap_complex_sdtop.in"
            if not self._ignore_charges_completely:
                raise Exception(
                    "Charges have to be ignored completely when using hybrid single-dual topology."
                )
        else:
            self._complex_tleap_in = "leap_complex.in"

        self._use_hybrid_single_dual_top = boolean

    @property
    def ligand_tleap_in(self):
        """
        The name of the tleap input file for ambertools for the ligand.

        :return: Default ('leap_ligand.in')
        :rtype: string
        """
        if self._ligand_tleap_in is not None:
            # return the user provided filename
            return self._ligand_tleap_in
        if self.use_hybrid_single_dual_top:
            # return the default option for the hybrid
            return "leap_ligand_sdtop.in"

        # return the default
        return "leap_ligand.in"

    @property
    def complex_tleap_in(self):
        """
        The tleap input file for the complex.

        :return: Default 'leap_complex.in'
        :type: string
        """
        if self._complex_tleap_in is None:
            # assume that hybrid single-dual topology is not used
            # this will initiate the standard leap_ligand.in
            # self._complex_tleap_in = 'leap_complex_sdtop.in'
            self.use_hybrid_single_dual_top = False

            self._complex_tleap_in = "leap_complex.in"
        return self._complex_tleap_in

    # PAIR constants configuration
    @property
    def prep_dir(self):
        """
        Path to the `prep` directory. Currently in the `workdir`

        :return: Default (workdir/prep)
        """
        return self.workdir / pathlib.Path("prep")

    @property
    def pair_morphfrcmods_dir(self):
        """
        Path to the .frcmod files for the morph.

        :return: Default (workdir/prep/morph_frcmods)
        """
        return self.prep_dir / "morph_frcmods"

    @property
    def pair_morphfrmocs_tests_dir(self):
        """
        Path to the location where a test is carried out with .frcmod

        :return: Default (workdir/prep/morph_frcmods/tests)
        """
        return self.pair_morphfrcmods_dir / "tests"

    @property
    def pair_unique_atom_names_dir(self):
        """
        Location of the morph files with unique filenames.

        :return: Default (workdir/prep/morph_unique_atom_names)
        """
        return self.prep_dir / "morph_unique_atom_names"

    @property
    def lig_unique_atom_names_dir(self):
        """
        Directory location for files with unique atom names.

        :return: Default (workdir/prep/unique_atom_names)
        """
        return self.prep_dir / "unique_atom_names"

    @property
    def lig_frcmod_dir(self):
        """
        Directory location with the .frcmod created for each ligand.

        :return: Default (workdir/prep/ligand_frcmods)
        """
        return self.prep_dir / "ligand_frcmods"

    @property
    def lig_acprep_dir(self):
        """
        Directory location where the .ac charges are converted into the .mol2 format.

        :return: Default (workdir/prep/acprep_to_mol2)
        """
        return self.prep_dir / "acprep_to_mol2"

    @property
    def lig_dir(self):
        """
        Directory location with the .mol2 files.

        :return: Default (workdir/mol2)
        """
        return self.workdir / "mol2"

    @staticmethod
    def get_element_map():
        """


        :return:
        """
        # Get the mapping of atom types to elements
        element_map_filename = (
            pathlib.Path(os.path.dirname(__file__))
            / "data"
            / "element_atom_type_map.txt"
        )
        # remove the comments lines with #
        lines = filter(
            lambda line: not line.strip().startswith("#") and not line.strip() == "",
            open(element_map_filename).readlines(),
        )
        # convert into a dictionary

        element_map = {}
        for line in lines:
            element, atom_types = line.split("=")

            for atom_type in atom_types.split():
                element_map[atom_type.strip()] = element.strip()

        return element_map

    # fixme - this should be determined at the location where it is relevant rather than here in the conf
    # antechamber parameters, by default compute AM1-BCC charges
    antechamber_charge_type = ["-c", "bcc"]

    def get_serializable(self):
        """
        Get a JSON serializable structure of the config.

        pathlib.Path is not JSON serializable, so replace it with str

        todo - consider capturing all information about the system here,
        including each suptop.get_serializable() so that you can record
        specific information such as the charge changes etc.

        :return: Dictionary {key:value} with the settings
        :rtype: Dictionary
        """

        exclude = [
            "_workdir_tempdir",  # exists (during runtime) for cleaning up purposes
        ]

        host_specific = [
            "code_root",
            "script_dir0",
            "namd_script_dir",
            "ambertools_script_dir",
            "tleap_check_protein",
            "vmd_vis_script",
        ]

        ser = {}
        for k, v in self.__dict__.items():
            if k in host_specific or k in exclude:
                continue

            if type(v) is logging.Formatter:
                continue

            if type(v) is pathlib.PosixPath:
                v = str(v)

            # account for the ligands being pathlib objects
            if k == "ligands" and v is not None:
                # a list of ligands, convert to strings
                v = [str(lig) for lig in v]
            if k == "_ligand_files":
                continue

            ser[k] = v

        return ser

    def set_configs(self, **kwargs):
        # set all the configs one by one
        for k, v in kwargs.items():
            # skip None values, these come from the command line
            if v is None:
                continue

            setattr(self, k, v)
