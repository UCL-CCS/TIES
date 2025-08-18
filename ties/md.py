import os

from TIES_MD import TIES
from TIES_MD import ties_analysis
from openmm import unit
import json


class MD:
    """
    Class a wrapper around TIES_MD API that exposes a simplified interface.

    :param sim_dir: str, points to where the simulation is running i.e where the TIES.cfg file is.
    :param sim_name: str, the prefix to the input param and topo file i.e complex for complex.pdb/prmtop.
    :param fast: boolean, if True the setting is TIES.cfg will be overwritten with minimal TIES protocol.
    """

    def __init__(self, sim_dir, sim_name="complex", fast=False):
        cwd = os.getcwd()
        self.sim_dir = os.path.join(cwd, sim_dir)
        self.analysis_dir = os.path.join(self.sim_dir, "..", "..", "..")
        # This is the main TIES MD object we call it options as the user interacts with this object to change options
        self.options = TIES(cwd=self.sim_dir, exp_name=sim_name)

        if fast:
            # modify md to cut as many corners as possible i.e. reps, windows, sim length
            self.options.total_reps = 3
            self.options.global_lambdas = [
                0.00,
                0.05,
                0.1,
                0.3,
                0.5,
                0.7,
                0.9,
                0.95,
                1.00,
            ]
            self.options.sampling_per_window = 2 * unit.nanoseconds()

        self.options.setup()

    def run(self):
        """
        Wrapper for TIES_MD.TIES.run()

        :return: None
        """
        self.options.run()

    def analysis(self, legs, analysis_cfg="./analysis.cfg"):
        """
        Wrapper for TIES_MD.ties_analysis()

        :param legs: list of strings, these are the thermodynamic legs of the simulation i.e. ['lig', 'com'].
        :param analysis_cfg: str, for what the analysis config file is called.

        :return: None
        """
        os.chdir(self.analysis_dir)
        if not os.path.exists("exp.dat"):
            ties_analysis.make_exp(verbose=False)

        # read the experimental data
        with open("exp.dat") as f:
            data = f.read()
        exp_js = json.loads(data)

        ana_cfg = ties_analysis.Config(analysis_cfg)
        ana_cfg.simulation_legs = legs
        ana_cfg.exp_data = exp_js
        ana = ties_analysis.Analysis(ana_cfg)
        ana.run()
