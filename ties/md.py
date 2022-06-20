import os

from TIES_MD import TIES, cli
from ties_analysis.ties_analysis import Analysis
from ties_analysis.config import Config


class MD():
    def __init__(self, sim_dir, sim_name='complex', fast=False):

        self.sim_dir = sim_dir
        cfg_file = os.path.join(self.sim_dir, 'TIES.cfg')
        md_config = cli.read_config(cfg_file)
        self.md = TIES.TIES(cwd=self.sim_dir, exp_name=sim_name, **md_config)

        if fast:
            #modify md to cut as many corners as possible i.e. reps, windows, sim length
            self.md.reps_per_exec = 1
            self.md.total_reps = 3
            self.md.global_lambdas = [0.00, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 1.00]
            self.md.sampling_per_window = md.sampling_per_window*0.5
            self.md.update_cfg()

    def run(self):
        self.md.run()

    def setup(self):
        self.md.setup()

    def analysis(self, exp_data, legs, analysis_cfg='./analysis.cfg'):
        ana_cfg = Config(analysis_cfg)
        ana_cfg.simulation_legs = legs
        ana_cfg.exp_data = exp_data
        ana = Analysis(ana_cfg)
        ana.run()



