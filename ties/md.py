import os

from TIES_MD import TIES, cli
from ties_analysis.ties_analysis import Analysis
from ties_analysis.config import Config


class MD():
    def __init__(self, sim_dir, sim_name='complex', sim_type='run', fast=False):

        cfg_file = os.path.join(sim_dir, 'TIES.cfg')
        md_config = cli.read_config(cfg_file)
        md = TIES.TIES(cwd=sim_dir, exp_name=sim_name, **md_config)

        if fast:
            #modify md to cut as many corners as possible i.e. reps, windows, sim length
            md.reps_per_exec = 1
            md.total_reps = 1

        if sim_type == 'run':
            md.run()
        elif sim_type == 'setup':
            md.setup()
        else:
            raise ValueError('Unknown run type {}. Please select from [run/setup]'.format(type))

    @staticmethod
    def analysis(exp_data, legs, analysis_cfg='./analysis.cfg'):
        ana_cfg = Config(analysis_cfg)
        ana_cfg.simulation_legs = legs
        ana_cfg.exp_data = exp_data
        ana = Analysis(ana_cfg)
        ana.run()



