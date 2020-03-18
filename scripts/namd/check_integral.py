import os
import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import sem
import pandas as pd
from itertools import accumulate
from pymbar import timeseries
from collections import OrderedDict

# complex meta
cmeta = {'avdw':
            {'merged_mean':
                 {0.0: -9.220627690769744, 0.05: -9.25651148255814, 0.1: -9.748334155281574, 0.2: -10.72909420193269,
                  0.3: -11.12781176274575, 0.4: -11.165582972342552, 0.5: -10.663698433855382, 0.6: -9.185767299270072,
                  0.7: -4.473792335888038, 0.8: -1.3722448183938687, 0.9: 2.130893802065978, 0.95: 5.987857547484172,
                  1.0: 9.7847188937021}},
        'dvdw':
            {'merged_mean': {1.0: 3.7946149616794407, 0.95: 4.3926539970930225, 0.9: -1.2456096967677441,
                             0.8: -3.3225216594468514, 0.7: -5.45948863712096, 0.6: -6.185408397200932,
                             0.5: -7.775019493502167, 0.4: -8.116299197080293, 0.3: -8.38231872709097,
                             0.2: -8.721359680106632, 0.1: -9.179896201266244, 0.05: -9.288304731756082,
                             0.0: -9.153832989003666}},
        'aele': {
            'merged_mean': {0.0: -8.600337554148618, 0.0909091: -10.561911929356882, 0.272727: -13.426903649635037,
                            0.454545: -15.860555481506164, 0.636364: -16.906734488503833, 0.818182: -17.4291955014995,
                            0.909091: -17.54583662112629, 1.0: -17.88305388203932}},
        'dele': {
            'merged_mean': {1.0: 7.1268798733755405, 0.909091: 7.051049127906978, 0.818182: 7.729964911696102,
                            0.636364: 8.738021292902365, 0.454545: 9.377922425858047, 0.272727: 10.370538453848717,
                            0.0909091: 11.594650083305565, 0.0: 11.976245839416059}}}

# integrate disap q, which should be positive
dele = np.array(list(cmeta['dele']['merged_mean'].items())).T
dele_xs = dele[0][::-1]
dele_de = dele[1][::-1]
dele_int = np.trapz(dele_de, x=dele_xs)
print(dele_int)