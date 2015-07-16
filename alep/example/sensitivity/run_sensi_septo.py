""" Run sensitivity analysis of Morris of septoria model.
"""

from collections import OrderedDict
from multiprocessing import cpu_count
from openalea.multiprocessing.parallel import pymap
from disease_sensi_morris import (generate_parameter_set,
                                  param_values_to_dict)
from sensi_septo_tools import annual_loop, variety_code
import numpy as np

parameters = OrderedDict([('sporulating_fraction', [0., 2.e-4]),
                          ('degree_days_to_chlorosis', [130., 300.]),
                          ('Smin', [0.01, 0.09]),
                          ('Smax', [0.1, 1.0]),
                          ('growth_rate', [0.0001, 0.001]),
                          ('age_physio_switch_senescence', [0.01, 1.]),
                          ('sporulating_capacity', [0., 1.]),
                          ('density_dus_emitted', [1e3, 3e3]),
                          ('reduction_by_rain', [0., 1.]), 
                          ('proba_inf', [0., 1.]), 
                          ('loss_delay', [72., 168.]), 
                          ('temp_min', [0., 10.])])

v = variety_code()
#scenarios = [(yr, 'Tremie13') for yr in [1998, 2003, 2004]] + \
#            [(2004, var) for var in v.itervalues()]
            
scenarios = [(2012, v['Tremie12'])]

list_param_names = ['year', 'variety'] + parameters.keys()

generate_parameter_set(parameters,
                       scenarios,
                       parameter_range_file = 'septo_param_range.txt',
                       sample_file = 'septo_morris_input.txt',
                       num_trajectories = 10,
                       num_levels = 10)

nb_cpu = cpu_count()
filename = 'septo_morris_input_full.txt'
param_values = np.loadtxt(filename, delimiter=' ').tolist()
samples = map(lambda x: param_values_to_dict(x, list_param_names), param_values)

# Run disease simulation
if __name__ == '__main__':
    pymap(annual_loop, samples, nb_cpu)