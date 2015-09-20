""" Run sensitivity analysis of Morris of septoria model.
"""

from collections import OrderedDict
from disease_sensi_morris import generate_parameter_set
from sensi_septo_tools import variety_code

#parameters = OrderedDict([('sporulating_fraction', [0., 1.e-2]),
#                          ('degree_days_to_chlorosis', [100., 250.]),
#                          ('Smin', [0.001, 0.09]),
#                          ('Smax', [0.1, 1.0]),
#                          ('growth_rate', [0.0001, 0.001]),
#                          ('age_physio_switch_senescence', [0., 1.]),
#                          ('sporulating_capacity', [0., 1.]),
#                          ('density_dus_emitted', [1e3, 5e3]),
#                          ('reduction_by_rain', [0., 1.]), 
#                          ('proba_inf', [0., 1.]), 
#                          ('loss_delay', [72., 168.]), 
#                          ('temp_min', [0., 10.])])

parameters = OrderedDict([('tiller_probability', [0.5, 1.]),
                          ('proba_main_nff', [0.4, 1.]),
                          ('phyllochron', [85., 120.]),
                          ('nb_green_leaves', [3.8, 5.4]),
                          ('scale_leafDim', [0.5, 1.5]),
                          ('scale_stemDim', [0.5, 1.5]),
                          ('scale_leafRate', [0.5, 1.5]),
                          ('scale_stemRate', [0.5, 1.5]),
                          ('scale_fallingRate', [0.5, 1.5])])
                          
v = variety_code()

# scenarios = [(2011, v['Mercia']), 
             # (2011, v['Rht3']),
             # (2012, v['Tremie12']),
             # (2013, v['Tremie13'])]
scenarios = [(2003, v['Custom'])]

list_param_names = ['i_sample', 'i_boot', 'year', 'variety'] + parameters.keys()
nboots = 1
#generate_parameter_set(parameters,
#                       scenarios,
#                       parameter_range_file = './septoria/septo_param_range.txt',
#                       sample_file = './septoria/septo_morris_input.txt',
#                       num_trajectories = 10,
#                       num_levels = 5,
#                       grid_jump = 1,
#                       optimal_trajectories = None,
#                       nboots = nboots)
generate_parameter_set(parameters,
                       scenarios,
                       parameter_range_file = './septo_wheat/septo_param_range.txt',
                       sample_file = './septo_wheat/septo_morris_input.txt',
                       num_trajectories = 10,
                       num_levels = 5,
                       grid_jump = 1,
                       optimal_trajectories = None,
                       nboots = nboots)
    
