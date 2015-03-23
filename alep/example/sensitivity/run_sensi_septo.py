from collections import OrderedDict
from multiprocessing import cpu_count
from sensi_septo_morris import *

quantitative_parameters = OrderedDict([('sporulating_fraction', [0.e-4, 2.e-4]),
                                       ('degree_days_to_chlorosis', [130., 230.]),
                                       ('degree_days_to_necrosis', [70., 130.]),
                                       ('degree_days_to_sporulation', [12., 50.]),
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

#qualitative_parameters = OrderedDict([('year', {'default':2004., 'values':[1998., 2003., 2004.]}),
#                                      ('variety', {'default':1, 'values':[1, 2, 3, 4]})])

qualitative_parameters = OrderedDict([('year', {'default':2012, 'values':[2012]})])

variety_code = {1:'Mercia', 2:'Rht3', 3:'Tremie12', 4:'Tremie13'}

list_param_names = qualitative_parameters.keys() + quantitative_parameters.keys()

generate_parameter_set(quantitative_parameters,
                       qualitative_parameters,
                       parameter_range_file = 'param_range_SA.txt',
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