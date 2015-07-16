
# coding: utf-8

""" Procedure of sensitivity analysis for the septoria model: Method of Morris"""
import sys
sys.path.append("..")
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from SALib.sample.morris import sample
from SALib.util import read_param_file
from SALib.analyze import morris

### Generation of parameter set
def add_parameter(parameter_name = 'sporulating_fraction', interval = [1e-6, 1e-2], filename = 'param_range_SA.txt'):
    """ Add a new line with parameter name and interval of variation in parameter range file.
    """
    f = open(filename, "a")
    s = parameter_name + ' %f %f' % tuple(float(i) for i in interval)
    f.writelines(s + '\n')
    f.close()

def generate_parameter_range_file(parameters, filename = 'param_range_SA.txt'):
    """ Generate the file with range of variation for all tested parameters.
    """
    for parameter_name, interval in parameters.iteritems():
        add_parameter(parameter_name = parameter_name, interval = interval, filename = filename)

def generate_parameter_set(parameters,
                           scenarios,
                           parameter_range_file = 'param_range_SA.txt',
                           sample_file = 'morris_input.txt',
                           num_trajectories = 10,
                           num_levels = 10):
    """ Generate the file with samples for the analysis of Morris.

    Parameters
    ----------
    parameters: OrderedDict([('name', [min, max])])
        Names and variation range of quantitative parameters
    scenarios: list[tuple(year, variety)]
        Years and varieties of wheat on which SA is repeated
    """   
    # Reset parameter range file
    open(parameter_range_file, 'w').close()
    generate_parameter_range_file(parameters = parameters, 
                                  filename = parameter_range_file)
    
    # Generate samples
    problem = read_param_file(parameter_range_file)
    param_values = sample(problem, N = num_trajectories, 
                          num_levels = num_levels, grid_jump = 5)
    
    # For Method of Morris, save the parameter values in a file (they are needed in the analysis)
    np.savetxt(sample_file, param_values, delimiter=' ')
    
    # Repeat samples for scenario
    full_params = []
    for scen in scenarios:
        for i_set, param_set in enumerate(param_values):              
            full_params += [np.insert(param_set, 0, scen).tolist()]            
    
    # Add indices of sample
    full_params = [np.insert(param_set, 0, i_sample).tolist()
                    for i_sample, param_set in enumerate(full_params)]

    # Save full parameter values
    np.savetxt(sample_file[:-4]+'_full'+sample_file[-4:], 
               full_params, delimiter=' ')

def param_values_to_dict(values, list_param_names):
    keys = ['i_sample'] + list_param_names
    return dict(zip(keys, values))

def read_sensitivity_outputs(filename = 'septo_morris_output_mercia_2004.csv'):
     return pd.read_csv(filename)
    
def plot_morris_by_leaf(df_out, variable = 'normalized_audpc',
                        parameter_range_file = 'param_range_SA.txt',
                        input_file = 'morris_input.txt'):
    problem = read_param_file(parameter_range_file)
    param_values = np.loadtxt(input_file)
    fig, axs = plt.subplots(6,2, figsize=(15,30))
    for i, ax in enumerate(axs.flat):
        lf = i+1
        df = df_out[df_out['num_leaf_top']==lf]
        Y = df[variable]
        Si = morris.analyze(problem, param_values, Y, grid_jump = 5, 
                            num_levels = 10, conf_level=0.95, 
                            print_to_console=False)
        ax.plot(Si['mu_star'], Si['sigma'], 'b*')
        tags = iter(Si['names'])
        for x,y in zip(Si['mu_star'], Si['sigma']):
            ax.annotate(tags.next(), (x,y))
        ax.annotate('Leaf %d' % lf, xy=(0.05, 0.85), 
                    xycoords='axes fraction', fontsize=18)
        
def plot_morris_by_leaf_synthetic(df_out, variable = 'normalized_audpc',
                                  parameter_range_file = 'param_range_SA.txt',
                                  input_file = 'morris_input.txt', 
                                  leaves = [10, 4, 1]):
    problem = read_param_file(parameter_range_file)
    param_values = np.loadtxt(input_file)
    fig, axs = plt.subplots(1, 3, figsize=(15,6))
    leaves = iter(leaves)
    for i, ax in enumerate(axs.flat):
        lf = next(leaves)
        df = df_out[df_out['num_leaf_top']==lf]
        Y = df[variable]
        Si = morris.analyze(problem, param_values, Y, grid_jump = 5, 
                            num_levels = 10, conf_level=0.95, 
                            print_to_console=False)
        ax.plot(Si['mu_star'], Si['sigma'], 'b*')
        tags = iter(Si['names'])
        for x,y in zip(Si['mu_star'], Si['sigma']):
            ax.annotate(tags.next(), (x,y))
        ax.annotate('Leaf %d' % lf, xy=(0.05, 0.85), 
                    xycoords='axes fraction', fontsize=18)

def scatter_plot_by_leaf(df_out, variable = 'normalized_audpc', parameter = 'Smax'):
    fig, axs = plt.subplots(6,2, figsize=(15,30))
    for i, ax in enumerate(axs.flat):
        lf = i+1
        df = df_out[df_out['num_leaf_top']==lf]
        X = df[parameter]
        Y = df[variable]
        ax.plot(X, Y, 'o ')
        ax.annotate('Leaf %d' % lf, xy=(0.05, 0.85), 
                    xycoords='axes fraction', fontsize=18)