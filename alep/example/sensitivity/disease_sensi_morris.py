
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
                           num_levels = 10,
                           grid_jump = 5,
                           optimal_trajectories = None,
                           nboots = 5):
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
    for boot in range(nboots):
        param_values = sample(problem, N = num_trajectories, 
                              num_levels = num_levels, grid_jump = grid_jump, 
                              optimal_trajectories=optimal_trajectories)
        
        # For Method of Morris, save the parameter values in a file (they are needed in the analysis)
        sample_file = sample_file[:-4]+'_boot'+str(boot)+sample_file[-4:]
        np.savetxt(sample_file, param_values, delimiter=' ')
        
        # Repeat samples for scenario
        full_params = []
        for scen in scenarios:
            for i_set, param_set in enumerate(param_values):              
                full_params += [np.insert(param_set, 0, scen).tolist()]            

        # Add boot number
        full_params = [np.insert(param_set, 0, boot).tolist()
                        for i_sample, param_set in enumerate(full_params)]
                
        # Add indices of sample
        full_params = [np.insert(param_set, 0, i_sample).tolist()
                        for i_sample, param_set in enumerate(full_params)]

        # Save full parameter values
        np.savetxt(sample_file[:-4]+'_full'+sample_file[-4:], 
                   full_params, delimiter=' ')

def param_values_to_dict(values, list_param_names):
    return dict(zip(list_param_names, values))

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
                    
def plot_morris_3_leaves(df_out, leaves = [10, 5, 1],
                        variable = 'normalized_audpc',
                        parameter_range_file = 'param_range_SA.txt',
                        input_file = 'morris_input.txt'):
    problem = read_param_file(parameter_range_file)
    param_values = np.loadtxt(input_file)
    fig, axs = plt.subplots(1,3, figsize=(18,6))
    leaves = iter(leaves)
    for i, ax in enumerate(axs.flat):
        lf = leaves.next()
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
        ax.set_ylabel('Sigma', fontsize=18)
        ax.set_xlabel('Mu_star', fontsize=18)
    plt.tight_layout()
        
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
                    
def plot_sample_levels(df_out, variable = 'normalized_audpc', leaf = 1, cmap='YlOrRd'):
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    def level_by_column(x, unique_values):
        return unique_values.tolist().index(x)
        
    fig, ax = plt.subplots()
    df = df_out[df_out['num_leaf_top']==leaf]
    df = df.drop(['Unnamed: 0', 'year', 'variety', 'num_leaf_top'], 1)
    for col in df.columns:
        if not col in ['normalized_audpc', 'max_severity']:
            df[col] = df[col].apply(lambda x: level_by_column(x, np.unique(df[col])))
    df = df.sort(variable)
    
    cm = plt.get_cmap(cmap)
    cNorm = colors.Normalize()
    cNorm.autoscale(df[variable])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    nb_output_model = 2
    for i, row in df.iterrows():
        color = scalarMap.to_rgba(row[variable])
        ax.plot(range(len(row)-nb_output_model), row.values[:-nb_output_model], color=color)
    ax.set_xticks(range(len(row)-nb_output_model))
    ax.set_xticklabels(row.keys()[:-nb_output_model])
    ax.set_xlim([-1,len(row)-nb_output_model+1])
    
def scatter_plot_inputs(parameter_range_file = 'param_range_SA.txt',
                        input_file = 'morris_input.txt'):
    import matplotlib.gridspec as gridspec
    df = pd.read_csv(input_file, sep=' ', header=None)
    df_pars = pd.read_csv(parameter_range_file, sep=' ', header=None)
    cols = df_pars[0].values
    df.columns = cols
    gs = gridspec.GridSpec(len(cols), len(cols))
    for i, i_col in enumerate(cols):
        for j, j_col in enumerate(cols):
            if j<=i:
                ax = plt.subplot(gs[len(cols)-i-1, len(cols)-j-1])
                ax.plot(df[i_col], df[j_col], 'bo')
                if j==0:
                    ax.annotate(i_col, xy=(1.05, 0.5), xycoords='axes fraction')
                if len(cols)-i-1==0:
                    ax.set_title(j_col, fontsize=10)
    plt.show(False)