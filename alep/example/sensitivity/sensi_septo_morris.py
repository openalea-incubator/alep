
# coding: utf-8

""" Procedure of sensitivity analysis for the septoria model: Method of Morris"""
import sys
sys.path.append("..")
import numpy as np
import random as rd
import pandas as pd
import matplotlib.pyplot as plt
import itertools
from collections import OrderedDict
from SALib.sample.morris import sample
from SALib.util import scale_samples, read_param_file
from SALib.analyze import morris
from openalea.multiprocessing.parallel import pymap
from septo_decomposed import run_disease, make_canopy
from alinea.alep.disease_outputs import get_recorder, mean_by_leaf, mean_audpc_by_leaf
try:
    import cPickle as pickle
except:
    import pickle

from alinea.alep.disease_outputs import *
from alinea.echap.disease.septo_data_reader import *
from alinea.echap.disease.septo_data_treatment import *
from alinea.echap.disease.alep_septo_evaluation import *
    
variety_code = {1:'Mercia', 2:'Rht3', 3:'Tremie12', 4:'Tremie13'}

### Generation of parameter set
def add_parameter(parameter_name = 'sporulating_fraction', interval = [1e-6, 1e-2], filename = 'param_range_SA.txt'):
    """ Add a new line with parameter name and interval of variation in parameter range file.
    """
    f = open(filename, "a")
    s = parameter_name + ' %f %f' % tuple(float(i) for i in interval)
    f.writelines(s + '\n')
    f.close()

def generate_parameter_range_file(quantitative_parameters, filename = 'param_range_SA.txt'):
    """ Generate the file with range of variation for all tested parameters.
    """
    for parameter_name, interval in quantitative_parameters.iteritems():
        add_parameter(parameter_name = parameter_name, interval = interval, filename = filename)

def generate_parameter_set(quantitative_parameters,
                           qualitative_parameters,
                           parameter_range_file = 'param_range_SA.txt',
                           sample_file = 'septo_morris_input.txt',
                           num_trajectories = 10,
                           num_levels = 10):
    """ Generate the file with samples for the analysis of Morris. 
    """   
    # Reset parameter range file
    open(parameter_range_file, 'w').close()
    generate_parameter_range_file(quantitative_parameters = quantitative_parameters, 
                                    filename = parameter_range_file)
    
    # Generate samples
    problem = read_param_file(parameter_range_file)
    param_values = sample(problem, N = num_trajectories, num_levels = num_levels, grid_jump = 5)
    
    # For Method of Morris, save the parameter values in a file (they are needed in the analysis)
    np.savetxt(sample_file, param_values, delimiter=' ')
    
    # Repeat samples for each value of qualitative parameters (avoid repetition of combination with default values)
    full_params = []
    defaults = OrderedDict([(k,v['default']) for k,v in qualitative_parameters.iteritems()])

    for i_set, param_set in enumerate(param_values):
        full_params += [np.insert(param_set, 0, defaults.values()).tolist()]
        
        for param, val in qualitative_parameters.iteritems():
            for value in val['values']:
                if value != val['default']:
                    param_combination = [v for k,v in defaults.iteritems() if k!=param]
                    param_combination.insert(defaults.keys().index(param), value)
                    full_params += [np.insert(param_set, 0, param_combination).tolist()]
    
    # Add indices of sample
    full_params = [np.insert(param_set, 0, i_sample).tolist() for i_sample, param_set in enumerate(full_params)]

    # Save full parameter values
    np.savetxt(sample_file[:-4]+'_full'+sample_file[-4:], full_params, delimiter=' ')

### Run and save simulation
def wheat_path((year, variety, nplants, nsect)):
    if variety.lower().startswith('tremie'):
        variety = 'tremie'
    return '../adel/'+variety.lower()+'_'+str(int(year))+'_'+str(nplants)+'pl_'+str(nsect)+'sect'

def make_canopies((yr, day, variety, nplants, nsect, wheat_path)):
    make_canopy(start_date = str(int(yr-1))+"-10-"+str(int(day))+" 12:00:00", end_date = str(int(yr))+"-08-01 00:00:00",
            variety = variety, nplants = nplants, nsect = nsect, disc_level = 5, dir = wheat_path)
    
def reconstruct_wheat(qualitative_parameters, nb_plants = 6, nb_sects = 5):   
    combinations = list(itertools.product(*[qualitative_parameters['year']['values'], variety_code.values(), [nb_plants], [nb_sects]]))
    combinations = map(lambda x: x + (wheat_path(x),), combinations)
    make_canopies(combinations)

def param_values_to_dict(values, list_param_names):
    keys = ['i_sample'] + list_param_names
    return dict(zip(keys, values))

def annual_loop(sample):
    try:
        # Get indice of simulation
        i_sample = sample.pop('i_sample')

        # Get year of simulation
        if 'year' in sample:
            year = int(sample.pop('year'))
            if year == 2011:
                start_date = "2010-10-15"
                end_date = "2011-06-20"
            elif year == 2012:
                start_date = "2011-10-21"
                end_date = "2012-07-18"
                variety = 'Tremie12'
            elif year == 2013:
                start_date = "2012-10-29"
                end_date = "2013-08-01"
                variety = 'Tremie13'
            else:
                start_date = str(year-1)+"-10-15"
                end_date = str(year)+"-08-01"
        else:
            year = 2005
            start_date = str(year-1)+"-10-15"
            end_date = str(year)+"-08-01"
        start_date += " 12:00:00"
        end_date += " 00:00:00"
            
        # Get variety
        if 'variety' in sample:
            variety_code = {1:'Mercia', 2:'Rht3', 3:'Tremie12', 4:'Tremie13'}
            variety = variety_code[sample.pop('variety')]
        #else:
         #   variety = 'Mercia'

        # Get wheat path
        nplants = 10
        #nplants = 1
        nsect = 7
        w_path = wheat_path((year, variety, nplants, nsect))

        # Run and save simulation
        g, recorder = run_disease(start_date = start_date, 
                         end_date = end_date, 
                         variety = variety, nplants = nplants, nsect = nsect,
                         dir = w_path, reset_reconst = True, 
                         degree_days_to_necrosis = 70.,
                         degree_days_to_sporulation = 70., **sample)
        stored_rec = './'+variety.lower()+'/recorder_'+str(int(i_sample))+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
    except:
        print 'evaluation failed'

def annual_loop_smin(sample):
    try:
        # Get indice of simulation
        i_sample = sample.pop('i_sample')

        # Get year of simulation
        if 'year' in sample:
            year = int(sample.pop('year'))
            if year == 2011:
                start_date = "2010-10-15"
                end_date = "2011-06-20"
            elif year == 2012:
                start_date = "2011-10-21"
                end_date = "2012-07-18"
                variety = 'Tremie12'
            elif year == 2013:
                start_date = "2012-10-29"
                end_date = "2013-08-01"
                variety = 'Tremie13'
            else:
                start_date = str(year-1)+"-10-15"
                end_date = str(year)+"-08-01"
        else:
            year = 2005
            start_date = str(year-1)+"-10-15"
            end_date = str(year)+"-08-01"
        start_date += " 12:00:00"
        end_date += " 00:00:00"
            
        # Get variety
        if 'variety' in sample:
            variety_code = {1:'Mercia', 2:'Rht3', 3:'Tremie12', 4:'Tremie13'}
            variety = variety_code[sample.pop('variety')]
        #else:
         #   variety = 'Mercia'

        # Get wheat path
        nplants = 10
        #nplants = 1
        nsect = 7
        w_path = wheat_path((year, variety, nplants, nsect))

        # Run and save simulation
        g, recorder = run_disease(start_date = start_date, 
                         end_date = end_date, 
                         variety = variety, nplants = nplants, nsect = nsect,
                         dir = w_path, reset_reconst = True, 
                         degree_days_to_necrosis = 70.,
                         degree_days_to_sporulation = 70.,
                         Smin = 0.03,
                         **sample)
        stored_rec = './'+variety.lower()+'/smin_cst/recorder_'+str(int(i_sample))+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
    except:
        print 'evaluation failed'
        
def save_sensitivity_outputs(quantitative_parameters, qualitative_parameters, weather,
                             scenario = {'year':2012},
                             variety = 'Tremie12',
                             filename = 'output_morris_tremie_12.csv'):
    list_param_names = ['i_sample'] + qualitative_parameters.keys() + quantitative_parameters.keys()
    df_in = pd.read_csv('septo_morris_input_full.txt', sep=' ', index_col=0, names=list_param_names)
    df_in = df_in[df_in[scenario.keys()[0]]==scenario.values()[0]]
    for k,v in qualitative_parameters.iteritems():
        if k!=scenario.values()[0]:
            df_in = df_in[df_in[k]==v['default']]
    i_samples = df_in.index

    df_out = pd.DataFrame(columns = ['num_leaf_top'] + list(df_in.columns) + ['normalized_audpc_pycnidia_coverage', 
                                                                              'normalized_audpc_necrosis_percentage',
                                                                              'max_pycnidia_coverage', 
                                                                              'max_necrosis_percentage',
                                                                              'speed', 'date_thr'])
    for i_sample in i_samples:
        stored_rec = './'+variety.lower()+'/recorder_'+str(int(i_sample))+'.pckl'
        recorder = get_recorder(stored_rec)
        df_reco = recorder.data
        df_reco = df_reco.rename(columns = {'num_plant':'plant', 'date':'datetime'})
        df_dates = get_date_threshold(df_reco, weather, variable = 'necrosis_percentage', threshold = 0.2).mean()
        df_speed = get_speed(df_reco, weather, variable = 'necrosis_percentage').mean()
        for lf in np.unique(df_reco['num_leaf_top']):
            df_reco_lf = df_reco[(df_reco['num_leaf_top']==lf)]
            output = {}
            output['num_leaf_top'] = lf
            for col in df_in.columns:
                output[col] = df_in.loc[i_sample, col]
            output['normalized_audpc_pycnidia_coverage'] = df_reco_lf.normalized_audpc_pycnidia_coverage[df_reco_lf['normalized_audpc_pycnidia_coverage'].apply(np.isreal)].mean()
            output['normalized_audpc_necrosis_percentage'] = df_reco_lf.normalized_audpc_necrosis_percentage[df_reco_lf['normalized_audpc_necrosis_percentage'].apply(np.isreal)].mean()
            output['max_pycnidia_coverage'] = min(df_reco_lf.groupby('plant').max()['pycnidia_coverage'].mean(), 100.)
            output['max_necrosis_percentage'] = min(df_reco_lf.groupby('plant').max()['necrosis_percentage'].mean(), 100.)
            output['speed'] = df_speed[lf] if lf in df_speed.index else np.nan
            output['date_thr'] = df_dates[lf] if lf in df_dates.index else np.nan
            df_out = df_out.append(output, ignore_index = True)
        del recorder
    df_out.to_csv(filename)

def morris_analysis(parameter_range_file = 'param_range_SA.txt',
                    input_file = 'septo_morris_input.txt', 
                    output_file = 'septo_morris_output_year_2004.txt'):
    # Perform the sensitivity analysis using the model output
    # Specify which column of the output file to analyze (zero-indexed)
    Si = morris.analyze(parameter_range_file, input_file, output_file,
                        column=0, conf_level=0.95, print_to_console=False)
    # Returns a dictionary with keys 'mu', 'mu_star', 'sigma', and 'mu_star_conf'
    # e.g. Si['mu_star'] contains the mu* value for each parameter, in the
    # same order as the parameter file
    return Si

def read_sensitivity_outputs(filename = 'output_morris_tremie_12.csv'):
     return pd.read_csv(filename)
    
def plot_morris_by_leaf(df_out, variable = 'normalized_audpc_pycnidia_coverage'):
    problem = read_param_file('param_range_SA.txt')
    param_values = np.loadtxt('septo_morris_input.txt')
    fig, axs = plt.subplots(6,2, figsize=(15,30))
    for i, ax in enumerate(axs.flat):
        lf = i+1
        df = df_out[df_out['num_leaf_top']==lf]
        Y = df[variable]
        Si = morris.analyze(problem, param_values, Y, grid_jump = 5, num_levels = 10, conf_level=0.95, print_to_console=False)
        ax.plot(Si['mu_star'], Si['sigma'], 'b*')
        tags = iter(Si['names'])
        for x,y in zip(Si['mu_star'], Si['sigma']):
            ax.annotate(tags.next(), (x,y))
        ax.annotate('Leaf %d' % lf, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)

def scatter_plot_by_leaf(df_out, variable = 'max_sev', parameter = 'degree_days_to_chlorosis'):
    fig, axs = plt.subplots(6,2, figsize=(15,30))
    for i, ax in enumerate(axs.flat):
        lf = i+1
        df = df_out[df_out['num_leaf_top']==lf]
        X = df[parameter]
        Y = df[variable]
        ax.plot(X, Y, 'o ')
        ax.annotate('Leaf %d' % lf, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)