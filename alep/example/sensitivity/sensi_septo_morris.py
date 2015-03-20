
# coding: utf-8

""" Procedure of sensitivity analysis for the septoria model: Method of Morris"""
import sys
sys.path.append("..")
import numpy as np
import random as rd
import pandas as pd
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


### Generation of parameter set
def add_parameter(parameter_name = 'sporulating_fraction', interval = [1e-6, 1e-2], filename = 'param_range_SA.txt'):
    """ Add a new line with parameter name and interval of variation in parameter range file.
    """
    f = open(filename, "a")
    s = parameter_name + ' %f %f' % tuple(float(i) for i in interval)
    f.writelines(s + '\n')
    f.close()

def generate_parameter_range_file(filename = 'param_range_SA.txt'):
    """ Generate the file with range of variation for all tested parameters.
    """
    for parameter_name, interval in quantitative_parameters.iteritems():
        add_parameter(parameter_name = parameter_name, interval = interval, filename = filename)

def generate_parameter_set(parameter_range_file = 'param_range_SA.txt',
                           sample_file = 'septo_morris_input.txt',
                           num_trajectories = 10,
                           num_levels = 10):
    """ Generate the file with samples for the analysis of Morris. 
    """   
    # Reset parameter range file
    open(parameter_range_file, 'w').close()
    generate_parameter_range_file(filename = parameter_range_file)
    
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
    
def reconstruct_wheat(nb_plants = 6, nb_sects = 5):   
    combinations = list(itertools.product(*[qualitative_parameters['year']['values'], variety_code.values(), [nb_plants], [nb_sects]]))
    combinations = map(lambda x: x + (wheat_path(x),), combinations)
    make_canopies(combinations)

def param_values_to_dict(values):
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
        nplants = 30
        #nplants = 1
        nsect = 7
        w_path = wheat_path((year, variety, nplants, nsect))

        # Run and save simulation
        g, recorder = run_disease(start_date = start_date, 
                         end_date = end_date, 
                         variety = variety, nplants = nplants, nsect = nsect,
                         dir = w_path, **sample)
        stored_rec = './'+variety.lower()+'/recorder_'+str(int(i_sample))+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
    except:
        print 'evaluation failed'

def get_results_audpc(input_file = 'septo_morris_input_full.txt',
                      qualitative_parameter = 'year',
                      value = 2004):   
    # Get i_sample numbers corresponding to qualitative_parameter==value in parameters set
    df = pd.read_csv('septo_morris_input_full.txt', sep=' ', index_col=0, names=list_param_names)
    df = df[df[qualitative_parameter]==value]
    for k,v in qualitative_parameters.iteritems():
        if k!=qualitative_parameter:
            df = df[df[k]==v['default']]
    i_samples = df.index
        
    # Get variety
    if qualitative_parameter == 'variety':
        variety = variety_code[value]
    else:
        variety = 'Mercia'
        
    # Get outputs of simulation for i_samples
    audpcs = []
    for i_sample in i_samples:
        stored_rec = './'+variety.lower()+'/recorder_'+str(int(i_sample))+'.pckl'
        recorder = get_recorder(stored_rec)
        # Output variable: mean audpc on all plants and for leaves 1 to 3
        mean_audpc_f1_to_f3 = mean_audpc_by_leaf(recorder, normalized=True)[['F%d' % lf for lf in range(1, 4)]].mean()
        audpcs.append(mean_audpc_f1_to_f3)
        del recorder
        
    # Save output
    output_file = 'septo_morris_output_'+qualitative_parameter+'_'+str(value)+'.txt'
    np.savetxt(output_file, audpcs, delimiter=' ')

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
