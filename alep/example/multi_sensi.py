import numpy as np
import random as rd
from SALib.sample import morris_oat
from SALib.util import scale_samples, read_param_file
from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import run_disease
try:
    import cPickle as pickle
except:
    import pickle

def create_full_params():
    # Set random seed (does not affect quasi-random Sobol sampling)
    seed = 1
    np.random.seed(seed)
    rd.seed(seed)

    # Read the parameter range file and generate samples
    param_file = 'params_sensi.txt'
    pf = read_param_file(param_file)

    # Generate samples 
    param_values = morris_oat.sample(10, pf['num_vars'], num_levels = 10, grid_jump = 5)

    # Samples are given in range [0, 1] by default. Rescale them to your parameter bounds. (If using normal distributions, use "scale_samples_normal" instead)
    scale_samples(param_values, pf['bounds'])

    # For Method of Morris, save the parameter values in a file (they are needed in the analysis)
    np.savetxt('SGInput.txt', param_values, delimiter=' ')

    # Add years, number of plants and number of sectors as non crossed discrete parameters
    years = [1998, 2001, 2003, 2010]
    nb_plants = [3, 6, 10]
    nb_sects = [1, 3, 5, 7]

    default_yr = 2001
    default_nplants = 6
    default_nsect = 5

    full_params = []
    for param_set in param_values:
        full_params+=[np.insert(param_set, 0, [yr, default_nplants, default_nsect]).tolist() for yr in years]
        full_params+=[np.insert(param_set, 0, [default_yr, nplants, default_nsect]).tolist() for nplants in nb_plants]
        full_params+=[np.insert(param_set, 0, [default_yr, default_nplants, nsect]).tolist() for nsect in nb_sects]

    np.savetxt('SGInputFull.txt', full_params, delimiter=' ')
    np.savetxt('SGInputFull2.txt', full_params, delimiter=' ')
    
def annual_loop(X):
    try:
        g, recorder = run_disease(start_date = str(int(X[0]))+"-10-15 12:00:00", 
                                 end_date = str(int(X[0]+1))+"-08-01 00:00:00", 
                                 nplants = int(X[1]),
                                 nsect = int(X[2]),
                                 disc_level = 30, 
                                 dir = './adel/adel_'+str(int(X[0]))+'_'+str(int(X[1]))+'pl_'+str(int(X[2]))+'sect',
                                 sporulating_fraction=X[3],
                                 degree_days_to_chlorosis=X[4],
                                 Smin = X[5],
                                 Smax = X[6],
                                 growth_rate = X[7],
                                 age_physio_switch_senescence=X[8],
                                 density_dus_emitted_ref = X[9],
                                 reduction_by_rain=X[10])
        stored_rec = '.\mercia\\SA\\recorder_'+str(int(X[0]))+'_'+str(int(X[1]))+'pl_'+str(int(X[2]))+'sect'+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
        np.savetxt('SGInputFull2.txt', full_params[1:], delimiter=' ')
        print 'evaluation succeded'
    except:
        print 'evaluation failed'

if __name__ == '__main__':
    nb_cpu = cpu_count()
    try:
        param_values = np.loadtxt('SGInputFull2.txt', delimiter=' ')
    except:
        create_full_params()
        param_values = np.loadtxt('SGInputFull2.txt', delimiter=' ')
    pymap(annual_loop, param_values.tolist(), max(1, nb_cpu-2))