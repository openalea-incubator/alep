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
from time import sleep

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
    years = [1998, 2003, 2004, 2010]
    # nb_plants = [3, 6, 10]
    # nb_sects = [1, 3, 5, 7]

    default_yr = 2004
    default_nplants = 6
    default_nsect = 5

    nb_plants = [default_nplants]
    nb_sects = [default_nsect]
    
    full_params = []
    for i_set, param_set in enumerate(param_values):
        full_params+=[np.insert(param_set, 0, [i_set, default_yr, default_nplants, default_nsect]).tolist()]
        full_params+=[np.insert(param_set, 0, [i_set, yr, default_nplants, default_nsect]).tolist() for yr in years if yr!=default_yr]
        full_params+=[np.insert(param_set, 0, [i_set, default_yr, nplants, default_nsect]).tolist() for nplants in nb_plants if nplants!=default_nplants]
        full_params+=[np.insert(param_set, 0, [i_set, default_yr, default_nplants, nsect]).tolist() for nsect in nb_sects if nsect!=default_nsect]

    full_params = [np.insert(param_set, len(param_set), i_simu).tolist() for i_simu, param_set in enumerate(full_params)]
        
    np.savetxt('SGInputFull.txt', full_params, delimiter=' ')
    np.savetxt('SGInputFull2.txt', full_params, delimiter=' ')

def annual_loop(X):
    try:
        for i_rep in range(5):
            g, recorder = run_disease(start_date = str(int(X[1]))+"-10-15 12:00:00", 
                             end_date = str(int(X[1]+1))+"-08-01 00:00:00", 
                             nplants = int(X[2]),
                             nsect = int(X[3]),
                             disc_level = 30, 
                             dir = './adel/adel_'+str(int(X[1]))+'_'+str(int(X[2]))+'pl_'+str(int(X[3]))+'sect',
                             sporulating_fraction=X[4],
                             degree_days_to_chlorosis=X[5],
                             Smin = X[6],
                             Smax = X[7],
                             growth_rate = X[8],
                             age_physio_switch_senescence=X[9],
                             density_dus_emitted_ref = X[10],
                             reduction_by_rain=X[11])
            # recorder = 1
            # sleep(1)
            stored_rec = '.\mercia\\SA\\recorder_'+str(int(X[1]))+'_'+str(int(X[2]))+'pl_'+str(int(X[3]))+'sect'+'_rep'+str(i_rep)+'_sample'+str(int(X[0]))+'.pckl'
            f_rec = open(stored_rec, 'w')
            pickle.dump(recorder, f_rec)
            f_rec.close()
            del recorder
        print 'evaluation successful'
    except:
        print 'evaluation failed'

    i_simu = X[-1]
    try:
        new_values = filter(lambda x: x[-1]!=i_simu, np.loadtxt('SGInputFull2.txt', delimiter=' ').tolist())
        np.savetxt('SGInputFull2.txt', new_values, delimiter=' ')
    except:
        # Last line of the file
        if np.loadtxt('SGInputFull2.txt', delimiter=' ').tolist()[-1] == i_simu:
            np.savetxt('SGInputFull2.txt', [], delimiter=' ')

    
try:
    param_values = np.loadtxt('SGInputFull2.txt', delimiter=' ')
except:
    create_full_params()
    param_values = np.loadtxt('SGInputFull2.txt', delimiter=' ')

if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(annual_loop, param_values.tolist(), max(1, nb_cpu-1))
    # pymap(annual_loop, param_values.tolist(), max(1, nb_cpu))
    # pymap(annual_loop, param_values.tolist(), 1)

