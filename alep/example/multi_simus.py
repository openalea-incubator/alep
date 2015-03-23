""" Example of multiprocessing for repetitions of the same simulation of septoria."""
from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import run_disease
from itertools import product
try:
    import cPickle as pickle
except:
    import pickle

reps = range(5)
nb_plants = [5, 10, 20, 30, 50]
samples = product(nb_plants, reps)
    
def annual_loop((nb_plants, i_rep)):
    g, recorder = run_disease_and_canopy(start_date="2011-10-21 12:00:00", end_date="2012-08-01 00:00:00",
                                            variety = 'Tremie12', nplants=nb_plants, nsect=7,
                                            dir = dir_path, sporulating_fraction = 1.e-4,
                                            temp_min = 0., degree_days_to_chlorosis = 150,
                                            age_physio_switch_senescence = 0.01,
                                            distri_chlorosis = {'mu':150., 'sigma':30.})
    stored_rec = './tremie/recorder_stabtest_tremie12_'+str(nb_plants)+'pl'+str(i_rep)+'.pckl'
    f_rec = open(stored_rec, 'w')
    pickle.dump(recorder, f_rec)
    f_rec.close()
    del recorder
    del g
        
if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(annual_loop, samples, nb_cpu)
