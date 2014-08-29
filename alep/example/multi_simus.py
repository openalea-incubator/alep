""" Example of multiprocessing for repetitions of the same simulation of septoria."""
from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import run_disease
try:
    import cPickle as pickle
except:
    import pickle

num_simus = range(966)
yr = 2004
nplants = 6
nsect = 5
frac = 1e-3

def annual_loop(i_rep):
    try:
        print '-----------------------------------------------'
        print i_rep
        g, recorder = run_disease(start_date = str(yr)+"-10-15 12:00:00", 
                             end_date = str(yr+1)+"-08-01 00:00:00",
                             nplants = nplants,
                             nsect = nsect,
                             disc_level = 30, 
                             dir = './adel/adel_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect',
                             sporulating_fraction = frac)
        stored_rec = '.\mercia\\stability\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep)+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
        print 'simulation succeded'
    except:
        print 'simulation failed'
        
if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(annual_loop, num_simus, nb_cpu-1)
    # pymap(annual_loop, num_simus, 1)