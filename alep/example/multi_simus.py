""" Example of multiprocessing for repetitions of the same simulation of septoria."""
from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import run_disease
from itertools import product
try:
    import cPickle as pickle
except:
    import pickle

num_simus = range(1,5)
yr = 2004
nplants = 30
nsect = 7
frac = 1e-3
dh = 0.1

# for i_rep in num_simus:
    # g, recorder = run_disease(start_date = str(yr)+"-10-15 12:00:00", 
                         # end_date = str(yr+1)+"-08-01 00:00:00",
                         # nplants = nplants,
                         # nsect = nsect,
                         # disc_level = 5, 
                         # dir = './adel/adel_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect',
                         # sporulating_fraction = frac)
    # stored_rec = '.\mercia\\stability\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep)+'_timetest.pckl'
    # f_rec = open(stored_rec, 'w')
    # pickle.dump(recorder, f_rec)
    # f_rec.close()
    # del recorder

# scenarios = list(product([1e-3], [5], [1, 3], num_simus))
# scenarios = list(product([1e-3], [7], [5], range(5)))
# scenarios = list(product([1e-3], [7], [15], [0]))
# scenarios = list(product([1e-3], [7], [25], [0]))
# dates = [10, 15, 20, 25, 30]
# scenarios = list(product(dates, [1998, 2003, 2004, 2010], [1e-3], [7], [6], range(5)))
# scenarios = list(product([15, 20, 25, 30], [2004], [1e-3], [7], [6], [0]))
# years = [1998, 2003, 2004]
# nb_plants = [6]
# nb_sects = [1, 3, 5, 7, 10, 15]
# layer_thicknesses = [1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]
# scenarios = list(product(years, nb_plants, nb_sects, layer_thicknesses, range(5)))

# def annual_loop((yr, nplants, nsect, dh, i_rep)):
# def annual_loop((frac, nsect, nplants, i_rep)):
def annual_loop(i_rep):
    try:
        print '-----------------------------------------------'
        print i_rep
        g, recorder = run_disease(start_date = str(yr)+"-10-15 12:00:00", 
                             end_date = str(yr+1)+"-08-01 00:00:00",
                             nplants = nplants,
                             nsect = nsect,
                             disc_level = 5, 
                             dir = './adel/adel_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect',
                             sporulating_fraction = frac,
                             layer_thickness = dh)
        stored_rec = '.\mercia\\stability\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep)+'_dh'+str(dh)+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
        print 'simulation succeded'
    except:
        print 'simulation failed'
        
if __name__ == '__main__':
    nb_cpu = cpu_count()
    # pymap(annual_loop, scenarios, nb_cpu-1)
    # pymap(annual_loop, scenarios, 1)
    pymap(annual_loop, num_simus, nb_cpu)
    # pymap(annual_loop, num_simus, 1)