""" Run simulations with default parameters for:
    - varied years
    - varied number of plants in canopies
    - varied number of phyto-elements by leaf. """

from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import run_disease
try:
    import cPickle as pickle
except:
    import pickle
    
years = [1998, 2003, 2004, 2010]
nb_plants = [3, 6, 10, 20, 30]
nb_sects = [1, 3, 5, 7]
fracs = [1e-2, 1e-3, 1e-4]

default_yr = 2004
default_nplants = 6
default_nsect = 5
default_frac = 1e-3

scen_default = [(default_yr, default_nplants, default_nsect, default_frac)]
scen_yr = [(yr, default_nplants, default_nsect, default_frac) for yr in years if yr!=default_yr]
scen_pl = [(default_yr, npl, default_nsect, default_frac) for npl in nb_plants if npl!=default_nplants]
scen_ns = [(default_yr, default_nplants, ns, default_frac) for ns in nb_sects if ns!=default_nsect]
scen_fr = [(default_yr, default_nplants, default_nsect, fr) for fr in fracs if fr!=default_frac]
scenarios = scen_default + scen_yr + scen_pl + scen_ns + scen_fr
# scenarios = scen_default

def annual_loop((yr, nplants, nsect, frac)):
    try:
        for i_rep in range(5):
            g, recorder = run_disease(start_date = str(yr)+"-10-15 12:00:00", 
                             end_date = str(yr+1)+"-08-01 00:00:00",
                             nplants = nplants,
                             nsect = nsect,
                             disc_level = 30, 
                             dir = './adel/adel_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect',
                             sporulating_fraction = frac)
            stored_rec = '.\mercia\\default\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep)+'.pckl'
            f_rec = open(stored_rec, 'w')
            pickle.dump(recorder, f_rec)
            f_rec.close()
            del recorder
        print 'evaluation successful'
    except:
        print 'evaluation failed'
        
if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(annual_loop, scenarios, nb_cpu-1)