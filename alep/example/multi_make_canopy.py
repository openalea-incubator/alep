from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import make_canopy
import itertools
try:
    import cPickle as pickle
except:
    import pickle

# years = [1998, 2003, 2004, 2010]
years = [2004]
# nb_plants = [3, 6, 10]
nb_plants = [20, 30, 50]
# nb_sects = [1, 3, 5, 7]
nb_sects = [5]
combinations = list(itertools.product(*[years, nb_plants, nb_sects]))
# combinations = [{'yr':i[0], 'nplants':i[1], 'nsect':i[2]} for i in combinations]

def make_canopies((yr, nplants, nsect)):
# def make_canopies(**kwds):
    # import pdb
    # pdb.set_trace()
    # if len(kwds)>0:
        # if 'yr' in kwds.keys:
            # yr = kwds.yr
        # if 'nplants' in kwds.keys:
            # nplants = kwds.nplants
        # if 'nsect' in kwds.keys:
            # nsect = kwds.nsect
        
        make_canopy(start_date = str(yr)+"-10-15 12:00:00", end_date = str(yr+1)+"-08-01 00:00:00",
                nplants = nplants, nsect = nsect, disc_level = 30, 
                dir = './adel/adel_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect')
        
if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(make_canopies, combinations, nb_cpu-1)