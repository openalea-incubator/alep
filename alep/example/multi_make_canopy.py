from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import make_canopy
import itertools
try:
    import cPickle as pickle
except:
    import pickle

# years = [1998, 2003, 2004, 2010]
years = [2011]
# starts = [10, 15, 20, 25, 30]
starts = [21]
# nb_plants = [3, 6, 10, 20, 25, 30]
nb_plants = [30]
# nb_plants = [15]
# nb_sects = [5]
nb_sects = [7]
combinations = list(itertools.product(*[years, nb_plants, nb_sects, starts]))
# combinations = [(2010,6,1)]

# for yr, nplants, nsect in combinations:
    # make_canopy(start_date = str(yr)+"-10-15 12:00:00", end_date = str(yr+1)+"-08-01 00:00:00",
            # nplants = nplants, nsect = nsect, disc_level = 5, 
            # dir = './adel/adel_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect')
        
def make_canopies((yr, nplants, nsect, start)):
    make_canopy(start_date = str(yr)+"-10-"+str(start)+" 12:00:00", end_date = str(yr+1)+"-08-01 00:00:00",
            nplants = nplants, nsect = nsect, disc_level = 5, 
            dir = './adel/adel_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_start'+str(start))
        
if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(make_canopies, combinations, nb_cpu)