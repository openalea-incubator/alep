from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import make_canopy
import itertools
import os
try:
    import cPickle as pickle
except:
    import pickle
        
def make_canopies((yr, day, variety, nplants, nsect, wheat_path)):
    make_canopy(start_date = str(int(yr-1))+"-10-"+str(int(day))+" 12:00:00", end_date = str(int(yr))+"-08-01 00:00:00",
            variety = variety, nplants = nplants, nsect = nsect, disc_level = 5, dir = wheat_path)
        
if __name__ == '__main__':
    combinations = [(2012, 21, 'Tremie12', 30, 7, './adel/tremie_2012_30pl_7sect'), (2013, 29, 'Tremie13', 30, 7, './adel/tremie_2013_30pl_7sect')]
    nb_cpu = cpu_count()
    pymap(make_canopies, combinations, min(len(combinations), nb_cpu))
