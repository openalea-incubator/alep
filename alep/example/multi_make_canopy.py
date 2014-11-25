from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import make_canopy
import itertools
import os
try:
    import cPickle as pickle
except:
    import pickle
        
def make_canopies((yr, variety, nplants, nsect, wheat_path)):
    make_canopy(start_date = str(int(yr-1))+"-10-15 12:00:00", end_date = str(int(yr))+"-08-01 00:00:00",
            variety = variety, nplants = nplants, nsect = nsect, disc_level = 5, dir = wheat_path)
        
def multi_make_canopies(combinations):
    nb_cpu = cpu_count()
    # pymap(make_canopies, combinations, min(len(combinations), nb_cpu))
    pymap(make_canopies, combinations, 1)
