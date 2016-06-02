from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from alinea.alep.simulation_tools.simulation_tools import make_canopy
from itertools import product
import os

nplants = 15
nsect = 7

# scenarios = [{'year':2011, 'variety':'Mercia', 'sowing_date':'10-15', 
            # 'nplants':nplants, 'nsect':nsect, 'fixed_rep':rep} for rep in range(30)]
# scenarios += [{'year':2011, 'variety':'Rht3', 'sowing_date':'10-15', 
            # 'nplants':nplants, 'nsect':nsect, 'fixed_rep':rep} for rep in range(30)]
# scenarios += [{'year':2012, 'variety':'Tremie12', 'sowing_date':'10-21', 
            # 'nplants':nplants, 'nsect':nsect, 'fixed_rep':rep} for rep in range(30)]
scenarios = [{'year':2013, 'variety':'Tremie13', 'sowing_date':'10-29', 
            'nplants':nplants, 'nsect':nsect, 'fixed_rep':rep} for rep in range(30)]
scenarios += [{'year':1999, 'variety':'Tremie13', 'sowing_date':'10-29', 
            'nplants':nplants, 'nsect':nsect, 'fixed_rep':rep} for rep in range(30)]

def make_canopy_helper(kwargs):
    return make_canopy(**kwargs)
            
if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(make_canopy_helper, scenarios, nb_cpu-1)
