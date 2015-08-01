from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from alinea.alep.simulation_tools.simulation_tools import make_canopy
from itertools import product
import os

scenarios = [{'year':2011, 'variety':'Mercia', 'sowing_date':'10-15', 
            'nplants':15, 'nsect':7, 'fixed_rep':rep} for rep in range(30)]
scenarios += [{'year':2011, 'variety':'Rht3', 'sowing_date':'10-15', 
            'nplants':15, 'nsect':7, 'fixed_rep':rep} for rep in range(30)]
scenarios += [{'year':2012, 'variety':'Tremie12', 'sowing_date':'10-21', 
            'nplants':15, 'nsect':7, 'fixed_rep':rep} for rep in range(30)]
scenarios += [{'year':2013, 'variety':'Tremie13', 'sowing_date':'10-29', 
            'nplants':15, 'nsect':7, 'fixed_rep':rep} for rep in range(30)]

if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(make_canopy, scenarios, nb_cpu-1)
