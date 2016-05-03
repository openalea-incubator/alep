# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:21:59 2015

@author: ggarin
"""

from alinea.alep.simulation_tools.septo_decomposed import *

def call1():
    run_reps_septo(year=2011, variety='Mercia', sowing_date='10-15',
               nplants=15, sporulating_fraction=5e-5, leaf_duration=3.,
               suffix='201600401_test', nreps=10)

def call2():
    run_reps_septo(year=2011, variety='Rht3', sowing_date='10-15',
               nplants=15, sporulating_fraction=5e-5, leaf_duration=3.,
               suffix='201600401_test', nreps=10)

def call3():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, sporulating_fraction=4e-3, leaf_duration=3.,
               single_nff=True,
               suffix='20160502_singleNFF_variabilityTrue', nreps=10)

def call4():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, sporulating_fraction=1.5e-3, leaf_duration=3.,
               suffix='201600401_test', nreps=10)
