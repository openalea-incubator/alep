# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:21:59 2015

@author: ggarin
"""

from alinea.alep.simulation_tools.septo_decomposed import *
import numpy as np

def wed1():
    run_reps_septo(year=2003, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=0.005, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=3., keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               suffix='wednesday', nreps=1)
def wed2():
    run_reps_septo(year=2011, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=0.005, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=3., keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               suffix='wednesday', nreps=1)
def mon3():
    run_reps_septo(year=2011, variety='Mercia', sowing_date='10-15',
               nplants=15, proba_inf=1, septo_delay_dday=5, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=2e-4, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=2., keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               suffix='monday_night', nreps=10)
def mon4():
    run_reps_septo(year=2011, variety='Rht3', sowing_date='10-15',
               nplants=15, proba_inf=1, septo_delay_dday=5, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=2e-4, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=2., keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               suffix='monday_night', nreps=10)
