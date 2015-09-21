# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:21:59 2015

@author: ggarin
"""

from alinea.alep.simulation_tools.septo_decomposed import *
import numpy as np

def awf1():
    leaf_durations = [2., 2.2, 2.4, 2.6, 2.8, 3.]
    for lf_dur in leaf_durations:
        run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=lf_dur,
               suffix='inoc01_lf_dur_'+str(lf_dur), nreps=1)
        run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=lf_dur,
               suffix='inoc001_lf_dur_'+str(lf_dur), nreps=1)
               
def awf2():
    leaf_durations = [2., 2.2, 2.4, 2.6, 2.8, 3.]
    for lf_dur in leaf_durations:
        run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=lf_dur,
               suffix='inoc01_lf_dur_'+str(lf_dur), nreps=1)
        run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=lf_dur,
               suffix='inoc001_lf_dur_'+str(lf_dur), nreps=1)
               
def awf3():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='all', rh_effect=True,
               suffix='inoc01_rh_all', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='all', rh_effect=True,
               suffix='inoc01_rh_all', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='chlorosis', rh_effect=True,
               suffix='inoc01_rh_chlorosis', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='chlorosis', rh_effect=True,
               suffix='inoc01_rh_chlorosis', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='necrosis', rh_effect=True,
               suffix='inoc01_rh_necrosis', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='chlorosis', rh_effect=True,
               suffix='inoc01_rh_necrosis', nreps=1)
               
def awf4():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='all', rh_effect=True,
               suffix='inoc001_rh_all', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='all', rh_effect=True,
               suffix='inoc001_rh_all', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='chlorosis', rh_effect=True,
               suffix='inoc001_rh_chlorosis', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='chlorosis', rh_effect=True,
               suffix='inoc001_rh_chlorosis', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='necrosis', rh_effect=True,
               suffix='inoc001_rh_necrosis', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
               degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
               apply_rh='chlorosis', rh_effect=True,
               suffix='inoc001_rh_necrosis', nreps=1)
               
def awf5():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
           nplants=15, proba_inf=1, age_infection=False, 
           growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=0.1, reduction_by_rain=0.,
           rain_events_to_empty=5, keep_leaves=False, leaf_duration=3.,
           suffix='inoc01_5rain', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
           nplants=15, proba_inf=1, age_infection=False, 
           growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=0.01, reduction_by_rain=0.,
           rain_events_to_empty=5, keep_leaves=False, leaf_duration=3.,
           suffix='inoc001_5rain', nreps=1)
           
def awf6():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
           nplants=15, proba_inf=1, age_infection=False, 
           growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=0.1, reduction_by_rain=0.,
           rain_events_to_empty=5, keep_leaves=False, leaf_duration=3.,
           suffix='inoc01_HDrop3', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
           nplants=15, proba_inf=1, age_infection=False, 
           growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=0.01, reduction_by_rain=0.,
           rain_events_to_empty=5, keep_leaves=False, leaf_duration=3.,
           suffix='inoc001_HDrop3', nreps=1) 

def bbich1():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
           nplants=15, proba_inf=1, age_infection=False, 
           growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=0.1, reduction_by_rain=0.,
           rain_events_to_empty=10, keep_leaves=False, leaf_duration=2.,
           suffix='inoc01_lf_dur_2', nreps=5)
           
def bbich2():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
           nplants=15, proba_inf=1, age_infection=False, 
           growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=0.01, reduction_by_rain=0.,
           rain_events_to_empty=10, keep_leaves=False, leaf_duration=2.,
           suffix='inoc001_lf_dur_2', nreps=5)
           
def bbich3():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
           nplants=15, proba_inf=1, age_infection=False, 
           growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=0.1, reduction_by_rain=0.,
           rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
           suffix='inoc01_lf_dur_3', nreps=5)
           
def bbich4():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
           nplants=15, proba_inf=1, age_infection=False, 
           growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=0.01, reduction_by_rain=0.,
           rain_events_to_empty=10, keep_leaves=False, leaf_duration=3.,
           suffix='inoc001_lf_dur_3', nreps=5)