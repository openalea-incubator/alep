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
           nplants=15, proba_inf=1, age_infection=False, growth_rate=0.0006,
           Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=5e-3, reduction_by_rain=0.,
           rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
           suffix='low_inoc_high_disp', nreps=5)
           
def bbich2():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
           nplants=15, proba_inf=1, age_infection=False, growth_rate=0.0006,
           Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=5e-3, reduction_by_rain=0.,
           rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=False,
           suffix='low_inoc_high_disp_keep', nreps=5)
           
def bbich3():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
           nplants=15, proba_inf=1, age_infection=False, growth_rate=0.0006,
           Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=5e-3, reduction_by_rain=0.,
           rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
           suffix='low_inoc_high_disp', nreps=5)
           
def bbich4():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
           nplants=15, proba_inf=1, age_infection=False, growth_rate=0.0006,
           Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=5e-3, reduction_by_rain=0.,
           rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=False,
           suffix='low_inoc_high_disp_keep', nreps=5)
           
def bbich5():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
           nplants=15, proba_inf=1, age_infection=False, growth_rate=0.0006,
           Smin=0.02, degree_days_to_chlorosis=140., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=5e-3, reduction_by_rain=0.,
           rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
           suffix='low_inoc_high_disp_370_140', nreps=5)
           
def bbich6():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-21',
           nplants=15, proba_inf=1, age_infection=False, growth_rate=0.0006,
           Smin=0.02, degree_days_to_chlorosis=140., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
           sporulating_fraction=5e-3, reduction_by_rain=0.,
           rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
           suffix='low_inoc_high_disp_370_140', nreps=5)
           
def bbich7():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-21',
           nplants=15, proba_inf=1, age_infection=False, growth_rate=0.0006,
           Smin=0.02, degree_days_to_chlorosis=160., 
           degree_days_to_necrosis=160., degree_days_to_sporulation=70., 
           sporulating_fraction=5e-3, reduction_by_rain=0.,
           rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
           suffix='low_inoc_high_disp_390_160', nreps=5)
                  
def thu1():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21', nplants=15,
                   proba_inf=1, age_infection=False, growth_rate=0.0006, Smin=0.02, 
                   degree_days_to_chlorosis=160., degree_days_to_necrosis=160., 
                   degree_days_to_sporulation=50., sporulating_fraction=0.01, 
                   reduction_by_rain=0., rain_events_to_empty=5, leaf_duration=2.5, 
                   keep_leaves=True, rh_effect=True, apply_rh='chlorosis',
                   rh_max=50, rh_min=49., suffix='inoc01_em10_rh50_Hrop3', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29', nplants=15,
               proba_inf=1, age_infection=False, growth_rate=0.0006, Smin=0.02, 
               degree_days_to_chlorosis=160., degree_days_to_necrosis=160., 
               degree_days_to_sporulation=50., sporulating_fraction=0.01, 
               reduction_by_rain=0., rain_events_to_empty=5, leaf_duration=2.5, 
               keep_leaves=True, rh_effect=True, apply_rh='chlorosis',
               rh_max=50, rh_min=49., suffix='inoc01_em10_rh50_Hrop3', nreps=5)

def thu2():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21', nplants=15,
                   proba_inf=1, age_infection=False, growth_rate=0.0006, Smin=0.02, 
                   degree_days_to_chlorosis=160., degree_days_to_necrosis=160., 
                   degree_days_to_sporulation=50., sporulating_fraction=0.01, 
                   reduction_by_rain=0., rain_events_to_empty=5, leaf_duration=2.5, 
                   keep_leaves=True, rh_effect=True, apply_rh='chlorosis',
                   rh_max=50, rh_min=49., suffix='inoc01_em10_rh50_Hrop28', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29', nplants=15,
               proba_inf=1, age_infection=False, growth_rate=0.0006, Smin=0.02, 
               degree_days_to_chlorosis=160., degree_days_to_necrosis=160., 
               degree_days_to_sporulation=50., sporulating_fraction=0.01, 
               reduction_by_rain=0., rain_events_to_empty=5, leaf_duration=2.5, 
               keep_leaves=True, rh_effect=True, apply_rh='chlorosis',
               rh_max=50, rh_min=49., suffix='inoc01_em10_rh50_Hrop28', nreps=5)
               
def thu3():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21', nplants=15,
                   proba_inf=1, age_infection=False, growth_rate=0.0006, Smin=0.02, 
                   degree_days_to_chlorosis=160., degree_days_to_necrosis=160., 
                   degree_days_to_sporulation=50., sporulating_fraction=0.01, 
                   reduction_by_rain=0., rain_events_to_empty=5, leaf_duration=2.5, 
                   keep_leaves=True, rh_effect=True, apply_rh='chlorosis',
                   rh_max=50, rh_min=49., suffix='inoc01_em10_rh50_Hrop26', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29', nplants=15,
               proba_inf=1, age_infection=False, growth_rate=0.0006, Smin=0.02, 
               degree_days_to_chlorosis=160., degree_days_to_necrosis=160., 
               degree_days_to_sporulation=50., sporulating_fraction=0.01, 
               reduction_by_rain=0., rain_events_to_empty=5, leaf_duration=2.5, 
               keep_leaves=True, rh_effect=True, apply_rh='chlorosis',
               rh_max=50, rh_min=49., suffix='inoc01_em10_rh50_Hrop26', nreps=5)               
               
def thu4():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21', nplants=15,
               proba_inf=1, age_infection=False, growth_rate=0.0006, Smin=0.02, 
               degree_days_to_chlorosis=160., degree_days_to_necrosis=160., 
               degree_days_to_sporulation=50., sporulating_fraction=0.01, 
               reduction_by_rain=0., rain_events_to_empty=5, leaf_duration=2.5, 
               keep_leaves=True, rh_effect=True, apply_rh='all',
               rh_max=50, rh_min=49., suffix='inoc01_em10_rh50all_Hrop3', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29', nplants=15,
               proba_inf=1, age_infection=False, growth_rate=0.0006, Smin=0.02, 
               degree_days_to_chlorosis=160., degree_days_to_necrosis=160., 
               degree_days_to_sporulation=50., sporulating_fraction=0.01, 
               reduction_by_rain=0., rain_events_to_empty=5, leaf_duration=2.5, 
               keep_leaves=True, rh_effect=True, apply_rh='all',
               rh_max=50, rh_min=49., suffix='inoc01_em10_rh50all_Hrop3', nreps=5)
               
def sun1():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.01, reduction_by_rain=0.,
                   rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='sunday2', nreps=1)

def sun2():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.01, reduction_by_rain=0.,
                   rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='sunday2', nreps=1)
                   
def mon1():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.01, reduction_by_rain=0.,
                   rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='monday', nreps=3)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.03, reduction_by_rain=0.,
                   rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='monday', nreps=3)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.01, reduction_by_rain=0.,
                   rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='monday', nreps=1)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.03, reduction_by_rain=0.,
                   rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='monday', nreps=1)

def mon2():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.01, reduction_by_rain=0.,
                   rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='monday_rain10', nreps=3)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.03, reduction_by_rain=0.,
                   rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='monday_rain10', nreps=3)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.01, reduction_by_rain=0.,
                   rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='monday_rain10', nreps=3)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
                   sporulating_fraction=0.03, reduction_by_rain=0.,
                   rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=True,
                   rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
                   suffix='monday_rain10', nreps=3)
                   
def tue1():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
               suffix='tuesday', nreps=3)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
               suffix='tuesday_rain10', nreps=3)
def tue2():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
               suffix='tuesday', nreps=3)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=2.5, keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=45., rh_min=44.9,
               suffix='tuesday_rain10', nreps=3)
               
def sat1():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=5e-3, reduction_by_rain=0.,
               rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               suffix='sat_modif1_compens', nreps=1)
def sat2():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=160.,
               degree_days_to_necrosis=160., degree_days_to_sporulation=50.,
               sporulating_fraction=0.01, reduction_by_rain=0.,
               rain_events_to_empty=5, leaf_duration=2.5, keep_leaves=False,
               rh_effect=False, apply_rh='all', rh_max=35., rh_min=35.,
               suffix='sat_modif1_compens', nreps=1)