# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:21:59 2015

@author: ggarin
"""

from alinea.alep.simulation_tools.septo_decomposed import *
import numpy as np

def awf1():
    inocs = np.arange(0.02, 0.2, 0.02)
    inocs = np.insert(inocs, 0, 0.005)
    for inoc in inocs:
        run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                       nplants=15, proba_inf=1, age_infection=False,
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                       degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                       sporulating_fraction=inoc, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='inoc_'+str(inoc), nreps=1)
def awf11():
    inocs = np.arange(0.02, 0.2, 0.02)
    inocs = np.insert(inocs, 0, 0.005)
    for inoc in inocs:
        run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                       nplants=15, proba_inf=1, age_infection=False, 
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                       degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                       sporulating_fraction=inoc, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='inoc_'+str(inoc), nreps=1)
                       
def awf2():
    ref_states = [160,160]
    states = []
    for i in np.arange(-40, 50, 20):
        states.append([ref_states[0]+i, ref_states[1]-i])
    for i_state, state in enumerate(states):
        run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                       nplants=15, proba_inf=1, age_infection=False,
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=state[0], 
                       degree_days_to_necrosis=state[1], degree_days_to_sporulation=50., 
                       sporulating_fraction=0.1, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='states_370_'+str(i_state), nreps=1)

def awf22():
    ref_states = [160,160]
    states = []
    for i in np.arange(-40, 50, 20):
        states.append([ref_states[0]+i, ref_states[1]-i])
    for i_state, state in enumerate(states):
        run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                       nplants=15, proba_inf=1, age_infection=False, 
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=state[0], 
                       degree_days_to_necrosis=state[1], degree_days_to_sporulation=50., 
                       sporulating_fraction=0.1, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='states_370_'+str(i_state), nreps=1)
                       
def awf3():
    ref_states = [160,160]
    states = []
    for i in np.arange(-40, 50, 20):
        states.append([ref_states[0]+i, ref_states[1]-i])
    for i_state, state in enumerate(states):
        run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                       nplants=15, proba_inf=1, age_infection=False,
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=state[0], 
                       degree_days_to_necrosis=state[1], degree_days_to_sporulation=30., 
                       sporulating_fraction=0.1, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='states_390_'+str(i_state), nreps=1)
def awf33():
    ref_states = [160,160]
    states = []
    for i in np.arange(-40, 50, 20):
        states.append([ref_states[0]+i, ref_states[1]-i])
    for i_state, state in enumerate(states):
        run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                       nplants=15, proba_inf=1, age_infection=False, 
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=state[0], 
                       degree_days_to_necrosis=state[1], degree_days_to_sporulation=70., 
                       sporulating_fraction=0.1, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='states_390_'+str(i_state), nreps=1)
                       
def awf4():
    ref_states = [160,160]
    states = []
    for i in np.arange(-40, 50, 20):
        states.append([ref_states[0]+i, ref_states[1]-i])
    for i_state, state in enumerate(states):
        run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                       nplants=15, proba_inf=1, age_infection=False,
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=state[0], 
                       degree_days_to_necrosis=state[1], degree_days_to_sporulation=30., 
                       sporulating_fraction=0.1, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='states_350_'+str(i_state), nreps=1)
def awf44():
    ref_states = [160,160]
    states = []
    for i in np.arange(-40, 50, 20):
        states.append([ref_states[0]+i, ref_states[1]-i])
    for i_state, state in enumerate(states):
        run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                       nplants=15, proba_inf=1, age_infection=False, 
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=state[0], 
                       degree_days_to_necrosis=state[1], degree_days_to_sporulation=30., 
                       sporulating_fraction=0.1, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='states_350_'+str(i_state), nreps=1)
                       
def old():
    run_reps_septo(year=2011, variety='Mercia', sowing_date='10-15',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='calib', nreps=5)
    run_reps_septo(year=2011, variety='Rht3', sowing_date='10-15',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='calib', nreps=5)
                   
def bbich1():
    inocs = np.arange(0.02, 0.2, 0.02)
    inocs = np.insert(inocs, 0, 0.005)
    for inoc in inocs:
        run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                       nplants=15, proba_inf=1, age_infection=False,
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
                       degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                       sporulating_fraction=inoc, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='inoc_180_'+str(inoc), nreps=1)
def bbich11():
    inocs = np.arange(0.02, 0.2, 0.02)
    inocs = np.insert(inocs, 0, 0.005)
    for inoc in inocs:
        run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                       nplants=15, proba_inf=1, age_infection=False, 
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
                       degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                       sporulating_fraction=inoc, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='inoc_180_'+str(inoc), nreps=1)
                       
def bbich2():
    inocs = np.arange(0.02, 0.2, 0.02)
    inocs = np.insert(inocs, 0, 0.005)
    for inoc in inocs:
        run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                       nplants=15, proba_inf=1, age_infection=False,
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                       degree_days_to_necrosis=160., degree_days_to_sporulation=70., 
                       sporulating_fraction=inoc, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='inoc_70_'+str(inoc), nreps=1)

def bbich22():
    inocs = np.arange(0.02, 0.2, 0.02)
    inocs = np.insert(inocs, 0, 0.005)
    for inoc in inocs:
        run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                       nplants=15, proba_inf=1, age_infection=False, 
                       growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                       degree_days_to_necrosis=160., degree_days_to_sporulation=70., 
                       sporulating_fraction=inoc, reduction_by_rain=0.,
                       rain_events_to_empty=10, suffix='inoc_70_'+str(inoc), nreps=1)
                       
def bbich3():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.01, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='smin', nreps=5)
                   
def bbich4():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='ref', nreps=5)
def bbich44():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='ref', nreps=5)
                   
                   
def friday1():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
               degree_days_to_necrosis=140., degree_days_to_sporulation=70., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, suffix='ref_390_HDrop', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
               degree_days_to_necrosis=140., degree_days_to_sporulation=70., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=5, suffix='ref_390_HDrop', nreps=5)
def friday2():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
               degree_days_to_necrosis=140., degree_days_to_sporulation=70., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, suffix='ref_390', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
               degree_days_to_necrosis=140., degree_days_to_sporulation=70., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, suffix='ref_390', nreps=5)
def saturday1():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
               degree_days_to_necrosis=140., degree_days_to_sporulation=70., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=10, keep_leaves=True, suffix='ref_390', nreps=1)
def saturday2():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
               degree_days_to_necrosis=140., degree_days_to_sporulation=70., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=5, keep_leaves=True, suffix='ref_390_5rain', nreps=1)
def saturday3():               
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=True, 
               growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=180., 
               degree_days_to_necrosis=140., degree_days_to_sporulation=70., 
               sporulating_fraction=0.1, reduction_by_rain=0.,
               rain_events_to_empty=5, keep_leaves=True, suffix='ref_390_5_rain_ageinf', nreps=1)
