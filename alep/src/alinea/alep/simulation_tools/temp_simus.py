# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:45:27 2015

@author: ggarin
"""

from alinea.alep.simulation_tools.septo_decomposed import *

# TEMP
def temp1():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=1.,
                   growth_rate=0.0006, density_dus_emitted_ref=1.79e3,
                   nreps=5, suffix='ref')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=1e-2,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=1.,
                   growth_rate=0.0006, density_dus_emitted_ref=1.79e3, 
                   nreps=5, suffix='inoc')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.03, proba_inf=1.,
                   growth_rate=0.0006, density_dus_emitted_ref=1.79e3, 
                   nreps=5, suffix='smin')
                   
def temp2():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=1.,
                   growth_rate=0.001, density_dus_emitted_ref=1.79e3, 
                   nreps=5, suffix='rate')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=1,
                   growth_rate=0.0006, density_dus_emitted_ref=3e3, 
                   nreps=5, suffix='emission')

def temp3():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=0.5,
                   growth_rate=0.0006, density_dus_emitted_ref=1.79e3,
                   nreps=5, suffix='proba')
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=90., degree_days_to_necrosis=90., 
                   degree_days_to_sporulation=70, Smin=0.01, 
                   rh_max=85, rh_min=75, nreps=5, rh_effect=True, suffix='rh_effect')
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50, Smin=0.01, 
                   rh_max=85, rh_min=75, nreps=5, rh_effect=True, suffix='ref')
                   
def temp4():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=90., degree_days_to_necrosis=90., 
                   degree_days_to_sporulation=70, Smin=0.01, 
                   rh_max=85, rh_min=75, nreps=5, suffix='rh_effect')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50., Smin=0.01, 
                   rh_max=85, rh_min=75, nreps=5, suffix='states')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50., Smin=0.1, 
                   rh_max=85, rh_min=75, nreps=5, suffix='big_smin')
                   
def temp5():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50, Smin=0.01, 
                   rh_max=85, rh_min=75, rh_effect=True, apply_rh='necrosis',
                   nreps=5, suffix='rh_1')
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50, Smin=0.01, 
                   rh_max=85, rh_min=75, rh_effect=True, apply_rh='chlorosis',
                   nreps=5, suffix='rh_2')
                   
def temp6():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50, Smin=0.01, 
                   rh_max=85, rh_min=75, rh_effect=True, apply_rh='all',
                   nreps=5, suffix='rh_3')
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=60., degree_days_to_necrosis=95., 
                   degree_days_to_sporulation=95., Smin=0.01, 
                   rh_max=85, rh_min=75, rh_effect=True, apply_rh='all',
                   nreps=5, suffix='rh_4')

def temp7():                
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=2e-2,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50., Smin=0.01, 
                   rh_max=85, rh_min=75, rh_effect=True, apply_rh='all',
                   nreps=5, suffix='rh_5')
                   
def temp8():                
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=1e-2,
                   degree_days_to_chlorosis=120., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=30., nreps=5, suffix='rh_6')
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=1e-2,
                   degree_days_to_chlorosis=120., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=30., nreps=5, 
                   competition='simple', suffix='rh_7')

def temp9():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50., nreps=5, 
                   competition='simple', suffix='rh_8')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=100., degree_days_to_necrosis=100., 
                   degree_days_to_sporulation=50., nreps=5, 
                   competition='simple', suffix='rh_9')
                   
def temp10():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., degree_days_to_necrosis=110., 
                   degree_days_to_sporulation=50, 
                   proba_inf=0.5, nreps=5, suffix='elong_ref')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., degree_days_to_necrosis=110., 
                   degree_days_to_sporulation=50, Smin=0.01,
                   proba_inf=0.5, nreps=5, suffix='elong_smin')
                   
def temp11():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=1e-3,
                   degree_days_to_chlorosis=120., degree_days_to_necrosis=110., 
                   degree_days_to_sporulation=50, proba_inf=0.5, nreps=5, 
                   suffix='elong_inoc_1')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=1e-3,
                   degree_days_to_chlorosis=120., degree_days_to_necrosis=110., 
                   degree_days_to_sporulation=50, nreps=5, 
                   suffix='elong_inoc_2')
                   
def temp12():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=150., degree_days_to_necrosis=110., 
                   degree_days_to_sporulation=20, proba_inf=0.5, nreps=5, 
                   suffix='elong_states_1')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=1e-3,
                   degree_days_to_chlorosis=80., degree_days_to_necrosis=110., 
                   degree_days_to_sporulation=80, nreps=5, 
                   suffix='elong_states_2')
                   
def temp13():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,proba_inf=0.3,
                   density_dus_emitted_ref=5e4, nreps=5, 
                   suffix='big_emission_1')
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3,proba_inf=0.3,
                   density_dus_emitted_ref=8e4, nreps=5, 
                   suffix='big_emission_2')

def temp14():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, sporulating_fraction=5e-3, proba_inf=0.3,
                   density_dus_emitted_ref=5e4, nreps=5, age_infection=True,
                   suffix='age_infection')
                   
def temp15():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.001,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_age_logi', nreps=5)

def temp16():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.001,
                   Smin=0.02, degree_days_to_chlorosis=150., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_rate', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.001,
                   Smin=0.02, degree_days_to_chlorosis=150., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_rate', nreps=5)
                   
def temp17():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.03, degree_days_to_chlorosis=150., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_smin', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.03, degree_days_to_chlorosis=150., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_smin', nreps=5)

def temp18():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.02, degree_days_to_chlorosis=170., 
                   degree_days_to_necrosis=90., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_states', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.02, degree_days_to_chlorosis=170., 
                   degree_days_to_necrosis=90., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_states', nreps=5)
                   
def temp19():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_states_2', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_states_2', nreps=1)
                   
def temp20():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.001,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_age_lin', nreps=1)
                   
def temp21():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.001,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='new_calib_age_lin', nreps=1)
                   
def temp22():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                  nplants=15, proba_inf=1, age_infection=True,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=100., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='new_calib_states_3', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=True,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=100., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='new_calib_states_3', nreps=5)
                   
def temp23():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                  nplants=15, proba_inf=1, age_infection=False,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=100., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='new_calib_states_4', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=100., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='new_calib_states_4', nreps=5)
                   
def temp24():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                  nplants=15, proba_inf=1, age_infection=True,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=20., 
                   degree_days_to_necrosis=280., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='new_calib_states_5', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=True,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=20., 
                   degree_days_to_necrosis=280., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='new_calib_states_5', nreps=5)
                   
def temp25():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                  nplants=15, proba_inf=1, age_infection=False,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=20., 
                   degree_days_to_necrosis=280., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='new_calib_states_6', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=20., 
                   degree_days_to_necrosis=280., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='new_calib_states_6', nreps=5)
                   
def temp26():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                  nplants=15, proba_inf=1, age_infection=True,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=20., 
                   degree_days_to_necrosis=280., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='small_1', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=True,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=20., 
                   degree_days_to_necrosis=280., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='small_1', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                  nplants=15, proba_inf=1, age_infection=False,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=20., 
                   degree_days_to_necrosis=280., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='small_2', nreps=5)

def temp27():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False,
                   density_dus_emitted_ref=1e5, growth_rate=0.0006, 
                   Smin=0.05, degree_days_to_chlorosis=20., 
                   degree_days_to_necrosis=280., degree_days_to_sporulation=60.,
                   sporulating_fraction=7e-2, suffix='small_2', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='small_3', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=0.5, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='small_3', nreps=1)

def temp28():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=0.5, age_infection=False, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='small_4', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=0.5, age_infection=False, 
                   density_dus_emitted_ref=1e5, growth_rate =0.0006,
                   Smin=0.02, degree_days_to_chlorosis=130., 
                   degree_days_to_necrosis=110., degree_days_to_sporulation=30., 
                   sporulating_fraction=7e-2, suffix='small_4', nreps=1)
                   
def temp29():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False, 
                   density_dus_emitted_ref=1e5, growth_rate=0.0006,
                   Smin=0.03, degree_days_to_chlorosis=120., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60., 
                   sporulating_fraction=7e-2, suffix='120_200_60_F', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   density_dus_emitted_ref=1e5, growth_rate=0.0006,
                   Smin=0.03, degree_days_to_chlorosis=120., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60., 
                   sporulating_fraction=7e-2, suffix='120_200_60_F', nreps=1)
                   
def temp30():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate=0.0006,
                   Smin=0.03, degree_days_to_chlorosis=120., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60., 
                   sporulating_fraction=7e-2, suffix='120_200_60_T', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate=0.0006,
                   Smin=0.03, degree_days_to_chlorosis=120., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60., 
                   sporulating_fraction=7e-2, suffix='120_200_60_T', nreps=1)
                   
def temp31():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate=0.0006,
                   Smin=0.03, degree_days_to_chlorosis=100., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60., 
                   sporulating_fraction=7e-2, suffix='120_200_60_T2', nreps=1)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=True, 
                   density_dus_emitted_ref=1e5, growth_rate=0.0006,
                   Smin=0.03, degree_days_to_chlorosis=100., 
                   degree_days_to_necrosis=200., degree_days_to_sporulation=60., 
                   sporulating_fraction=7e-2, suffix='120_200_60_T2', nreps=1)
                   
def temp32():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='debug_1', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='debug_1', nreps=5)
                   
def temp33():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=True,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='debug_2', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=True, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='debug_2', nreps=5)
                   
def temp34():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=30, proba_inf=1, age_infection=True, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='debug_2', nreps=1)
                   
def temp35():
    run_reps_septo(year=2012, variety='Custom', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=True, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='debug_2', nreps=5)
                   
def temp36():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=3, suffix='rain_1', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=3, suffix='rain_1', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=5, suffix='rain_2', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=5, suffix='rain_2', nreps=5)
                   
def temp37():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.01, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='inoc_1', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.01, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='inoc_1', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.05, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='inoc_2', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.05, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='inoc_2', nreps=5)
                   
def temp38():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=250., 
                   degree_days_to_necrosis=50., degree_days_to_sporulation=100., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='new_states_1', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=250., 
                   degree_days_to_necrosis=50., degree_days_to_sporulation=100., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='new_states_1', nreps=5)
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=200., 
                   degree_days_to_necrosis=100., degree_days_to_sporulation=100., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='new_states_2', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=200., 
                   degree_days_to_necrosis=100., degree_days_to_sporulation=100., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='new_states_2', nreps=5)
                   
def temp39():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False,
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=80., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='new_states_3', nreps=5)
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=80., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='new_states_3', nreps=5)
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

def temp40():                 
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, proba_inf=1, age_infection=False, 
                   growth_rate=0.0006, Smin=0.02, degree_days_to_chlorosis=160., 
                   degree_days_to_necrosis=160., degree_days_to_sporulation=50., 
                   sporulating_fraction=0.1, reduction_by_rain=0.,
                   rain_events_to_empty=10, suffix='rain_33', nreps=1)
