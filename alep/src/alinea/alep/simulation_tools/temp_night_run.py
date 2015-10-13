# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:21:59 2015

@author: ggarin
"""

from alinea.alep.simulation_tools.septo_decomposed import *
import numpy as np

def mon1():
    run_reps_septo(year=2011, variety='Mercia', sowing_date='10-15',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=220.,
               degree_days_to_necrosis=60., degree_days_to_sporulation=50.,
               sporulating_fraction=2e-4, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=3., keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               apply_sen='incubation', age_physio_switch_senescence=100/220.,
               suffix='monday_100_220_330', nreps=10)
def mon2():
    run_reps_septo(year=2011, variety='Rht3', sowing_date='10-15',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=220.,
               degree_days_to_necrosis=60., degree_days_to_sporulation=50.,
               sporulating_fraction=2e-4, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=3., keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               apply_sen='incubation', age_physio_switch_senescence=100/220.,
               suffix='monday_100_220_330', nreps=10)
def sun3():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=220.,
               degree_days_to_necrosis=60., degree_days_to_sporulation=50.,
               sporulating_fraction=2.5e-3, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=3., keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               apply_sen='incubation', age_physio_switch_senescence=100/250.,
               suffix='sunday_100_220_330', nreps=10)
def sun4():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-29',
               nplants=15, proba_inf=1, age_infection=False,
               growth_rate=0.0006, Smin=0.03, degree_days_to_chlorosis=220.,
               degree_days_to_necrosis=60., degree_days_to_sporulation=50.,
               sporulating_fraction=7e-3, reduction_by_rain=0.,
               rain_events_to_empty=10, leaf_duration=3., keep_leaves=False,
               rh_effect=True, apply_rh='all', rh_max=35., rh_min=35.,
               apply_sen='incubation', age_physio_switch_senescence=100/220.,
               suffix='sunday_100_220_330', nreps=10)

#suffixes = ['friday_100_250_350', 'friday_100_220_330']
#from matplotlib.pyplot import savefig
##for variety, inoc in zip(['Tremie12', 'Tremie13'], [7e-3, 5e-3]):
#for variety, inoc in zip(['Tremie13'], [5e-3]):
#    for suffix in suffixes:
#        data_sim = get_aggregated_data_sim(variety=variety, nplants = 15,
#                                           sporulating_fraction=inoc, 
#                                           suffix=suffix)
#        temp_plot_by_leaf(data_sim, multiply_sev = False, xlims=[800,2200])
#        savefig(variety+'_'+suffix+'_inoc'+str(inoc)+'.png')
#        plot_by_leaf(data_sim, 'severity')
#        savefig(variety+'_'+suffix+'_inoc'+str(inoc)+'_severity.png')
