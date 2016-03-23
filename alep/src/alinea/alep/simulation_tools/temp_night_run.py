# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:21:59 2015

@author: ggarin
"""

from alinea.alep.simulation_tools.septo_decomposed import *
import numpy as np

              
def mon1():
    run_reps_septo(year=2011, variety='Mercia', sowing_date='10-15',
               nplants=15, sporulating_fraction=5e-5, leaf_duration=3.,
               suffix='20160119_new_wheat_inoc', nreps=10)
               
def mon2():
    run_reps_septo(year=2011, variety='Rht3', sowing_date='10-15',
               nplants=15, sporulating_fraction=5e-5, leaf_duration=3.,
               suffix='20160119_new_wheat_inoc', nreps=10)
               
def mon3():
    run_reps_septo(year=2012, variety='Tremie12', sowing_date='10-21',
               nplants=15, sporulating_fraction=1.5e-3, leaf_duration=3.,
               suffix='20160126_new_inoc_Hdrop4', nreps=10)
                             
def mon4():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
               nplants=15, sporulating_fraction=1.5e-3, leaf_duration=3.,
               suffix='20160126_new_inoc_Hdrop4', nreps=10)

from alinea.alep.simulation_tools.septo_rust_decomposed import *
def wed1():
    example_climate(years=[2003, 2012, 2013], variety='Tremie13', sowing_date='10-29',
           nplants=15, inoc_septo = 2.5e-3, inoc_rust = 300.,
           suffix='20160205_incub', nreps=3)

def wed2():
    example_climate(years=[2003, 2012, 2013], variety='Tremie13', sowing_date='10-29',
           nplants=15, inoc_septo = 2.5e-3, inoc_rust = 1e5, force_inoc_flag=True,
           suffix='20160205_incub_strong_inoc', nreps=3)
           
def wed3():
    example_climate(years=[2003, 2012, 2013], variety='Tremie13', sowing_date='10-29',
           nplants=15, inoc_septo = 2.5e-3, inoc_rust = 300.,
           x0_rust=250., latency_rust=100,
           suffix='20160205_incub_x0250_latency100', nreps=3)
           
def wed4():
    example_climate(years=[2003, 2012, 2013], variety='Tremie13', sowing_date='10-29',
           nplants=15, inoc_septo = 2.5e-3, inoc_rust = 300.,
           Smax_rust = 0.3,
           suffix='20160205_incub_Smax03', nreps=3)


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
