from septo_decomposed import *
import time

start_time = time.clock()

g, recorder = run_disease(start_date = "2011-10-21 12:00:00", end_date = "2012-07-18 01:00:00", 
                variety = 'Tremie12', nplants = 5, nsect = 7,
                disc_level = 5, dir = './adel/tremie_2012_5pl_7sect', 
                sporulating_fraction = 2e-4, layer_thickness = 0.01, record = True,
                degree_days_to_chlorosis = 150., distri_chlorosis = {'mu':150., 'sigma':30.},
                age_physio_switch_senescence = 0.01, temp_min = 0.)

elapsed_time = time.clock() - start_time
print('\n')
print('----------------------------------------------')
print("Time elapsed: {} seconds".format(elapsed_time))
print('----------------------------------------------')

start_time = time.clock()

g, recorder = run_disease(start_date = "2011-10-21 12:00:00", end_date = "2012-07-18 01:00:00", 
                variety = 'Tremie12', nplants = 5, nsect = 7,
                disc_level = 5, dir = './adel/tremie_2012_5pl_7sect', 
                sporulating_fraction = 1e-4, layer_thickness = 0.01, record = False,
                degree_days_to_chlorosis = 150., distri_chlorosis = {'mu':150., 'sigma':30.},
                age_physio_switch_senescence = 0.01, temp_min = 0.)

elapsed_time = time.clock() - start_time
print('\n')
print('----------------------------------------------')
print("Time elapsed: {} seconds".format(elapsed_time))
print('----------------------------------------------')
