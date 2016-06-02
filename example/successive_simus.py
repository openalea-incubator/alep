from septo_decomposed import *

def simu(nb_plants):
    dir_path = './adel/tremie_2012_'+str(nb_plants)+'pl_7sect'
    g, recorder = run_disease(start_date="2011-10-21 12:00:00", end_date="2012-08-01 00:00:00",
                                variety = 'Tremie12', nplants=nb_plants, nsect=7,
                                dir = dir_path, sporulating_fraction = 1.5e-4,
                                temp_min = 0., degree_days_to_chlorosis = 150,
                                age_physio_switch_senescence = 0.01,
                                distri_chlorosis = {'mu':150., 'sigma':30.})
    stored_rec = './tremie/recorder_stabtest_tremie12_'+str(nb_plants)+'pl.pckl'
    f_rec = open(stored_rec, 'w')
    pickle.dump(recorder, f_rec)
    f_rec.close()
    del recorder
    del g

#for nb_pl in [1, 5, 10, 15, 20, 30, 50, 100]:
 #   simu(nb_pl)
simu(10)
