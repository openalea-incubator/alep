from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import run_disease
try:
    import cPickle as pickle
except:
    import pickle

num_simus = range(200)
adel, domain, domain_area, convUnit, weather, seq, rain_timing, canopy_timing, septo_timing = setup(
    start_date = "2010-10-15 12:00:00", end_date = "2011-08-31 01:00:00", nplants = 3, nsect = 5, disc_level = 20)

def save_simu(i_sim):
    try:
        print '-----------------------------------------------'
        print i_sim
        g, recorder = run_disease(sporulating_fraction = 1e-4, adel = adel, domain = domain,
                                  domain_area = domain_area, convUnit = convUnit, weather = weather,
                                  seq = seq, rain_timing = rain_timing,
                                  canopy_timing = canopy_timing, septo_timing = septo_timing)
        stored_rec = '.\mercia\\recorder_'+str(i_sim)+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
        print 'simulation succeded'
    except:
        print 'simulation failed'
        
if __name__ == '__main__':
    nb_cpu = cpu_count()
    pymap(save_simu, num_simus, nb_cpu-1)