from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import run_disease

num_simus = range(200)

def save_simu(i_sim):
    try:
        print '-----------------------------------------------'
        print i_sim
        g, recorder = run_disease(end_date = "2011-08-31 01:00:00", sporulating_fraction = 1e-4)
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