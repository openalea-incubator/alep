""" Script to work on performance of typical septoria annual loop """
from septo_decomposed import make_canopy, run_disease
import cProfile
import pstats

# make_canopy(start_date = "2011-10-21 12:00:00", end_date = "2012-07-18 01:00:00",
            # variety = 'Tremie12', nplants = 5, nsect = 7, disc_level = 5, 
            # dir = './adel/tremie_2012_5pl_7sect', reset_reconst = True)
            
def annual_loop(mutable = False):
    if mutable == True:
        # Signifie que chaque instance de lesion possede ses propres parametres
        # (Bon pour simu de variabilite, mais couteux en memoire)
        distri_chlorosis = {'mu':150., 'sigma':30}
    else:
        # Signifie que toutes les lesions ont les memes parametres
        distri_chlorosis = None
    g, recorder = run_disease(start_date = "2011-10-21 12:00:00", end_date = "2012-07-18 01:00:00", 
                                variety = 'Tremie12', nplants = 5, nsect = 7, disc_level = 5, 
                                dir = './adel/tremie_2012_5pl_7sect', sporulating_fraction = 1e-4,
                                distri_chlorosis = distri_chlorosis, degree_days_to_chlorosis = 150.)
                                
def stat_profiler(call='annual_loop()', filename = 'septo_stats.pstats'):
    cProfile.run(call, filename)
    
if __name__ == '__main__':
    stat_profiler(call='annual_loop(mutable = True)', filename = 'septo_stats_mutable_params.pstats')
    stat_profiler(call='annual_loop(mutable = False)', filename = 'septo_stats_fixed_params.pstats')