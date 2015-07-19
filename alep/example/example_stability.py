""" Test the stability of simulation of epidemics
"""

from septo_decomposed import annual_loop_septo
from brown_rust_decomposed import annual_loop_rust
import matplotlib.pyplot as plt

def run_and_save_septo(variety='Tremie13', 
                       year = 2013,
                       sowing_date = '10-15',
                       sporulating_fraction = 1e-3,
                       nplants=[1, 5, 10, 20, 30, 50, 100]):
    for nb_pl in nplants:
        output_file = './stability/septoria/'+variety.lower()+'_'+\
                    str(year)+'_'+str(nb_pl)+'pl.csv'
        annual_loop_septo(variety=variety, year=year, sowing_date=sowing_date,
                          sporulating_fraction=sporulating_fraction,
                          nb_plants=nb_pl, output_file=output_file)
                          
def run_and_save_rust(variety='Tremie13', 
                       year = 2013,
                       sowing_date = '10-15',
                       density_dispersal_units = 500,
                       nplants=[1, 5, 10, 20, 30, 50, 100]):
    for nb_pl in nplants:
        output_file = './stability/'+variety.lower()+'_'+\
                    str(year)+'_'+str(nb_pl)+'pl.csv'
        annual_loop_rust(variety=variety, year=year, sowing_date=sowing_date,
                          density_dispersal_units=density_dispersal_units,
                          nb_plants=nb_pl, output_file=output_file)
