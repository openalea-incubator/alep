""" Test the stability of simulation of epidemics
"""

from septo_decomposed import annual_loop_septo
from brown_rust_decomposed import annual_loop_rust
import matplotlib.pyplot as plt
import pandas as pd
import os
from alinea.alep.disease_outputs import plot_by_leaf

def get_output_path(fungus='rust', variety='Tremie13', 
                    year = 2013, nb_pl=1):
    return './stability/'+fungus+'/'+variety.lower()+'_'+ \
            str(year)+'_'+str(nb_pl)+'pl.csv'

def run_and_save_septo(variety='Tremie13', 
                       year = 2013,
                       sowing_date = '10-15',
                       sporulating_fraction = 1e-3,
                       nplants=[1, 5, 10, 20, 30, 50, 100]):
    for nb_pl in nplants:
        output_file = get_output_path(fungus='septoria', variety=variety,
                                      year=year, nb_pl=nb_pl)
        annual_loop_septo(variety=variety, year=year, sowing_date=sowing_date,
                          sporulating_fraction=sporulating_fraction,
                          nb_plants=nb_pl, output_file=output_file)
                          
def run_and_save_rust(variety='Tremie13', 
                       year = 2013,
                       sowing_date = '10-15',
                       density_dispersal_units = 300,
                       nplants=[1, 5, 10, 20, 30, 50, 100]):
    for nb_pl in nplants:
        output_file = get_output_path(fungus='rust', variety=variety,
                                      year=year, nb_pl=nb_pl)
        annual_loop_rust(variety=variety, year=year, sowing_date=sowing_date,
                          density_dispersal_units=density_dispersal_units,
                          nplants=nb_pl, output_file=output_file)
                          
def plot_stability(fungus='rust', variety='Tremie13',
                   year = 2013, sowing_date = '10-15',  
                   nplants = [1, 5, 10, 20, 30, 50, 100],
                   leaves = [10, 5, 1]):
    fig, axs = plt.subplots(3,1)
    for nb_pl in nplants:
        output_file = get_output_path(fungus='rust', variety=variety,
                                      year=year, nb_pl=nb_pl)
        if os.path.exists(output_file):
            df = pd.read_csv(output_file, sep=',')
            for i, ax in enumerate(axs):
                lf = leaves[i]
                plot_by_leaf(df, variable='severity', xaxis='degree_days',
                             leaves=[lf], fixed_color='k', alpha=0.3, ax=ax)
                ax.annotate('Leaf %d' % lf, xy=(0.05, 0.95), 
                            xycoords='axes fraction', fontsize=14)
            
