""" Test the stability of simulation of epidemics
"""

from septo_decomposed import annual_loop_septo
from brown_rust_decomposed import annual_loop_rust
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pandas as pd
import os
from alinea.alep.disease_outputs import plot_by_leaf

def get_output_path(fungus='rust', variety='Tremie13', 
                    year = 2013, nb_pl=1, inoc=300):
    inoc = str(inoc)
    inoc = inoc.replace('.', '_')
    return './stability/'+fungus+'/'+variety.lower()+'_'+ \
            str(year)+'_'+str(nb_pl)+'pl_inoc'+inoc+'.csv'

def run_and_save_septo(variety='Tremie13', 
                       year = 2013,
                       sowing_date = '10-15',
                       sporulating_fraction = 1e-3,
                       nplants=[1, 5, 10, 20, 30, 50],
                       **kwds):
    for nb_pl in nplants:
        output_file = get_output_path(fungus='septoria', variety=variety,
                                      year=year, nb_pl=nb_pl,
                                      inoc=sporulating_fraction)
        annual_loop_septo(variety=variety, year=year, sowing_date=sowing_date,
                          sporulating_fraction=sporulating_fraction,
                          nb_plants=nb_pl, output_file=output_file, **kwds)
                          
def run_and_save_rust(variety='Tremie13', 
                       year = 2013,
                       sowing_date = '10-15',
                       density_dispersal_units = 300,
                       nplants=[1, 5, 10, 20, 30, 50],
                       **kwds):
    for nb_pl in nplants:
        output_file = get_output_path(fungus='rust', variety=variety,
                                      year=year, nb_pl=nb_pl,
                                      inoc=density_dispersal_units)
        annual_loop_rust(variety=variety, year=year, sowing_date=sowing_date,
                          density_dispersal_units=density_dispersal_units,
                          nplants=nb_pl, output_file=output_file,**kwds)
                          
def plot_stability(fungus='rust', variety='Tremie13',
                   year = 2013, sowing_date = '10-15', 
                   inoc = 300,
                   nplants = [1, 5, 10, 20, 30, 50],
                   leaves = [10, 5, 1]):
    fig, axs = plt.subplots(1,3)
    cm = plt.get_cmap('hot') 
    cNorm  = colors.Normalize(vmin=0, vmax=nplants[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    for nb_pl in nplants:
        output_file = get_output_path(fungus='rust', variety=variety,
                                      year=year, nb_pl=nb_pl, inoc=inoc)
        if os.path.exists(output_file):
            df = pd.read_csv(output_file, sep=',')
            for i, ax in enumerate(axs):
                lf = leaves[i]
                color = scalarMap.to_rgba(nb_pl)
                plot_by_leaf(df, variable='severity', xaxis='degree_days',
                             leaves=[lf], fixed_color=color, alpha=0.5, ax=ax,
                             ylims=[0, 1], legend=False)
                ax.annotate('Leaf %d' % lf, xy=(0.05, 0.95), 
                            xycoords='axes fraction', fontsize=14)
            
