""" Test the stability of simulation of epidemics
"""

from septo_decomposed import annual_loop_septo
from brown_rust_decomposed import annual_loop_rust
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pandas as pd
import numpy as np
import os
from alinea.alep.disease_outputs import (plot_by_leaf, conf_int,
                                         get_synthetic_outputs_by_leaf)

# Working with full recorder ##################################################
def get_output_path_full_recorder(fungus='rust', variety='Tremie13', 
                                  year = 2013, nb_pl=1, inoc=300, rep=0):
    inoc = str(inoc)
    inoc = inoc.replace('.', '_')
    return './stability/'+fungus+'/'+variety.lower()+'_'+ \
            str(year)+'_'+str(nb_pl)+'pl_inoc'+inoc+'_rep_'+str(rep)+'.csv'

def run_and_save_septo_full_recorder(variety='Tremie13', 
                                   year = 2013,
                                   sowing_date = '10-15',
                                   sporulating_fraction = 1e-3,
                                   nplants=[1, 5, 10, 20, 30, 50],
                                   nreps = 5,
                                   **kwds):
    for rep in range(nreps):
        for nb_pl in nplants:
            output_file = get_output_path_full_recorder(fungus='septoria',
                                                        variety=variety,
                                                        year=year, nb_pl=nb_pl,
                                                        inoc=sporulating_fraction,
                                                        rep = rep)
            annual_loop_septo(variety=variety, year=year, sowing_date=sowing_date,
                              sporulating_fraction=sporulating_fraction,
                              nplants=nb_pl, output_file=output_file, **kwds)
                          
def run_and_save_rust_full_recorder(variety='Tremie13', 
                                   year = 2013,
                                   sowing_date = '10-15',
                                   density_dispersal_units = 300,
                                   nplants=[1, 5, 10, 20, 30, 50],
                                   nreps = 5,
                                   **kwds):
    for rep in range(nreps):    
        for nb_pl in nplants:
            output_file = get_output_path_full_recorder(fungus='rust',
                                                        variety=variety,
                                                        year=year, nb_pl=nb_pl,
                                                        inoc=density_dispersal_units,
                                                        rep=rep)
            annual_loop_rust(variety=variety, year=year, sowing_date=sowing_date,
                              density_dispersal_units=density_dispersal_units,
                              nplants=nb_pl, output_file=output_file,**kwds)
                          
def plot_stability_full_recorder(fungus='rust', variety='Tremie13',
                                   year = 2013, sowing_date = '10-15', 
                                   inoc = 300,
                                   nplants = [1, 5, 10, 20, 30, 50],
                                   leaves = [10, 5, 1],
                                   nreps = 5, cmap='hot'):
    fig, axs = plt.subplots(1,3)
    cm = plt.get_cmap(cmap) 
    cNorm  = colors.Normalize(vmin=0, vmax=nplants[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    for rep in range(nreps):
        for nb_pl in nplants:
            output_file = get_output_path_full_recorder(fungus=fungus,
                                                        variety=variety,
                                                        year=year, nb_pl=nb_pl,
                                                        inoc=inoc,
                                                        rep=rep)
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

# Working with audpc only #####################################################
def get_output_path_synthetic(fungus='rust', variety='Tremie13', 
                          year = 2013, inoc=300):
    inoc = str(inoc)
    inoc = inoc.replace('.', '_')
    return './stability/'+fungus+'/'+variety.lower()+'_'+ \
            str(year)+'_pl_inoc'+inoc+'_synthetic.csv'

def run_and_save_septo_synthetic(variety = 'Tremie13', 
                                 year = 2013,
                                 sowing_date = '10-15',
                                 sporulating_fraction = 1e-3,
                                 nplants = [1, 5, 10, 15, 20, 25, 30, 40, 50],
                                 nb_reps = 5, **kwds):
    output_file = get_output_path_synthetic(fungus='septoria',
                                            variety=variety,
                                            year=year, 
                                            inoc=sporulating_fraction)
    df_out = pd.DataFrame()
    for nb_pl in nplants:
        for rep in nb_reps:
            g, reco = annual_loop_septo(variety=variety, year=year,
                                        sowing_date=sowing_date,
                                        sporulating_fraction=sporulating_fraction,
                                        nplants=nb_pl, 
                                        output_file=output_file, **kwds)
            df = get_synthetic_outputs_by_leaf(reco.data)
            df['rep'] = rep
            df['nb_plants'] = nb_pl
            df_out = pd.concat([df_out, df])
    df_out.to_csv(output_file, sep = ',')

def run_and_save_rust_synthetic(variety = 'Tremie13', 
                             year = 2013,
                             sowing_date = '10-15',
                             density_dispersal_units = 300,
                             nplants = [1, 5, 10, 15, 20, 25, 30, 40],
                             nb_reps = 5, **kwds):
    output_file = get_output_path_synthetic(fungus='rust',
                                            variety=variety,
                                            year=year, 
                                            inoc=density_dispersal_units)
    df_out = pd.DataFrame()
    for nb_pl in nplants:
        for rep in range(nb_reps):
            g, reco = annual_loop_rust(variety=variety, year=year,
                                        sowing_date=sowing_date,
                                        density_dispersal_units=density_dispersal_units,
                                        nplants=nb_pl, **kwds)
            df = get_synthetic_outputs_by_leaf(reco.data)
            df['rep'] = rep
            df['nb_plants'] = nb_pl
            df_out = pd.concat([df_out, df])
    df_out.to_csv(output_file, sep = ',')
    
def plot_stability_synthetic(fungus='rust', 
                             variety = 'Tremie13',
                             year = 2013,
                             inoc = 300,
                             nplants = [1, 5, 10, 15, 20, 25, 30, 40],
                             nb_reps = 5, leaves = [10, 5, 1], 
                             variable = 'normalized_audpc', 
                             ylabel = 'Normalized AUDPC', 
                             ylims = None):
    fig, axs = plt.subplots(1,3, figsize=(16,10))
    output_file = get_output_path_synthetic(fungus=fungus,
                                            variety=variety,
                                            year=year, 
                                            inoc=inoc)
    df = pd.read_csv(output_file, sep=',')
    for i, lf in enumerate(leaves):
        ax = axs[i]
        df_lf = df[df['num_leaf_top']==lf]
        df_mean = df_lf.groupby('nb_plants').agg(np.mean)
        df_conf = df_lf.groupby('nb_plants').agg(conf_int)
        ax.errorbar(df_mean.index, df_mean[variable], yerr=df_conf[variable],
                    linestyle='', marker='o')
        ax.set_xlabel('Number of plants', fontsize=16)
        if ax==axs[0]:
            ax.set_ylabel(ylabel, fontsize=16)
        ax.set_xlim([0, max(nplants)+5])
        if ylims is not None:
            ax.set_ylim(ylims)
        ax.annotate('Leaf %d' % lf, xy=(0.05, 0.95), 
                    xycoords='axes fraction', fontsize=14)
