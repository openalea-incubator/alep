""" Test the stability of simulation of epidemics
"""

from alinea.alep.simulation_tools.septo_decomposed import annual_loop_septo
from alinea.alep.simulation_tools.brown_rust_decomposed import annual_loop_rust
from alinea.alep.simulation_tools.simulation_tools import count_available_canopies
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pandas as pd
import numpy as np
import random as rd
import os
from alinea.alep.disease_outputs import (plot_by_leaf, conf_int,
                                         get_synthetic_outputs_by_leaf)

import collections
def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)
    
def get_iterable(nplants=[1, 5, 10, 20, 30, 50], nsect = 7, 
                 layer_thickness = 0.01):
    kwargs = locals()
    iterable = {k:v for k,v in kwargs.items() if is_iterable(v)}
    assert len(iterable)==1,"Must provide one and only iterable as argument"
    return iterable, nplants, nsect, layer_thickness

def get_value_iteration(it, iterable, nb_pl, ns, dh):
    if list(iterable.keys())[0]=='nplants':
        nb_pl = it
    elif list(iterable.keys())[0]=='nsect':
        ns = it
    elif list(iterable.keys())[0]=='layer_thickness':
        dh = it
    return nb_pl, ns, dh

# Working with full recorder ##################################################
def get_output_path_full_recorder(fungus='rust', variety='Tremie13', 
                                  year = 2013, nb_pl=1, nsect=5,
                                  dh=0.01, inoc=300, rep=0):
    inoc = str(inoc)
    inoc = inoc.replace('.', '_')
    dh = str(dh)
    dh = dh.replace('.', '_')
    return './stability/'+fungus+'/'+variety.lower()+'_'+ \
            str(year)+'_'+str(nb_pl)+'pl_'+str(nsect)+'sect_dh_'+ \
            dh+'_inoc'+inoc+'_rep_'+str(rep)+'.csv'

def run_and_save_septo_full_recorder(variety='Tremie13', 
                                   year = 2013,
                                   sowing_date = '10-15',
                                   sporulating_fraction = 1e-3,
                                   nplants=[1, 5, 10, 20, 30, 50],
                                   nreps = 5,
                                   nsect = 7,
                                   layer_thickness = 0.01,
                                   **kwds):
    iterable, nb_pl, ns, dh = get_iterable(nplants, nsect, layer_thickness)
    for rep in range(nreps):
        for it in list(iterable.values())[0]:
            nb_pl, ns, dh = get_value_iteration(it, iterable, nb_pl, ns, dh)
            output_file = get_output_path_full_recorder(fungus='septoria',
                                                        variety=variety,
                                                        year=year, nb_pl=nb_pl,
                                                        nsect = ns, dh = dh,
                                                        inoc=sporulating_fraction,
                                                        rep = rep)
            annual_loop_septo(variety=variety, year=year, sowing_date=sowing_date,
                              sporulating_fraction=sporulating_fraction,
                              nplants=nb_pl, nsect=ns, layer_thickness=dh,
                              output_file=output_file, **kwds)
                          
def run_and_save_rust_full_recorder(variety='Tremie13', 
                                   year = 2013,
                                   sowing_date = '10-15',
                                   density_dispersal_units = 300,
                                   nplants=[1, 5, 10, 20, 30, 50],
                                   nreps = 5,
                                   nsect = 7,
                                   layer_thickness = 1.,
                                   **kwds):
    iterable, nb_pl, ns, dh = get_iterable(nplants, nsect, layer_thickness)
    for rep in range(nreps):
        for it in list(iterable.values())[0]:
            nb_pl, ns, dh = get_value_iteration(it, iterable, nb_pl, ns, dh)
            output_file = get_output_path_full_recorder(fungus='rust',
                                                        variety=variety,
                                                        year=year, nb_pl=nb_pl,
                                                        nsect = ns, dh = dh,
                                                        inoc=density_dispersal_units,
                                                        rep=rep)
            annual_loop_rust(variety=variety, year=year, sowing_date=sowing_date,
                              density_dispersal_units=density_dispersal_units,
                              nplants=nb_pl, nsect=ns, layer_thickness=dh,
                              output_file=output_file,**kwds)
                          
def plot_stability_full_recorder(fungus='rust', variety='Tremie13',
                                   year = 2013, sowing_date = '10-15', 
                                   inoc = 300,
                                   nplants = [1, 5, 10, 20, 30, 50],
                                   nsect = 7,
                                   layer_thickness = 1.,
                                   leaves = [10, 5, 1],
                                   nreps = 5, cmap='hot'):
    iterable, nb_pl, ns, dh = get_iterable(nplants, nsect, layer_thickness)
    fig, axs = plt.subplots(1,3)
    cm = plt.get_cmap(cmap) 
    cNorm  = colors.Normalize(vmin=0, vmax=nplants[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    for rep in range(nreps):
        for it in list(iterable.values())[0]:
            nb_pl, ns, dh = get_value_iteration(it, iterable, nb_pl, ns, dh)
            output_file = get_output_path_full_recorder(fungus=fungus,
                                                        variety=variety,
                                                        year=year, nb_pl=nb_pl,
                                                        nsect = ns, dh = dh,
                                                        inoc=inoc, rep=rep)
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
                              year = 2013, inoc=300, iterable='nplants'):
    inoc = str(inoc)
    inoc = inoc.replace('.', '_')
    return './stability/'+fungus+'/'+variety.lower()+'_'+ \
            str(year)+'_inoc'+inoc+'_synthetic_'+iterable+'.csv'

def get_rep_wheats(variety = 'Tremie13', year = 2013, 
                    nplants = 15, nsect = 7, nreps = 5):
    nb_can = count_available_canopies(year, variety, nplants, nsect)
    if nb_can >= nreps:
        rep_wheats = rd.sample(list(range(nb_can)), nreps)
    else:
        rep_wheats = [None for rep in range(nreps)]
    return rep_wheats
            
def run_and_save_septo_synthetic(variety = 'Tremie13', 
                                 year = 2013,
                                 sowing_date = '10-15',
                                 sporulating_fraction = 1e-3,
                                 nplants = [1, 5, 10, 15, 20, 25, 30, 40],
                                 nb_reps = 5, nsect = 7,
                                 layer_thickness = 0.01, **kwds):
    iterable, nb_pl, ns, dh = get_iterable(nplants, nsect, layer_thickness)
    output_file = get_output_path_synthetic(fungus='septoria',
                                            variety=variety,
                                            year=year, 
                                            inoc=sporulating_fraction,
                                            iterable=list(iterable.keys())[0])
    rep_wheats = get_rep_wheats(variety, year, nplants, nsect, nreps)
    df_out = pd.DataFrame()
    for it in list(iterable.values())[0]:
        nb_pl, ns, dh = get_value_iteration(it, iterable, nb_pl, ns, dh)
        for irep, rep_wheat in enumerate(rep_wheats):
            g, reco = annual_loop_septo(variety=variety, year=year,
                                        sowing_date=sowing_date,
                                        sporulating_fraction=sporulating_fraction,
                                        nplants=nb_pl,
                                        nsect = ns, dh = dh,
                                        rep_wheat=rep_wheat, **kwds)
            df = get_synthetic_outputs_by_leaf(reco.data)
            df['rep'] = irep
            df['nb_plants'] = nb_pl
            df['nb_sects'] = ns
            df['layer_thickness'] = dh
            df_out = pd.concat([df_out, df])
    df_out.to_csv(output_file, sep = ',')

def run_and_save_rust_synthetic(variety = 'Tremie13', 
                             year = 2013,
                             sowing_date = '10-15',
                             density_dispersal_units = 300,
                             nplants = [1, 5, 10, 15, 20, 25, 30, 40],
                             nb_reps = 5, nsect = 7, layer_thickness = 1., 
                             **kwds):
    iterable, nb_pl, ns, dh = get_iterable(nplants, nsect, layer_thickness)
    output_file = get_output_path_synthetic(fungus='rust',
                                            variety=variety,
                                            year=year, 
                                            inoc=density_dispersal_units,
                                            iterable=list(iterable.keys())[0])
    rep_wheats = get_rep_wheats(variety, year, nplants, nsect, nreps)
    df_out = pd.DataFrame()
    for it in list(iterable.values())[0]:
        nb_pl, ns, dh = get_value_iteration(it, iterable, nb_pl, ns, dh)
        for irep, rep_wheat in enumerate(rep_wheats):
            g, reco = annual_loop_rust(variety=variety, year=year,
                                        sowing_date=sowing_date,
                                        density_dispersal_units=density_dispersal_units,
                                        nplants=nb_pl, nsect = ns, dh = dh,
                                        rep_wheat=rep_wheat, **kwds)
            df = get_synthetic_outputs_by_leaf(reco.data)
            df['rep'] = irep
            df['nb_plants'] = nb_pl
            df['nb_sects'] = ns
            df['layer_thickness'] = dh
            df_out = pd.concat([df_out, df])
    df_out.to_csv(output_file, sep = ',')
    
def plot_stability_synthetic(fungus='rust', 
                             variety = 'Tremie13',
                             year = 2013,
                             inoc = 300,
                             nplants = [1, 5, 10, 15, 20, 25, 30, 40],
                             nb_reps = 5, nsect = 7, layer_thickness = 1.,
                             leaves = [10, 5, 1], 
                             variable = 'normalized_audpc', 
                             ylabel = 'Normalized AUDPC', 
                             ylims = None):
    iterable, nb_pl, ns, dh = get_iterable(nplants, nsect, layer_thickness)
    fig, axs = plt.subplots(1,3, figsize=(16,6))
    output_file = get_output_path_synthetic(fungus=fungus,
                                            variety=variety,
                                            year=year, 
                                            inoc=inoc,
                                            iterable=list(iterable.keys())[0])
    df = pd.read_csv(output_file, sep=',')
    if list(iterable.keys())[0]=='nplants':
        group_by = 'nb_plants'
        xlabel = 'Number of plants'
    elif list(iterable.keys())[0]=='nsect':
        group_by = 'nb_sects'
        xlabel = 'Number of sectors'
    elif list(iterable.keys())[0]=='layer_thickness':
        group_by = 'layer_thickness'
        xlabel = 'Layer thickness'
    
    for i, lf in enumerate(leaves):
        ax = axs[i]
        df_lf = df[df['num_leaf_top']==lf]
        df_mean = df_lf.groupby(group_by).agg(np.mean)
        df_conf = df_lf.groupby(group_by).agg(conf_int)
        ax.errorbar(df_mean.index, df_mean[variable].values, yerr=df_conf[variable].values,
                    linestyle='', marker='o')
        ax.set_xlabel(xlabel, fontsize=16)
        if ax==axs[0]:
            ax.set_ylabel(ylabel, fontsize=16)
        ax.set_xlim([0, max(list(iterable.values())[0])*1.1])
        if ylims is not None:
            ax.set_ylim(ylims)
        ax.annotate('Leaf %d' % lf, xy=(0.05, 0.95), 
                    xycoords='axes fraction', fontsize=14)
