""" Examples to test consistency of septoria model.
"""
from alinea.echap.weather_data import read_weather_year
from alinea.alep.septoria import plugin_septoria
from alinea.alep.disease_outputs import AdelSeptoRecorder, plot_by_leaf, conf_int
from alinea.alep.growth_control import SeptoRustCompetition
from alinea.alep.inoculation import AirborneContamination
from alinea.alep.protocol import infect, update, disperse, external_contamination
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum
from alinea.popdrops.alep_interface import PopDropsSoilContamination, PopDropsEmission, PopDropsTransport
from alinea.alep.growth_control import SeptoRustCompetition, NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.echap.architectural_reconstructions import (EchapReconstructions,
                                                        echap_reconstructions,
                                                        reconstruction_parameters,
                                                        soisson_reconstruction)
from alinea.adel.data_samples import adel_two_metamers_stand
from alinea.alep.architecture import get_leaves, set_properties
from alinea.adel.newmtg import adel_labels
from alinea.alep.simulation_tools.simulation_tools import get_weather
from alinea.alep.alep_time_control import CustomIterWithDelays
from alinea.alep.alep_weather import plot_wetness_and_temp
from alinea.astk.TimeControl import (time_filter, IterWithDelays,
                                     thermal_time_filter, DegreeDayModel,
                                     time_control)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from datetime import datetime

def get_small_g():
    g, domain_area, domain, convunit = adel_two_metamers_stand(leaf_sectors = 1,
                                                               density = 350.,
                                                               interleaf = 10.,
                                                               leaf_length = 20,
                                                               leaf_width = 1, 
                                                               Einc = 0)
    return g

def get_g_and_one_leaf():
    g = get_small_g()
    leaves = get_leaves(g)
    return g, g.node(leaves[0])

def get_one_leaf():
    g = get_small_g()
    leaves = get_leaves(g)
    return g.node(leaves[0])

def example_surface(nb_steps = 4500, density_lesions = 1, with_compet = False,
                    leaf_area = 1., **kwds):
    weather = read_weather_year(2012)
    df = weather.data
    df = df.reset_index()
    df = df[df['degree_days']>0]
    g, leaf = get_g_and_one_leaf()
    leaf.area = leaf_area
    leaf.green_area = leaf_area
    septoria = plugin_septoria()
    lesion = septoria.lesion(mutable = False, group_dus = True, **kwds)
    nb_lesions = density_lesions * leaf_area
    lesion.set_position([[1] for i in range(int(nb_lesions))])
    leaf.lesions = [lesion]
    growth_controler = SeptoRustCompetition()
    surfs = []
    surfs_inc = []
    surfs_chlo = []
    surfs_nec = []
    surfs_spo = []
    surfs_empty = []
    tt = []
    cum_tt = 0.
    count = 0.
    for i, row in df.iterrows():
         if count < nb_steps:
            count += 1
            leaf.temperature_sequence = [row['temperature_air']]
            lesion.update(leaf = leaf)
            if with_compet == False:
                lesion.control_growth(growth_offer = lesion.growth_demand)
            else:
                growth_controler.control(g)
            surfs.append(lesion.surface)
            surfs_inc.append(lesion.surface_inc)
            surfs_chlo.append(lesion.surface_chlo)
            surfs_nec.append(lesion.surface_nec)
            surfs_spo.append(lesion.surface_spo)
            surfs_empty.append(lesion.surface_empty)
            cum_tt += lesion.ddday
            tt.append(cum_tt)

    fig, ax = plt.subplots()
    ax.plot(tt, surfs, 'b')
    ax.plot(tt, surfs_inc, 'g')
    ax.plot(tt, surfs_chlo, 'm')
    ax.plot(tt, surfs_nec, 'k')
    ax.plot(tt, surfs_spo, 'r')
    ax.set_xlabel("Thermal time (Teff)", fontsize = 18)
    ax.set_ylabel("Surface d'une lesion (cm2)", fontsize = 18)
    if sum(surfs_empty)>0:
        ax.plot(tt, surfs_empty, 'k--')
        labels = ['Total', 'Incubation', 'Chlorose', 'Necrose', 'Sporulant', 'Vide']
    else:
        labels = ['Total', 'Incubation', 'Chlorose', 'Necrose', 'Sporulant']
    ax.legend(labels, loc = 'best')