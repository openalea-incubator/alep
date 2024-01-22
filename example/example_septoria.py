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
                    leaf_area = 1., senescence_date = 5000, 
                    return_lesion=False,**kwds):
    weather = read_weather_year(2012)
    df = weather.data
    df = df.reset_index()
    df = df[df['degree_days']>0]
    df['date'] = df['datetime']
    g, leaf = get_g_and_one_leaf()
    leaf.area = leaf_area
    leaf.green_area = leaf_area
#    leaf.green_area = leaf_area-3
    septoria = plugin_septoria()
    lesion = septoria.lesion(mutable = False, group_dus = True, **kwds)
    nb_lesions = density_lesions * leaf_area
    lesion.set_position([[1] for i in range(int(nb_lesions))])
#    lesion.set_position([[1] for i in range(int(nb_lesions/2))]+[[10.] for i in range(int(nb_lesions/2))])
    leaf.lesions = [lesion]
    growth_controler = SeptoRustCompetition()
#    growth_controler = NoPriorityGrowthControl()
    surfs = []
    surfs_inc = []
    surfs_chlo = []
    surfs_nec = []
    surfs_spo = []
    surfs_empty = []
    surfs_sen = []
    tt = []
    cum_tt = 0.
    count = 0.
    for i, row in df.iterrows():
        if count == senescence_date:
#            leaf.green_area = 0.
            leaf.senesced_length = 3.
        if count < nb_steps:
            count += 1
            leaf.temperature_sequence = [row['temperature_air']]
            leaf.relative_humidity_sequence = [row['relative_humidity']]
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
            surfs_sen.append(lesion.surface_senescent)
            cum_tt += lesion.ddday
            tt.append(cum_tt)

    fig, ax = plt.subplots()
    inc = (141/255., 200/255., 2/255.)
    chlo = (1., 220/255., 63/255.)
    nec = (190/255., 128/255., 60/255.)
    spo = (228/255., 108/255., 10/255.)
    ax.plot(tt, surfs, 'b', linewidth=2)
    ax.plot(tt, surfs_inc, color=inc, linewidth=2)
    ax.plot(tt, surfs_chlo, color=chlo, linewidth=2)
    ax.plot(tt, surfs_nec, color=nec, linewidth=2)
    ax.plot(tt, surfs_spo, color=spo, linewidth=2)
    ax.plot(tt, surfs_sen, 'k')
#    ax.plot([tt[0], tt[-1]], [leaf_area, leaf_area], 'k--')
    ax.set_xlabel("Thermal time (Teff)", fontsize = 20)
    ax.set_ylabel("Surface d'une lesion (cm2)", fontsize = 20)
    ax.set_ylim([0,0.31*density_lesions*leaf_area])
    if sum(surfs_empty)>0:
        ax.plot(tt, surfs_empty, 'k--')
        labels = ['Total', 'Incubation', 'Chlorose', 'Necrose', 'Sporulant', 'Vide']
    else:
        labels = ['Total', 'Incubation', 'Chlorose', 'Necrose', 'Sporulant']
    ax.legend(labels, loc = 'best')
    ax.tick_params(axis='both', labelsize=18)
    if return_lesion==True:
        return lesion
    
def scenario_compet(nb_ddays = 1500, leaf_area = 22., 
                    deposits=[50,50,200], **kwds):
    weather = read_weather_year(2012)
    df = weather.data
    df = df.reset_index()
    df = df[df['degree_days']>1200]
    g, leaf = get_g_and_one_leaf()
    leaf.area = leaf_area
    leaf.green_area = leaf_area
    leaf.length = 1.
    leaf.senesced_length = 0.
    septoria = plugin_septoria()
    lesion = septoria.lesion(mutable = False, group_dus = True, **kwds)
    nb_lesions_init = deposits[0]
    lesion.set_position([[np.random.random()] for i in range(int(nb_lesions_init))])
    leaf.lesions = [lesion]
    growth_controler = SeptoRustCompetition()
#    growth_controler = NoPriorityGrowthControl()
    severity = []
    tot_severity = []
    severity_inc = []
    severity_chlo = []
    severity_nec = []
    severity_spo = []
    nb_lesions = []
    green_area = []
    tt = []
    count = 0
    for i, row in df.iterrows():
         if row['degree_days'] < 1060+nb_ddays:
            leaf.temperature_sequence = [row['temperature_air']]
            leaf.relative_humidity_sequence = [row['relative_humidity']]
            leaf.green_area = leaf_area - logistic(row['degree_days'], Kmax=leaf_area)
            green_area.append(leaf.green_area*100./leaf.area)
            leaf.senesced_length = logistic(row['degree_days'], Kmax=leaf.length)
            if row['degree_days'] >= 1160 and count == 0:
                count += 1
                nb_new_les = deposits[1]             
                new_les = septoria.lesion(mutable = False, group_dus = True, **kwds)
                new_les.set_position([[np.random.random()] for i in range(int(nb_new_les))])
                leaf.lesions.append(new_les)
            if row['degree_days'] >= 1400 and count == 1:
                count += 1
                nb_new_les = deposits[2]                
                new_les = septoria.lesion(mutable = False, group_dus = True, **kwds)
                new_les.set_position([[np.random.random()] for i in range(int(nb_new_les))])
                leaf.lesions.append(new_les)
            for lesion in leaf.lesions:
                lesion.update(leaf = leaf)
            growth_controler.control(g)
            surfs = sum([(l.surface_spo + l.surface_empty) for l in leaf.lesions])
            tot_surfs = sum([(l.surface) for l in leaf.lesions])
            inc_surfs = sum([l.surface_inc for l in leaf.lesions])
            chlo_surfs = sum([l.surface_chlo for l in leaf.lesions])
            nec_surfs = sum([l.surface_nec for l in leaf.lesions])
            spo_surfs = sum([l.surface_spo for l in leaf.lesions])
            severity.append(surfs*100./leaf.area)
            tot_severity.append(tot_surfs*100./leaf.area)
            severity_inc.append(inc_surfs*100./leaf.area)
            severity_chlo.append(chlo_surfs*100./leaf.area)
            severity_nec.append(nec_surfs*100./leaf.area)
            severity_spo.append(spo_surfs*100./leaf.area)
            nb_lesions.append(sum([l.nb_lesions for l in leaf.lesions]))
            tt.append(row['degree_days'])

    fig, ax = plt.subplots()
    ax.plot(tt, tot_severity, 'b')
    ax.plot(tt, severity_inc, 'g')
    ax.plot(tt, severity_chlo, 'm')
    ax.plot(tt, severity_nec, 'k')
    ax.plot(tt, severity_spo, 'r')
#    ax.plot(tt, severity, 'm')
#    ax.plot(tt, nb_lesions, 'b--')
    ax.plot(tt, green_area, 'g')
    ax.set_ylim([0,100])
    ax.set_xlabel("Thermal time (Teff)", fontsize = 18)
    ax.set_ylabel("Severity", fontsize = 18)
    labels = ['total severity', 'severity', 'nb_lesions']
    ax.legend(labels, loc = 'best')
    
def logistic(x, x0=1800, k=0.05, Kmax=22.):
    """ Calculate y for x with logistic curve """
    return Kmax / (1. + np.exp( -k * (x - x0)))
    
def gompertz(x, k = 0.003, B = 3.):
    """ Calculate y for x with gompertz curve """
    return np.exp(-B * np.exp(-k*x))
    
def linear_growth(x):
    """ Calculate y for x in linear growth in septoria lesion model. """
    if x <= 220:
        return x*0.03/220.
    elif x <= 670:
        return 0.03+0.0006*(x-220)
    else:
        return 0.3
        
def linear_growth_list(xs):
    return np.array([linear_growth(x) for x in xs])
    
def optimize_gompertz_on_linear():
    from scipy.optimize import minimize
    def rmse_gomp(params = (0.003, 3.)):
        ref = linear_growth_list(x)
        y = gompertz(x, k=params[0], B=params[1])*0.3
        rmse = np.sqrt((ref-y)**2).mean()
        return rmse
    x = np.arange(1000)
    res = minimize(rmse_gomp, (0.003, 3.))
    k, B = res.x
    plot(x, linear_growth_list(x))
    plot(x, gompertz(x, k, B)*0.3)
    return k, B
    
def optimize_logistic_on_linear():
    from scipy.optimize import minimize
    def rmse_logi(params = (300., 0.05)):
        ref = linear_growth_list(x)
        y = logistic(x, x0=params[0], k=params[1], Kmax=0.3)
        rmse = np.sqrt((ref-y)**2).mean()
        return rmse
    x = np.arange(1000)
    res = minimize(rmse_logi, (300., 0.05))
    x0, k = res.x
    plot(x, linear_growth_list(x))
    plot(x, logistic(x, x0=x0, k=k, Kmax=0.3))
    return x0, k

if __name__ == '__main__':
    example_surface()
    plt.show()