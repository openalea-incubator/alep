# -*- coding: latin1 -*-
""" Examples to test consistency of brown rust model.
"""
from alinea.echap.weather_data import read_weather_year
from alinea.alep.brown_rust import BrownRustFungus
from alinea.alep.septoria_age_physio import SeptoriaFungus
from alinea.alep.disease_outputs import BrownRustRecorder, plot_by_leaf
from alinea.alep.growth_control import (NoPriorityGrowthControl,
                                        GeometricCircleCompetition,
                                        SeptoRustCompetition)
from alinea.alep.inoculation import AirborneContamination
from alinea.alep.protocol import infect, update
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.dispersal_transport import BrownRustDispersal
from alinea.echap.architectural_reconstructions import EchapReconstructions
from alinea.adel.data_samples import adel_two_metamers_stand
from alinea.alep.architecture import get_leaves
from alinea.alep.architecture import set_properties
from alinea.adel.newmtg import adel_labels
from septo_decomposed import get_weather
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

def example_surface(**kwds):
    weather = read_weather_year(2012)
    df = weather.data
    df = df.reset_index()
    df = df[df['degree_days']>0]
    leaf = get_one_leaf()
    brown_rust = BrownRustFungus()
    lesion = brown_rust.lesion(mutable = False, group_dus = False, **kwds)
    surfs = []
    surfs_chlo = []
    surfs_spo = []
    surfs_sink = []
    surfs_empty = []
    tt = []
    cum_tt = 0.
    count = 0.
    for i, row in df.iterrows():
         if count < 4500:
            count += 1
            leaf.temperature_sequence = [row['temperature_air']]
            lesion.update(leaf = leaf)
            lesion.control_growth(growth_offer = lesion.growth_demand)
            surfs.append(lesion.surface)
            surfs_chlo.append(lesion.surface_chlo)
            surfs_spo.append(lesion.surface_spo)
            surfs_sink.append(lesion.surface_sink)
            surfs_empty.append(lesion.surface_empty)
            cum_tt += lesion.delta_thermal_time_growth(leaf_temperature = [row['temperature_air']])
            tt.append(cum_tt)

    fig, ax = plt.subplots()
    ax.plot(tt, surfs, 'b')
    ax.plot(tt, surfs_sink, 'k')
    ax.plot(tt, surfs_chlo, 'r')
    ax.plot(tt, surfs_spo, 'g')
    ax.plot(tt, surfs_empty, 'y')
    ax.set_xlabel("Thermal time (Teff)", fontsize = 18)
    ax.set_ylabel("Surface of 1 lesion (cm2)", fontsize = 18)
    ax.legend(['sink + Chlo + Spo + Nec', 
               'sinking', 'Chlorotic', 'Sporulating', 'Necrotic'], loc = 'center right')

def example_surface_competition(date_competition = 500, **kwds):
    weather = read_weather_year(2012)
    df = weather.data
    df = df.reset_index()
    df = df[df['degree_days']>0]
    leaf = get_one_leaf()
    brown_rust = BrownRustFungus()
    lesion = brown_rust.lesion(mutable = False, group_dus = False, **kwds)
    surfs = []
    surfs_chlo = []
    surfs_spo = []
    surfs_sink = []
    surfs_empty = []
    tt = []
    cum_tt = 0.
    count = 0.
    for i, row in df.iterrows():
        if count < 3000:
            count += 1
            leaf.temperature_sequence = [row['temperature_air']]
            lesion.update(leaf = leaf)
            if count < date_competition:
                lesion.control_growth(growth_offer = lesion.growth_demand)
            else:
                lesion.control_growth(growth_offer = 0.)
            surfs.append(lesion.surface)
            surfs_chlo.append(lesion.surface_chlo)
            surfs_spo.append(lesion.surface_spo)
            surfs_sink.append(lesion.surface_sink)
            surfs_empty.append(lesion.surface_empty)
            cum_tt += lesion.delta_thermal_time_growth(leaf_temperature = [row['temperature_air']])
            tt.append(cum_tt)

    fig, ax = plt.subplots()
    ax.plot(tt, surfs, 'b')
    ax.plot(tt, surfs_sink, 'k')
    ax.plot(tt, surfs_chlo, 'r')
    ax.plot(tt, surfs_spo, 'g')
    ax.plot(tt, surfs_empty, 'y')
    ax.set_xlabel("Thermal time (Teff)", fontsize = 18)
    ax.set_ylabel("Surface of 1 lesion (cm2)", fontsize = 18)
    ax.legend(['sink + Chlo + Spo + Nec', 
               'sinking', 'Chlorotic', 'Sporulating', 'Necrotic'], loc = 'best')
               
def example_geom_competition(density = 30., **kwds):
    g, leaf = get_g_and_one_leaf()
    leaf.area = 10.
    leaf.green_area = 10.    
    
    brown_rust = BrownRustFungus()
#    growth_controler = NoPriorityGrowthControl()
    growth_controler = GeometricCircleCompetition()
    
    df_temp = pd.DataFrame([24.], 
                            columns = ['temp'])
    leaf.temperature_sequence = df_temp['temp']
    les = brown_rust.lesion(mutable = False, group_dus = True,**kwds)
    les.set_position([[1] for i in range(int(density*leaf.area))])
    leaf.lesions = [les]
    sev = []
    for i in range(500):
        les.update(leaf = leaf)
        growth_controler.control(g)
        sev.append(les.surface/leaf.green_area)
    return sev
    
def example_density_robert_2002(**kwds):
    brown_rust = BrownRustFungus()
    growth_controler = NoPriorityGrowthControl()

    def example_single_density(d = 40):
        g, leaf = get_g_and_one_leaf()
        leaf.area = 10.
        leaf.green_area = 10.
        df_temp = pd.DataFrame([24.], 
                                columns = ['temp'])
        leaf.temperature_sequence = df_temp['temp']
        les = brown_rust.lesion(mutable = False, group_dus = True,**kwds)
        les.set_position([[1] for i in range(int(d*leaf.area))])
        leaf.lesions = [les]
        for i in range(350):
            les.update(leaf = leaf)
            growth_controler.control(g)
        surfs = les.surface / les.nb_lesions
        surfs_spo = les.surface_spo / les.nb_lesions
        return surfs, surfs_spo

    filename = 'rust_density.csv'
    df_obs = pd.read_csv(filename, sep = ';')
    df_sim = pd.DataFrame(columns = ['surf_mean', 'surf_spo_mean'])
    # x = np.sort(np.concatenate([df_obs['density'], np.arange(1, 41, 2)]))
    x = np.sort(np.concatenate([df_obs['density'], np.arange(1,1)]))
    for d in x:
        df_sim.loc[len(df_sim)+1] = example_single_density(d)
    df_sim['density'] = x

    rmse = np.sqrt(((df_sim[df_sim['density'].isin(df_obs['density'])]['surf_spo_mean'] -
                        df_obs['surface']) ** 2).mean())

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(df_sim.density, df_sim.surf_spo_mean)
    ax.plot(df_obs.density, df_obs.surface, 'ro')
#    ax.annotate('Smax = %.2f' % Smax, xy=(0.05, 0.95),
#                xycoords='axes fraction', fontsize=14)
#    ax.annotate('Tau = %.2f' % sporulating_capacity,
#                xy=(0.05, 0.85), xycoords='axes fraction', fontsize=14)
    ax.annotate('RMSE : %.4f' %rmse, xy=(0.05, 0.95), 
                xycoords='axes fraction', fontsize=14)
    ax.set_xlabel("Density (nb lesions/cm2)", fontsize = 18)
    ax.set_ylabel("Mean surface of 1 lesion (cm2)", fontsize = 18)

def example_density_robert_2004(**kwds):
    brown_rust = BrownRustFungus()
#    growth_controler = NoPriorityGrowthControl()
    growth_controler = GeometricCircleCompetition()

    def example_single_density(d = 40):
        g, leaf = get_g_and_one_leaf()
        leaf.area = 10.
        leaf.green_area = 10.
#        df_temp = pd.DataFrame([24.], 
#                                columns = ['temp'])
        df_temp = pd.DataFrame([18.], 
                                columns = ['temp'])
        leaf.temperature_sequence = df_temp['temp']
        les = brown_rust.lesion(mutable = False, group_dus = True,**kwds)
        les.set_position([[1] for i in range(int(d*leaf.area))])
        leaf.lesions = [les]
        for i in range(160+21*18):
            les.update(leaf = leaf)
            growth_controler.control(g)
        surfs = les.surface / les.nb_lesions
        surfs_spo = les.surface_spo / les.nb_lesions
        return surfs, surfs_spo

    filename = 'calibration_lesion_rust_2004.csv'
    df_obs = pd.read_csv(filename, sep = ';')
    df_sim = pd.DataFrame(columns = ['surf_mean', 'surf_spo_mean'])
    # x = np.sort(np.concatenate([df_obs['density'], np.arange(1, 41, 2)]))
    x = np.sort(np.concatenate([df_obs['density'], np.arange(1,1)]))
    for d in x:
        df_sim.loc[len(df_sim)+1] = example_single_density(d)
    df_sim['density'] = x
    df_sim['severity'] = df_sim['density'] * df_sim['surf_spo_mean']
    df_obs['severity'] = df_obs['density'] * df_obs['surface']

    rmse = np.sqrt(((df_sim[df_sim['density'].isin(df_obs['density'])]['surf_spo_mean'] -
                        df_obs['surface']) ** 2).mean())

    fig, ax = plt.subplots()
    ax.plot(df_sim.density, df_sim.surf_spo_mean)
    ax.plot(df_obs.density, df_obs.surface, 'ro')
#    ax.annotate('Smax = %.2f' % Smax, xy=(0.05, 0.95),
#                xycoords='axes fraction', fontsize=14)
#    ax.annotate('Tau = %.2f' % sporulating_capacity,
#                xy=(0.05, 0.85), xycoords='axes fraction', fontsize=14)
    ax.annotate('RMSE : %.4f' %rmse, xy=(0.05, 0.95), 
                xycoords='axes fraction', fontsize=14)
    ax.set_xlabel("Density (nb lesions/cm2)", fontsize = 18)
    ax.set_ylabel("Mean surface of 1 lesion (cm2)", fontsize = 18)
    
    fig, ax = plt.subplots()
    ax.plot(df_sim.density, df_sim.severity)
    ax.plot(df_obs.density, df_obs.severity, 'ro')
    ax.annotate('RMSE : %.4f' %rmse, xy=(0.05, 0.95), 
                xycoords='axes fraction', fontsize=14)
    ax.set_xlabel("Density (nb lesions/cm2)", fontsize = 18)
    ax.set_ylabel("Severity in sporulation", fontsize = 18)
    ax.set_ylim([0., 1.])

def example_density_robert_2004_complete(temperature = 14., **kwds):
    def _parse(date, hour):
        return datetime(int(date[6:10]),int(date[3:5]),int(date[0:2]),int(hour))
    weather = pd.read_csv('temperature_frezal.csv', sep = ';', 
                       parse_dates={'datetime':['date', 'hour']},
                       date_parser=_parse)
    weather['day'] = map(lambda x : x.dayofyear -  weather['datetime'][0].dayofyear, weather['datetime'])
    weather.set_index('datetime', inplace = True)
    
    brown_rust = BrownRustFungus()
#    growth_controler = NoPriorityGrowthControl()
    growth_controler = GeometricCircleCompetition()

    filename = 'calibration_lesion_rust_2004_complete.csv'
    df_obs = pd.read_csv(filename, sep = ';')
    df_obs['severity'] = df_obs['density'] * df_obs['area']
    dates_obs = np.unique(df_obs['date'])
    dates = {}
#    dates = {date:date*temperature for date in dates_obs}
    columns = ['date', 'density', 'area', 'severity']
    densities = np.arange(1, 50, 5)
    df_sim = pd.DataFrame([[np.nan for col in columns] 
                           for i in range(len(dates_obs)*len(densities))],
                           columns = columns)
    indx = 0
    for d in densities:
        g, leaf = get_g_and_one_leaf()
        leaf.area = 10.
        leaf.green_area = 10.
        les = brown_rust.lesion(mutable = False, group_dus = True, **kwds)
        les.set_position([[1] for i in range(int(d*leaf.area))])
        leaf.lesions = [les]
        i_dates = 0
        day = 0
        for date in weather.index:
            leaf.temperature_sequence = [weather.loc[date, 'temperature']]
            les.update(leaf = leaf)
            growth_controler.control(g)
            if leaf.lesions[0].age_sporulation > 0:
                if day == 0.:
                    print date
                    print les.age_tt
                    day = weather.loc[date, 'day']
                if (day > 0 and date.hour == 12 and i_dates < len(dates_obs) and
                    weather.loc[date, 'day'] - day + 1 >= dates_obs[i_dates]):
                    surf_spo = les.surface_spo / les.nb_lesions
                    spo = {'date':dates_obs[i_dates], 'density':d, 
                           'area':surf_spo, 'severity':d*surf_spo}
                    df_sim.loc[indx, :] = pd.Series(spo)
                    indx +=1
                    if d == densities[0]:
                        dates[dates_obs[i_dates]] = les.age_sporulation
                    i_dates += 1
                    
#            if i_dates < len(dates):
#                date = dates_obs[i_dates]
#                if leaf.lesions[0].age_sporulation >= dates[date]:
#                    surf_spo = les.surface_spo / les.nb_lesions
#                    spo = {'date':date, 'density':d, 'area':surf_spo}
#                    df_sim.loc[indx, :] = pd.Series(spo)
#                    indx +=1
#                    i_dates += 1
                    
    fig, axs = plt.subplots(2, 4)
    dates_iter = iter(dates_obs)
    rmses = []
    for i, ax in enumerate(axs.flat):
        date = next(dates_iter)
        df_o = df_obs[df_obs['date']==date]
        df_s = df_sim[df_sim['date']==date]

        # Plot
        ax.plot(df_o['density'], df_o['area'], 'ro')
        ax.plot(df_s['density'], df_s['area'], 'b')            
        
        # Get RMSE
        df_o = df_o.sort('density')
        df_s.loc[0, ['density', 'area']] = [0.,0.]
        df_s = df_s.sort('density')
        
        s = interpolate.interp1d(df_s['density'], df_s['area'])
        df_i = df_o.copy()
        df_i.loc[:,'area'] = [s(d) for d in df_i['density']]
        rmses.append(np.sqrt((df_i['area'].astype(float) - df_o['area'].astype(float)) ** 2).mean())
        
        if i<4:
            ax.set_ylim([0, 0.016])
        else:
            ax.set_ylim([0, 0.03])
            ax.set_xlabel('Lesion density', fontsize = 16)
        ax.set_xlim([0, 45])
        if i%4 == 0:
            ax.set_ylabel('Lesion size (in cm2)', fontsize = 16)
        ax.annotate('%d days\n%d Cd' %(date, dates[date]), xy=(0.7, 0.8), 
                    xycoords='axes fraction', fontsize=14)
    
    rmse = np.mean(rmses)
    axs[0][0].annotate('RMSE : %.4f' %rmse, xy=(0.05, 0.9), 
                xycoords='axes fraction', fontsize=14)

    fig, axs = plt.subplots(2, 4)
    dates_iter = iter(dates_obs)
    rmses = []
    for i, ax in enumerate(axs.flat):
        date = next(dates_iter)
        df_o = df_obs[df_obs['date']==date]
        df_s = df_sim[df_sim['date']==date]

        # Plot
        ax.plot(df_o['density'], df_o['severity'], 'ro')
        ax.plot(df_s['density'], df_s['severity'], 'b')            
        
        # Get RMSE
        df_o = df_o.sort('density')
        df_s.loc[0, ['density', 'severity']] = [0.,0.]
        df_s = df_s.sort('density')
        
        s = interpolate.interp1d(df_s['density'], df_s['severity'])
        df_i = df_o.copy()
        df_i.loc[:,'severity'] = [s(d) for d in df_i['density']]
        rmses.append(np.sqrt((df_i['severity'].astype(float) - df_o['severity'].astype(float)) ** 2).mean())
        
        if i>=4:
            ax.set_xlabel('Lesion density', fontsize = 16)
        ax.set_ylim([0, 1.])
        ax.set_xlim([0, 45])
        if i%4 == 0:
            ax.set_ylabel('Severity in sporulation (%)', fontsize = 16)
        ax.annotate('%d days\n%d Cd' %(date, dates[date]), xy=(0.7, 0.8), 
                    xycoords='axes fraction', fontsize=14)
    rmse = np.mean(rmses)
    axs[0][0].annotate('RMSE : %.4f' %rmse, xy=(0.05, 0.9), 
                xycoords='axes fraction', fontsize=14)

def example_density_robert_2005(temperature = 18., **kwds):
    brown_rust = BrownRustFungus()
    growth_controler = NoPriorityGrowthControl()

    filename = 'calibration_lesion_rust_2005.csv'
    df_obs = pd.read_csv(filename, sep = ';')
    natures = np.array(['spo', 'chlo', 'nec'])
    dates = np.unique(df_obs['date'])
    dates_tt = np.array([11, 16, 21, 28, 36])
    dates = {date:dates_tt[i]*temperature for i, date in enumerate(dates)}
    columns = ['nature', 'date', 'density', 'area']
    df_sim = pd.DataFrame([[np.nan for col in columns] 
                            for i in range(3*5*len(np.arange(1, 121, 10)))],
                            columns = columns)
    indx = 0
    densities = np.concatenate([np.arange(1, 10), np.arange(10, 130, 10)])
    for d in densities:
        g, leaf = get_g_and_one_leaf()
        leaf.area = 10.
        leaf.green_area = 10.
        df_temp = pd.DataFrame([temperature], columns = ['temp'])
        leaf.temperature_sequence = df_temp['temp']
        les = brown_rust.lesion(mutable = False, group_dus = True, **kwds)
        les.set_position([[1] for i in range(int(d*leaf.area))])
        leaf.lesions = [les]
        i_dates = 1
        for i in range(1500):
            les.update(leaf = leaf)
            growth_controler.control(g)
            if i_dates <= len(dates):
                if leaf.lesions[0].age_sporulation >= dates[i_dates]:
#                if leaf.lesions[0].age_tt >= dates[i_dates]:
                    surf_chlo = les.surface_chlo/leaf.area
                    surf_spo = les.surface_spo/leaf.area
                    surf_nec = les.surface_empty/leaf.area
                    chlo = {'date':i_dates, 'nature':'chlo', 'density':d,
                           'area':surf_chlo}
                    df_sim.loc[indx, :] =  pd.Series(chlo)
                    indx +=1
                    spo = {'date':i_dates, 'nature':'spo', 'density':d,
                           'area':surf_spo}
                    df_sim.loc[indx, :] = pd.Series(spo)
                    indx +=1
                    nec = {'date':i_dates, 'nature':'nec', 'density':d,
                           'area':surf_nec}
                    df_sim.loc[indx, :] = pd.Series(nec)
                    indx +=1
                    i_dates += 1
#    df_sim.loc[:, 'surface_visible'] = df_sim['chlo'] + df_sim['spo'] + \
#                                        df_sim['nec']
    
    fig, axs = plt.subplots(len(natures), len(dates))
    states = iter(['sporulation', 'chlorosis', 'necrosis'])
    dates_tt = iter(dates_tt)
    dates = dates.itervalues()
    rmses = []
    for i, axs_line in enumerate(axs):
        nature = natures[i]
        for j, ax in enumerate(axs_line):
            date = j+1
            df_o = df_obs[(df_obs['nature']==nature) &
                            (df_obs['date']==date)]
            df_s = df_sim[(df_sim['nature']==nature) &
                            (df_sim['date']==date)]

            # Plot
            ax.plot(df_o['density'], df_o['area'], 'ro')
            ax.plot(df_s['density'], df_s['area'], 'b')            
            
            # Get RMSE
            df_o = df_o.sort('density')
            df_s.loc[0, ['density', 'area']] = [0.,0.]
            df_s = df_s.sort('density')
            
            s = interpolate.interp1d(df_s['density'], df_s['area'])
            df_i = df_o.copy()
            df_i.loc[:,'area'] = [s(d) for d in df_i['density']]
            rmses.append(np.sqrt((df_i['area'].astype(float) - df_o['area'].astype(float)) ** 2).mean())
            

            ax.set_xlim([0,120])
            ax.set_ylim([0,1])
            if j == 0:
                ax.set_ylabel('Severity in '+ next(states), fontsize = 16)
            if i == 0:
                ax.set_title('%d days - %d Cd' %(next(dates_tt), next(dates)))
    rmse = np.mean(rmses)
    axs[0][0].annotate('RMSE : %.4f' %rmse, xy=(0.05, 0.9), 
                xycoords='axes fraction', fontsize=14)

def example_density_complete(**kwds):
    def _parse(date, hour):
        return datetime(int(date[6:10]),int(date[3:5]),int(date[0:2]),int(hour))

    def adapt_date_2005(i_date):
        return np.array([11, 16, 21, 28, 36])[i_date-1]
                
    weather = pd.read_csv('temperature_frezal.csv', sep = ';', 
                       parse_dates={'datetime':['date', 'hour']},
                       date_parser=_parse)
    weather['day'] = map(lambda x : x.dayofyear -  weather['datetime'][0].dayofyear, weather['datetime'])
    weather.set_index('datetime', inplace = True)
    
    brown_rust = BrownRustFungus()
#    growth_controler = NoPriorityGrowthControl()
    growth_controler = GeometricCircleCompetition()

    filename = 'calibration_lesion_rust_2004_complete.csv'
    df_obs_2004 = pd.read_csv(filename, sep = ';')
    df_obs_2004['severity'] = df_obs_2004['density'] * df_obs_2004['area']
    filename = 'calibration_lesion_rust_2005.csv'
    df_obs_2005 = pd.read_csv(filename, sep = ';')
    df_obs_2005 = df_obs_2005[df_obs_2005['nature']=='spo']
    df_obs_2005 = df_obs_2005.rename(columns = {'area':'severity'})
    df_obs_2005['area'] = df_obs_2005['severity'] / df_obs_2005['density']
    df_obs_2005['area'][np.isnan(df_obs_2005['area'])] = 0.
#    df_obs_2005 = df_obs_2005[df_obs_2005['density']<=45]
    
    dates_obs = np.unique(df_obs_2004['date'])
    dates = {}
#    dates = {date:date*temperature for date in dates_obs}
    columns = ['date', 'density', 'area', 'severity']
#    densities = np.arange(1, 50, 5)
    densities = np.arange(1, 112, 10)
    df_sim = pd.DataFrame([[np.nan for col in columns] 
                           for i in range(len(dates_obs)*len(densities))],
                           columns = columns)
    indx = 0
    for d in densities:
        g, leaf = get_g_and_one_leaf()
        leaf.area = 10.
        leaf.green_area = 10.
        les = brown_rust.lesion(mutable = False, group_dus = True, **kwds)
        les.set_position([[1] for i in range(int(d*leaf.area))])
        leaf.lesions = [les]
        i_dates = 0
        day = 0
        for date in weather.index:
            leaf.temperature_sequence = [weather.loc[date, 'temperature']]
            les.update(leaf = leaf)
            growth_controler.control(g)
            if leaf.lesions[0].age_sporulation > 0:
                if day == 0.:
                    print date
                    print les.age_tt
                    day = weather.loc[date, 'day']
                if (day > 0 and date.hour == 12 and i_dates < len(dates_obs) and
                    weather.loc[date, 'day'] - day + 1 >= dates_obs[i_dates]):
                    surf_spo = les.surface_spo / les.nb_lesions
                    spo = {'date':dates_obs[i_dates], 'density':d, 
                           'area':surf_spo, 'severity':d*surf_spo}
                    df_sim.loc[indx, :] = pd.Series(spo)
                    indx +=1
                    if d == densities[0]:
                        dates[dates_obs[i_dates]] = les.age_sporulation
                    i_dates += 1
                    
    fig, axs = plt.subplots(2, 4)
    dates_iter = iter(dates_obs)
    rmses = []
    for i, ax in enumerate(axs.flat):
        date = next(dates_iter)
        (df_o_2004, df_o_2005) = map(lambda df: df[df['date']==date],
                                    (df_obs_2004, df_obs_2005))
        df_s = df_sim[df_sim['date']==date]

        # Plot
        ax.plot(df_o_2004['density'], df_o_2004['area'], 'ro')
        if len(df_o_2005)>0:
            ax.plot(df_o_2005['density'], df_o_2005['area'], 'bo')
        ax.plot(df_s['density'], df_s['area'], 'k')
        
        # Get RMSE
        df_s.loc[0, ['density', 'area']] = [0.,0.]
        df_s = df_s.sort('density')
        
        s = interpolate.interp1d(df_s['density'], df_s['area'])
        df_o = pd.concat([df_o_2004, df_o_2005])
        df_o.sort('density')
        df_i = df_o.copy()
        df_i.loc[:,'area'] = [s(d) for d in df_i['density']]
        rmses.append(np.sqrt((df_i['area'].astype(float) - df_o['area'].astype(float)) ** 2).mean())
        
        if i<4:
            ax.set_ylim([0, 0.016])
        else:
            ax.set_ylim([0, 0.03])
            ax.set_xlabel('Lesion density', fontsize = 16)
        ax.set_xlim([0, max(densities)])
        if i%4 == 0:
            ax.set_ylabel('Lesion size (in cm2)', fontsize = 16)
        ax.annotate('%d days\n%d Cd' %(date, dates[date]), xy=(0.7, 0.8), 
                    xycoords='axes fraction', fontsize=14)
    
    rmse = np.mean(rmses)
    axs[0][0].annotate('RMSE : %.4f' %rmse, xy=(0.05, 0.9), 
                xycoords='axes fraction', fontsize=14)

    fig, axs = plt.subplots(2, 4)
    dates_iter = iter(dates_obs)
    rmses = []
    for i, ax in enumerate(axs.flat):
        date = next(dates_iter)
        (df_o_2004, df_o_2005) = map(lambda df: df[df['date']==date],
                                    (df_obs_2004, df_obs_2005))
        df_s = df_sim[df_sim['date']==date]

        # Plot
        ax.plot(df_o_2004['density'], df_o_2004['severity'], 'ro')
        if len(df_o_2005)>0:
            ax.plot(df_o_2005['density'], df_o_2005['severity'], 'bo')
        ax.plot(df_s['density'], df_s['severity'], 'k')            
        
        # Get RMSE
        df_s.loc[0, ['density', 'severity']] = [0.,0.]
        df_s = df_s.sort('density')
        
        s = interpolate.interp1d(df_s['density'], df_s['severity'])
        df_o = pd.concat([df_o_2004, df_o_2005])
        df_o.sort('density')
        df_i = df_o.copy()
        df_i.loc[:,'severity'] = [s(d) for d in df_i['density']]
        rmses.append(np.sqrt((df_i['severity'].astype(float) - df_o['severity'].astype(float)) ** 2).mean())
        
        if i>=4:
            ax.set_xlabel('Lesion density', fontsize = 16)
        ax.set_ylim([0, 1.])
        ax.set_xlim([0, max(densities)])
        if i%4 == 0:
            ax.set_ylabel('Severity in sporulation (%)', fontsize = 16)
        ax.annotate('%d days\n%d Cd' %(date, dates[date]), xy=(0.7, 0.8), 
                    xycoords='axes fraction', fontsize=14)
    rmse = np.mean(rmses)
    axs[0][0].annotate('RMSE : %.4f' %rmse, xy=(0.05, 0.9), 
                xycoords='axes fraction', fontsize=14)


def example_production_spo(**kwds):
    weather = read_weather_year(2012)
    df = weather.data
    df = df.reset_index()
    df = df[df['degree_days']>0]
    leaf = get_one_leaf()
    brown_rust = BrownRustFungus()
    lesion = brown_rust.lesion(mutable = False, group_dus = False, **kwds)
    stock = []
    daily_prod = []
    tt = []
    cum_tt = 0.
    count = 0.
    nb_spores = 0.
    for i, row in df.iterrows():
        if count < 5000:
            count += 1
            leaf.temperature_sequence = [20.]
#            leaf.temperature_sequence = [row['temperature_air']]
            lesion.update(leaf = leaf)
            lesion.control_growth(growth_offer = lesion.growth_demand)
            new_nb_spores = lesion.stock_spores - nb_spores
            stock.append(lesion.stock_spores)
#            daily_prod.append(new_nb_spores / lesion.surface_spo if lesion.surface_spo >0. else 0.)
            daily_prod.append(new_nb_spores)
            nb_spores = lesion.stock_spores

            cum_tt += lesion.delta_thermal_time_growth(leaf_temperature = [row['temperature_air']])
            tt.append(cum_tt)

    # daily_prod = np.concatenate([[0.], np.diff(stock)])
    # daily_prod = np.diff(stock)

    fig, ax = plt.subplots()
    ax.plot(tt, daily_prod)
#    ax.plot(tt, stock)
    ax.set_xlabel("Thermal time (Teff)", fontsize = 18)
    ax.set_ylabel("Daily spore production", fontsize = 18)

def example_dispersal(age_canopy = 1400., nb_dispersal_units = 1e5,
                      variety = 'Mercia', nplants = 200, 
                      density_factor = 1, vmax = 25,
                      k_dispersal = 0.16):
    reconst = EchapReconstructions()
    stand_density_factor = {'Mercia':1., 'Rht3':1, 'Tremie12':0.8, 'Tremie13':0.8}
    stand_density_factor[variety] = density_factor
    adel = reconst.get_reconstruction(name=variety, nplants = nplants, nsect = 7,
                                      stand_density_factor = stand_density_factor)
    g = adel.setup_canopy(age = age_canopy)
    fungus = BrownRustFungus()
    dispersor = BrownRustDispersal(fungus = fungus,
                                   domain = adel.domain,
                                   domain_area = adel.domain_area,
                                   k_dispersal = k_dispersal)
    dispersor.plot_distri_layers(g, nb_dispersal_units)
    dispersor.view_distri_layers(g, nb_dispersal_units, vmax = vmax)
    
def example_inoculation(age_canopy = 1400., density_inoculum = 2000.,
                        variety = 'Mercia', nplants = 200, 
                        density_factor = 1, vmax = None):    
    reconst = EchapReconstructions()
    stand_density_factor = {'Mercia':1., 'Rht3':1, 'Tremie12':0.8, 'Tremie13':0.8}
    stand_density_factor[variety] = density_factor
    adel = reconst.get_reconstruction(name=variety, nplants = nplants, nsect = 7,
                                      stand_density_factor = stand_density_factor)
    g = adel.setup_canopy(age = age_canopy)
    contaminator = AirborneContamination(fungus = BrownRustFungus(),
                                         domain_area = adel.domain_area)
    df = contaminator.plot_distri_layers(g, density_inoculum)
    contaminator.view_distri_layers(g, density_inoculum, vmax = vmax)
    return sum(df[1])/adel.domain_area

def disperse(g,
             emission_model=None,
             transport_model=None,
             fungus_name='', 
             label="LeafElement",
             weather_data=None):
    if weather_data is None:
        DU = emission_model.get_dispersal_units(g, fungus_name, label)
    else: 
        DU = emission_model.get_dispersal_units(g, fungus_name, label, weather_data)
    labels = g.property('label')
   
    # Transport of dispersal units
    if sum(DU.values())>0:
        # (C. Fournier, 22/11/2013 : keep compatibility with previous protocol and add a new one (used in septo3d/echap)
        if weather_data is not None:
            deposits = transport_model.disperse(g, DU, weather_data)
        else:
            deposits = transport_model.disperse(g, DU) # update DU in g , change position, status 
            
        # Allocation of new dispersal units
        for vid, dlist in deposits.iteritems():
            if len(dlist)>0 and labels[vid].startswith(label):
                leaf = g.node(vid)
                try:
                    leaf.dispersal_units += dlist
                except:
                    leaf.dispersal_units = dlist
    return g

def setup_simu(sowing_date="2000-10-15 12:00:00", 
               end_date="2001-05-25 01:00:00", start_date = None,
               variety = 'Mercia', nplants = 30, nsect = 7, age_canopy = 0.,
               TT_delay = 20, dispersal_delay = 24,
               k_dispersal = 1.85, **kwds):
    # Get weather
    weather = get_weather(start_date=sowing_date, end_date=end_date)
    
    # Set canopy
    reconst = EchapReconstructions()
    adel = reconst.get_reconstruction(name=variety, nplants=nplants, nsect=nsect)
    g = adel.setup_canopy(age=age_canopy)
    
    # Manage temporal sequence  
    if start_date is None:
        start_date = sowing_date
    seq = pd.date_range(start=start_date, end=end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase=0.)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay=TT_delay)
    every_dispersal = time_filter(seq, delay=dispersal_delay)
    rust_filter = thermal_time_filter(seq, weather, TTmodel, delay=TT_delay)
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')
    dispersal_timing = IterWithDelays(*time_control(seq, every_dispersal, weather.data))
    rust_timing = CustomIterWithDelays(*time_control(seq, rust_filter, weather.data), eval_time='end')
    
    # Set up models
    fungus = BrownRustFungus()
    fungus.parameters(**kwds)
    recorder = BrownRustRecorder()
    growth_controler = NoPriorityGrowthControl()
    infection_controler = BiotrophDUProbaModel()
    dispersor = BrownRustDispersal(fungus = fungus,
                                   domain = adel.domain,
                                   domain_area = adel.domain_area,
                                   k_dispersal = k_dispersal)
    return (g, adel, fungus,  canopy_timing, dispersal_timing, rust_timing, 
            recorder, growth_controler, infection_controler, dispersor)

def example_frezal(k_dispersal = 0.16, source = 2, position_source = 4, 
                   start_date = "2000-05-04 12:00:00",
                   end_date = "2000-05-29 12:00:00",
                   **kwds):
    # Cf table 3 frezal 2009 / revoir dates
    (g, adel, fungus, canopy_timing, dispersal_timing, rust_timing, recorder, 
     growth_controler, infection_controler, 
     dispersor) = setup_simu(sowing_date=str(int(start_date[:4])-1)+"-10-15 12:00:00", 
                   end_date = end_date,
                   start_date = start_date,
                   variety = 'Mercia', nplants = 30, nsect = 7, 
                   age_canopy = 1400., TT_delay = 20, dispersal_delay = 24, 
                   k_dispersal = k_dispersal, **kwds)
    
    # Inoculate all leaves 2 with 1 lesion
    labels = g.property('label')
    a_labels = adel_labels(g, scale = 5)
    plants = [p for p,l in labels.iteritems() if l.startswith('plant')]
    sources = []    
    for i_pl, pl in enumerate(plants):
        a_labs = {k:v for k,v in a_labels.iteritems()
                    if v.startswith('plant'+str(i_pl+1)+'_MS')
                    and 'LeafElement' in v}
        fnl = int(a_labs[max(a_labs.keys())].split('_')[2].split('metamer')[1])
        f_top = fnl-source+1
        sources += [k for k,l in a_labs.iteritems() 
                    if int(l.split('_')[2].split('metamer')[1])==f_top and
                    int(l.split('_')[4].split('LeafElement')[1])==position_source]
    for source in sources:
        les = fungus.lesion()
        les.set_position([0.5, 0.5])
        g.node(source).lesions = [les]
    
    # Simulation loop
    for i, controls in enumerate(zip(canopy_timing, 
                                     dispersal_timing, 
                                     rust_timing)):
        canopy_iter, dispersal_iter, rust_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            g = adel.grow(g, canopy_iter.value)
        
        # Get weather for date and add it as properties on leaves
        if rust_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = rust_iter.value.temperature_air,
                           wetness_sequence = rust_iter.value.wetness)

        # Develop disease (infect for dispersal units and update for lesions)
        if rust_iter:
            infect(g, rust_iter.dt, infection_controler, label='LeafElement')
            update(g, rust_iter.dt, growth_controler, label='LeafElement')
            
        # Disperse and wash
        geom = g.property('geometry')
        if dispersal_iter and len(geom)>0:
            disperse(g, dispersor, dispersor,
                     fungus_name = "brown_rust",
                     label='LeafElement', 
                     weather_data=dispersal_iter.value)
        
        # Save outputs
        if rust_iter:
            date = rust_iter.value.index[-1]
            print date
            recorder.record(g, date, 
                            degree_days = rust_iter.value.degree_days[-1])
    
    recorder.post_treatment(variety = 'Mercia')
    ax = plot_by_leaf(recorder.data, variable = 'nb_lesions', 
                      xaxis = 'degree_days', return_ax = True)
    ax.annotate('position on source : %d' %position_source, xy=(0.05, 0.9), 
                xycoords='axes fraction', fontsize=14)
    ax.set_ylabel("Number of lesions on leaf", fontsize = 18)
    ax.set_xlabel("Thermal time (Cd)", fontsize = 18)
    return g, recorder
    
def example_frezal_complete(k_dispersal = 0.16, position_source = 5, **kwds):
    leaves = [1, 2, 3]    
    starts = {2000 : {lf:"2000-05-04 12:00:00" for lf in leaves},
              2001 : {1: "2001-05-07 12:00:00",
                      2: "2001-04-27 12:00:00",
                      3: "2001-04-23 12:00:00"}}
    ends = {2000 : {lf:"2000-06-05 12:00:00" for lf in leaves},
            2001 : {1: "2001-06-11 12:00:00",
                    2: "2001-06-01 12:00:00",
                    3: "2001-05-31 12:00:00"}}
    df_out = pd.DataFrame(columns = ['year', 'source', 'target', 
                                 'percentage', 'multiplication_rate'])
    indx = 0
    for source in leaves:
        for year in [2000, 2001]:
            g, recorder = example_frezal(k_dispersal = k_dispersal,
                                         source = source,
                                         position_source = position_source,
                                         start_date=starts[year][source],
                                         end_date=ends[year][source])
            df = recorder.data
            deposits = {}
            for target in leaves:
                df_lf = df[(df['axis'].isin(['MS'])) & (df['num_leaf_top']==target)]
                df_lf = df_lf[pd.notnull(df_lf.loc[:,'nb_lesions'])].loc[:, ['degree_days', 'nb_lesions']]
                df_mean = df_lf.groupby('degree_days').mean()
                deposits[target] = df_mean['nb_lesions'].iloc[-1]
            
            total_deposits = sum(deposits.values())
            for target in leaves:
                indx +=1
                line = {'year':year, 'source':source, 'target':target, 
                        'percentage':deposits[target]/total_deposits,
                        'multiplication_rate':total_deposits}
                df_out.loc[indx, :] = pd.Series(line)
    
    df_percentage = df_out.set_index(['source', 'year', 'target'])
    df_percentage = df_percentage.unstack('target')['percentage']*100
    df_rate = df_out.loc[:, ['source', 'year', 'multiplication_rate']]
    df_rate = df_rate.drop_duplicates().set_index(['year', 'source'])
    return df_percentage, df_rate
    
def read_multiplication_frezal():
    filename = 'multiplication_frezal.csv'
    df = pd.read_csv(filename, sep = ';')
    return df.groupby(['year','leaf']).mean()
    
def read_percentage_frezal():
    df_percentage = pd.read_csv('percentage_frezal.csv', sep = ';')
    df_percentage = df_percentage.set_index(['source', 'year', 'target'])
    return df_percentage.unstack('target')['percentage']
    
def plot_sim_obs_frezal(df_sim, df_obs):
    import itertools 
    from math import floor
    fig, axs = plt.subplots(3,2)    
    df_sim = df_sim.reset_index()
    df_obs = df_obs.reset_index()
    colors = itertools.cycle(['b', 'r', 'g'])
    for i, ax in enumerate(axs.flat):
        x = [1, 2]
        src = int(floor(i/2)) + 1
        yr = 2000 + i%2
        (sim, obs) = map(lambda x: x[(x['source']==src) & (x['year']==yr)], 
                         (df_sim, df_obs))
        y_bottom = np.array([0., 0.])
        for target in np.arange(3., 0., -1.):
            y = [obs.loc[:, int(target)].values[0], sim.loc[:, str(target)].values[0]]            
            ax.bar(x, y, bottom=y_bottom, color = next(colors), align = 'center')
            y_bottom += np.array(y)
            ax.set_xlim([0, 3])
            ax.set_ylim([0, 100])
            ax.set_xticks(x)
            ax.set_xticklabels(['Obs', 'Sim'])
        ax.set_title('Source F%d - Year %d' %(src, yr))
        if i%2 == 0:
            ax.set_ylabel('Distribution\n by leaf (%)', fontsize = 14)
        if i == 1:
            proxy = [plt.Rectangle((0,0), 0,0, facecolor='g'),
                     plt.Rectangle((0,0), 0,0, facecolor='r'),
                     plt.Rectangle((0,0), 0,0, facecolor='b')]
            labels = ['F1', 'F2', 'F3']
            ax.legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    
def external_contamination(g, 
             contamination_source, 
             contamination_model,
             weather_data = None, 
             label = 'LeafElement', **kwds):
    """ Innoculates fungal objects (dispersal units) on elements of the MTG
        according to contamination_model.
    """
    stock = contamination_source.emission(g, weather_data, **kwds)
    labels = g.property('label')
    if stock > 0:
        # Allocation of stock of inoculum
        deposits = contamination_model.contaminate(g, stock, weather_data, label = label)
        stock = 0 # stock has been used (avoid uncontrolled future re-use)

        # Allocation of new dispersal units
        for vid, dlist in deposits.iteritems():
            if len(dlist)>0 and labels[vid].startswith(label):
                leaf = g.node(vid)
                try:
                    leaf.dispersal_units += dlist
                except:
                    leaf.dispersal_units = dlist
    return g
    
def example_annual_loop(variety = 'Mercia', nplants = 30, 
                        year = 2012, sowing_date = '10-15', 
                        density_dispersal_units = 1000):
    (g, adel, fungus, canopy_timing, dispersal_timing, rust_timing, recorder, 
     growth_controler, infection_controler, 
     dispersor) = setup_simu(sowing_date=str(year-1)+"-"+sowing_date+" 12:00:00", 
                   end_date=str(year)+"-07-01 00:00:00",
                   variety = variety, nplants = nplants, nsect = 7, 
                   TT_delay = 20, dispersal_delay = 24)
    
    contaminator = AirborneContamination(fungus = BrownRustFungus(),
                                         domain_area = adel.domain_area)
                                         
    # Simulation loop
    for i, controls in enumerate(zip(canopy_timing, 
                                     dispersal_timing, 
                                     rust_timing)):
        canopy_iter, dispersal_iter, rust_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            g = adel.grow(g, canopy_iter.value)
        
        # Get weather for date and add it as properties on leaves
        if rust_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = rust_iter.value.temperature_air,
                           wetness_sequence = rust_iter.value.wetness)

        # Simulate airborne contamination
        geom = g.property('geometry')
        if dispersal_iter and len(geom)>0:
            external_contamination(g, contaminator, contaminator, 
                                   density_dispersal_units = density_dispersal_units)
            
        # Develop disease (infect for dispersal units and update for lesions)
        if rust_iter:
            infect(g, rust_iter.dt, infection_controler, label='LeafElement')
            update(g, rust_iter.dt, growth_controler, label='LeafElement')
            
        # Disperse and wash
        if dispersal_iter and len(geom)>0:
            disperse(g, dispersor, dispersor,
                     fungus_name = "brown_rust",
                     label='LeafElement', 
                     weather_data=dispersal_iter.value)
        
        # Save outputs
        if rust_iter:
            date = rust_iter.value.index[-1]
            print date
            recorder.record(g, date, 
                            degree_days = rust_iter.value.degree_days[-1])
    
    recorder.post_treatment(variety = 'Mercia')
    ax = plot_by_leaf(recorder.data, variable = 'ratio_spo', 
                      xaxis = 'degree_days', return_ax = True)
    ax.set_title(variety+' '+str(year)+' sowing on '+sowing_date, fontsize = 18)
    ax.set_ylabel("Severity in sporulation (percent)", fontsize = 18)
    ax.set_xlabel("Thermal time (Cd)", fontsize = 18)
    return g, recorder

def plot_weather_annual_loop(year = 2012, sowing_date = '10-15', 
                             xlims = [0, 2500], title = None):
    sowing_date = str(year-1)+"-"+sowing_date+" 12:00:00"
    end_date = str(year)+"-07-01 00:00:00"
    weather = get_weather(start_date=sowing_date, end_date=end_date)
    plot_wetness_and_temp(weather, xlims = xlims, title = title)

def run_sensitivity_analysis():
    pass

def example_competition_complex_no_priority():    
    g, leaf = get_g_and_one_leaf()
    leaf.area = 10.
    leaf.green_area = 10.
    brown_rust = BrownRustFungus()
    rust = [brown_rust.lesion(mutable = False, group_dus = True) for i in range(150)]
    septoria = SeptoriaFungus()    
    septo = [septoria.lesion(mutable = False, group_dus = True) for i in range(150)]
    for l in septo:
        l.set_position([1.,1.])
    tot_les = rust + septo
    np.random.shuffle(tot_les)
    leaf.lesions = tot_les
    growth_controler = NoPriorityGrowthControl()
    df = pd.DataFrame(columns=['degree_days', 'sev_rust',
                               'sev_septo', 'sev_tot'])
    cum_tt = 0.
    for i in range(1000):
        leaf.temperature_sequence = [24.]
        cum_tt += 1.
        df.loc[i, 'degree_days'] = cum_tt
        update(g, 1., growth_controler, label='LeafElement')

        lesions = leaf.lesions
        r_surf = sum([l.surface for l in lesions if isinstance(l, brown_rust.lesion().__class__)])
        s_surf = sum([l.surface for l in lesions if isinstance(l, septoria.lesion().__class__)])
        df.loc[i, 'sev_rust'] = r_surf/leaf.area
        df.loc[i, 'sev_septo'] = s_surf/leaf.area
    df.loc[:, 'sev_tot'] = df['sev_rust'] + df['sev_septo']
    df = df.set_index('degree_days')
    df.plot()    

def example_competition_complex(nb_rust = 100, nb_septo = 100):
    g, leaf = get_g_and_one_leaf()
    leaf.area = 10.
    leaf.green_area = 10.
    brown_rust = BrownRustFungus()
    rust = [brown_rust.lesion(mutable = False, group_dus = True) for i in range(nb_rust)]
    septoria = SeptoriaFungus()    
    septo = [septoria.lesion(mutable = False, group_dus = True) for i in range(nb_septo)]
    for l in septo:
        l.set_position([1.,1.])
    tot_les = rust + septo
    np.random.shuffle(tot_les)
    leaf.lesions = tot_les
    growth_controler = SeptoRustCompetition(SeptoModel = septoria.lesion().__class__,
                                            RustModel = brown_rust.lesion().__class__)
    df = pd.DataFrame(columns=['degree_days', 
                               'sev_rust', 'sev_septo', 'sev_tot', 
                               'surf_rust', 'surf_septo', 'surf_tot',
                               'nb_les_rust', 'nb_les_septo', 'nb_les_tot'])
    cum_tt = 0.
    for i in range(500):
        leaf.temperature_sequence = [24.]
        cum_tt += 1.
        df.loc[i, 'degree_days'] = cum_tt
        update(g, 1., growth_controler, label='LeafElement')

        lesions = leaf.lesions
        r_surf = sum([l.surface for l in lesions if isinstance(l, brown_rust.lesion().__class__)])
        s_surf = sum([l.surface for l in lesions if isinstance(l, septoria.lesion().__class__)])
        df.loc[i, 'sev_rust'] = r_surf/leaf.area
        df.loc[i, 'sev_septo'] = s_surf/leaf.area
        df.loc[i, 'surf_rust'] = r_surf
        df.loc[i, 'surf_septo'] = s_surf
        df.loc[i, 'nb_les_rust'] = len([l for l in lesions if isinstance(l, brown_rust.lesion().__class__)])
        df.loc[i, 'nb_les_septo'] = len([l for l in lesions if isinstance(l, septoria.lesion().__class__)])
    df.loc[:, 'sev_tot'] = df['sev_rust'] + df['sev_septo']
    df.loc[:, 'surf_tot'] = df['surf_rust'] + df['surf_septo']
    df.loc[:, 'nb_les_tot'] = df['nb_les_rust'] + df['nb_les_septo']

    fig, ax = plt.subplots()
    ax.plot(df['degree_days'], df['sev_rust'])
    ax.plot(df['degree_days'], df['sev_septo'])
    ax.plot(df['degree_days'], df['sev_tot'])
    ax.legend(['sev_rust', 'sev_septo', 'sev_tot'], loc = 'best')
    
    fig, ax = plt.subplots()
    ax.plot(df['degree_days'], df['surf_rust'])
    ax.plot(df['degree_days'], df['surf_septo'])
    ax.plot(df['degree_days'], df['surf_tot'])
    ax.legend(['surf_rust', 'surf_septo', 'surf_tot'], loc = 'best')

    fig, ax = plt.subplots()
    ax.plot(df['degree_days'], df['nb_les_rust'])
    ax.plot(df['degree_days'], df['nb_les_septo'])
    ax.plot(df['degree_days'], df['nb_les_tot'])
    ax.legend(['nb_les_rust', 'nb_les_septo', 'nb_les_tot'], loc = 'best')