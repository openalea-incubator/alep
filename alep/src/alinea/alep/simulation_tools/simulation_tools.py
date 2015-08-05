# -*- coding: iso-8859-15  -*- 
"""
Main tools to run a simulation of disease epidemics on wheat
"""
import os

# Imports for weather
from alinea.alep.alep_weather import wetness_rapilly, linear_degree_days
from alinea.alep.alep_time_control import *
from alinea.astk.TimeControl import *
from alinea.echap.weather_data import *
import alinea.septo3d
import alinea.alep
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather

# Imports for wheat
from alinea.echap.architectural_reconstructions import (EchapReconstructions, pdict,
                                                        reconstruction_parameters)
from alinea.adel.newmtg import move_properties
from alinea.caribu.caribu_star import rain_and_light_star

def get_weather(start_date="2010-10-15 12:00:00", end_date="2011-06-20 01:00:00"):
    """ Get weather data for simulation. """
    start = datetime.datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    if start.year >= 2010:
        filename = 'Boigneville_0109'+str(start.year)+'_3108'+str(start.year+1)+'_h.csv'
        meteo_path = shared_data(alinea.echap, filename)
        weather = Weather(meteo_path, reader = arvalis_reader)
        weather.check(['temperature_air', 'PPFD', 'relative_humidity',
                       'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
        notation_dates_file = shared_data(alinea.alep, 'notation_dates/notation_dates_'+str(start.year+1)+'.csv')
        weather.check(varnames=['notation_dates'], models={'notation_dates':add_notation_dates}, notation_dates_file = notation_dates_file)
    else:
        start_yr = str(start.year)[2:4]
        end_yr = str(start.year+1)[2:4]
        filename = 'meteo'+ start_yr + '-' + end_yr + '.txt'
        meteo_path = shared_data(alinea.septo3d, filename)
        weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    weather.check(varnames=['degree_days'], models={'degree_days':linear_degree_days}, start_date=start_date, base_temp=0., max_temp=30.)
    # weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':linear_degree_days}, start_date=start_date, base_temp=0., max_temp=25.)
    weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':linear_degree_days}, start_date=start_date, base_temp=0., max_temp=30.)
    return weather

def count_available_canopies(year, variety, nplants, nsect):
    filepath = str(shared_data(alinea.alep)/'wheat_reconstructions')
    list_files = os.listdir(filepath)
    start = variety.lower()+'_'+str(int(year))+'_'+\
            str(nplants)+'pl_'+str(nsect)+'sect'
    return len([f for f in list_files if f.startswith(start)])

def wheat_path(year, variety, nplants, nsect, rep):
    if rep is None:
        rep = ''
    filepath = str(shared_data(alinea.alep)/'wheat_reconstructions')
    return filepath+'/'+variety.lower()+'_'+str(int(year))+'_'+\
            str(nplants)+'pl_'+str(nsect)+'sect_rep'+str(rep)

def init_canopy(adel, wheat_dir, rain_and_light=True):
    import pdb
    pdb.set_trace()
    if os.path.exists(wheat_dir):
        wheat_is_loaded = True
        it_wheat = 0
        g, TT = adel.load(it_wheat, dir=wheat_dir)
    else:
        wheat_is_loaded = False
        g = adel.setup_canopy(age=0.)
        if rain_and_light==True:
            rain_and_light_star(g, light_sectors = '1', 
                                domain=adel.domain, convUnit=adel.convUnit)
    return g, wheat_is_loaded

def grow_canopy(g, adel, canopy_iter, it_wheat,
                wheat_dir, wheat_is_loaded=True, rain_and_light=True):
    if wheat_is_loaded:
        newg, TT = adel.load(it_wheat, dir=wheat_dir)
        import pdb
        pdb.set_trace()
        move_properties(g, newg)
        return newg
    else:
        g = adel.grow(g, canopy_iter.value)
        if rain_and_light==True:
            rain_and_light_star(g, light_sectors = '1', 
                                domain=adel.domain, convUnit=adel.convUnit)
        return g
    
def alep_echap_reconstructions():
    pars = reconstruction_parameters()
    pars['density_tuning'] = pdict(None)
    pars['density_tuning']['Tremie12'] = 0.85
    pars['density_tuning']['Tremie13'] = 0.85
    reconst = EchapReconstructions(reset_data=True, pars=pars)
    reconst.axepop_fits['Tremie12'].MS_probabilities = {12:0.21, 13:0.79}
    reconst.axepop_fits['Tremie13'].MS_probabilities = {11:23./43, 12:20./43.}
    return reconst
    
def make_canopy(year = 2013, variety = 'Tremie13', sowing_date = '10-29',
                nplants = 15, nsect = 7, nreps=10, fixed_rep=None):
    """ Simulate and save canopy (prior to simulation). """    
    def make_canopy_one_rep(rep=0):
        reconst = alep_echap_reconstructions()
        adel = reconst.get_reconstruction(name=variety, nplants=nplants, nsect=nsect)    
        domain = adel.domain
        convUnit = adel.convUnit
        g = adel.setup_canopy(age=0.)
        rain_and_light_star(g, light_sectors = '1', domain = domain, convUnit = convUnit)
        it_wheat = 0
        wheat_dir = wheat_path(year, variety, nplants, nsect, rep)
        adel.save(g, it_wheat, dir=wheat_dir)
        for i, canopy_iter in enumerate(canopy_timing):
            if canopy_iter:
                it_wheat += 1
                g = adel.grow(g, canopy_iter.value)
                rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
                adel.save(g, it_wheat, dir=wheat_dir)
                
    # Manage weather and scheduling
    start_date=str(year-1)+"-"+sowing_date+" 12:00:00"
    end_date=str(year)+"-07-01 00:00:00"
    weather = get_weather(start_date=start_date, end_date=end_date)
    seq = pandas.date_range(start=start_date, end=end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase = 0)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 20.)
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')

    # Simulate and save wheat development
    if fixed_rep is None:
        for rep in range(nreps):
            make_canopy_one_rep(rep=rep)
    else:
        make_canopy_one_rep(rep=fixed_rep)