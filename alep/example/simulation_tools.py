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
from alinea.echap.architectural_reconstructions import echap_reconstructions
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
    
def wheat_path((year, variety, nplants, nsect)):
    if variety.lower().startswith('tremie'):
        variety = 'tremie'
    return './adel/'+variety.lower()+'_'+str(int(year))+'_'+str(nplants)+'pl_'+str(nsect)+'sect'

def init_canopy(adel, wheat_dir, rain_and_light=False):
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
                wheat_dir, wheat_is_loaded=True, rain_and_light=False):
    if wheat_is_loaded:
        newg, TT = adel.load(it_wheat, dir=wheat_dir)
        move_properties(g, newg)
        return newg
    else:
        g = adel.grow(g, canopy_iter.value)
        if rain_and_light==True:
            rain_and_light_star(g, light_sectors = '1', 
                                domain=adel.domain, convUnit=adel.convUnit)
        return g
        
def make_canopy(start_date = "2010-10-29 12:00:00", 
                end_date = "2011-07-01 01:00:00",
                variety = 'Tremie13', nplants = 30, nsect = 7, disc_level = 5, 
                wheat_dir = './adel/tremie_2013_30pl_7sect', 
                reset_reconst = True):
    """ Simulate and save canopy (prior to simulation). """        
    # Manage weather and scheduling
    weather = get_weather(start_date = start_date, end_date = end_date)
    seq = pandas.date_range(start = start_date, end = end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase = 0)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 20.)
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')

    # Simulate and save wheat development    
    reconst = echap_reconstructions(reset=True, reset_data=True)
    adel = reconst.get_reconstruction(name=variety, nplants=nplants, nsect=nsect)    
    domain = adel.domain
    convUnit = adel.convUnit
    g = adel.setup_canopy(age=0.)
    rain_and_light_star(g, light_sectors = '1', domain = domain, convUnit = convUnit)
    it_wheat = 0
    adel.save(g, it_wheat, dir=dir)
    for i, canopy_iter in enumerate(canopy_timing):
        if canopy_iter:
            it_wheat += 1
            g = adel.grow(g, canopy_iter.value)
            rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            adel.save(g, it_wheat, dir=dir)