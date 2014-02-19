""" Weather utilities for alep.

Gather methods concerning weather management specific to alep functioning.
"""
from pylab import *
import pandas
from alinea.weather.global_weather import *
from alinea.weather.mini_models import leaf_wetness_rapilly

def is_raining(rain_eval):
    """ Check if it is raining or not
    
    Parameters
    ----------
    rain_eval: TimeControl EvalValue instance
        Weather data divided according to rain occurences
    
    Returns
    -------
    True or False
    """
    if rain_eval:
        if rain_eval.value.rain.sum() > 0:
            return True
        else:
            return False
    else:
        return False

def get_septoria_weather(data_file='meteo00-01.txt', sep = ';'):
    """ Read weather data and adapt it to septoria model (add wetness)
    
    Parameters
    ----------
    data_file: str
        Name of weather data file (in txt)
    
    Returns
    -------
    weather: pandas dataframe
        Dataframe of weather indexed by date and with explicitely named columns
        See `astk.Weather`
    """
    meteo_path = shared_data(alinea.septo3d, 'meteo00-01.txt')
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    
    # weather = Weather(data_file=data_file, sep = sep)
    # weather = add_wetness(weather)
    # weather = add_rain_dispersal_events(weather)
    return weather

def wetness_rapilly(data):
    """ Compute leaf wetness as in Rapilly et Jolivet, 1976 as a 
        function of rain or relative humidity and PAR.
    
    Parameters
    ----------
    data: pandas dataframe 
        Weather data as read by Weather class
    
    Returns
    -------
    sequence of True or False
        True if the leaf is wet, False oherwise
    """
    return np.greater(data[['rain']].values, 0.) + np.less(data[['PPFD']].values,644.) * np.greater_equal(data[['relative_humidity']].values, 85.)

def mean_rain(data):
    """ Create a sequence of rain events with only mean values."""
    pass
    
def add_wetness(weather):
    """ Complete weather data with wetness.
    """
    wet = dict(wetness=[])
    for i_line in range(len(weather.data)):
        wet['wetness'].append(leaf_wetness_rapilly(weather.data.rain[i_line], 
                              weather.data.relative_humidity[i_line],
                              weather.data.PPFD[i_line]))
    wetness = pandas.DataFrame(wet)
    wetness.index = weather.data.rain.index
    weather.data = weather.data.join(wetness)
    return weather

def add_septoria_infection_risk(data):
    """ Add True or False if there is a risk of infection for septoria.
    
    Infection is possible if wetness duration > 10h and 10<temp<30 deg C. 
    """
    temp1 = data.temperature_air>=10
    temp2 = data.temperature_air<=30
    infect_cond = data.wetness * temp1 * temp2
    septo_infection_risk = np.zeros(len(data))
    counter = 0.
    for i_line in range(len(data)):
        if infect_cond[i_line]==True:
            counter += 1.
            if counter >= 10.:
                septo_infection_risk[i_line] = True
            else:
                septo_infection_risk[i_line] = False
        else:
            counter = 0.
            septo_infection_risk[i_line] = False
    return septo_infection_risk
    
    
def add_rain_dispersal_events(weather):
    """ Add a column to indicate the dispersal events at hours with max rain for each dispersal event,
        Add a column to indicate the value of average rain during dispersal event,
        Add a column to indicate the duration of each dispersal event.
    """
    max_rain = 0.
    ind_max = 0.
    rain_counter = 0.
    rain_amount = 0.
    
    dispersal = zeros(len(weather.data))
    rain = zeros(len(weather.data))
    rain_duration = zeros(len(weather.data))
    for i_line in range(len(weather.data)):
        if weather.data.rain[i_line] > 0.:
            rain_counter += 1.
            rain_amount += weather.data.rain[i_line]
            if weather.data.rain[i_line] > max_rain:
                max_rain = weather.data.rain[i_line]
                ind_max = i_line
        else:
            if ind_max > 0.:
                rain[ind_max] = rain_amount/rain_counter
                rain_duration[ind_max] = rain_counter
                dispersal[ind_max] = 1
                max_rain = 0.
                ind_max = 0.
                rain_counter = 0.
                rain_amount = 0.
    
    del weather.data['rain']
    rain = pandas.DataFrame(dict(rain=rain))
    rain.index = weather.data.temperature_air.index
    weather.data = weather.data.join(rain)
    
    rain_duration = pandas.DataFrame(dict(rain_duration=rain_duration))
    rain_duration.index = weather.data.temperature_air.index
    weather.data = weather.data.join(rain_duration)
    
    dispersal_event = pandas.DataFrame(dict(dispersal_event=dispersal))
    dispersal_event.index = weather.data.temperature_air.index
    weather.data = weather.data.join(dispersal_event)
    return weather