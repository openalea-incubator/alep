""" Weather utilities for alep.

Gather methods concerning weather management specific to alep functioning.
"""
from pylab import *
import pandas
from alinea.weather.global_weather import *
from alinea.weather.mini_models import leaf_wetness_rapilly

def get_septoria_weather(data_file='meteo01.csv'):
    """ Read weather data and adapt it to septoria model (add wetness)
    
    Parameters
    ----------
    data_file: str
        Name of weather data file (in csv)
    
    Returns
    -------
    weather: pandas dataframe
        Dataframe of weather indexed by date and with explicitely named columns
        See `alinea.weather.global_weather`
    """
    weather = Weather(data_file=data_file)
    weather = add_wetness(weather)
    weather = add_rain_dispersal_events(weather)
    return weather

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