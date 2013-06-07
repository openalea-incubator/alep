""" Weather utilities for alep.

Gather methods concerning weather management specific to alep functioning.
"""
from pylab import *
import pandas
from alinea.weather.mini_models import leaf_wetness_rapilly

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
    
def add_dispersal_events(weather, dispersal_event_model):
    """ Add a column to indicate the dispersal events,
        Add a column to indicate the hour with max rain and value of max rain,
        Add a column to indicate the duration of rain at same line.
    """
    max_rain = 0.
    ind_max = 0.
    dispersal = []
    rain_counter = 0.
    rain = zeros(len(weather.data))
    rain_duration = zeros(len(weather.data))
    for i_line in range(len(weather.data)):
        dispersal.append(dispersal_event_model(weather.data.rain[i_line], 
                                               weather.data.relative_humidity[i_line]))
        if dispersal[i_line] == True:
            rain_counter += 1.
            if weather.data.rain[i_line] > max_rain:
                max_rain = weather.data.rain[i_line]
                ind_max = i_line
        elif ind_max > 0.:
            rain[ind_max] = max_rain
            rain_duration[ind_max] = rain_counter
            max_rain = 0.
            ind_max = 0.
            rain_counter = 0.
    
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