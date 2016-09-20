""" Follow disease progress on a single leaf with various inoculum """


from alinea.alep.disease_outputs import plot_severity_by_leaf, count_dispersal_units, count_lesions

# General useful imports
import random as rd
import numpy as np
import sys

# Imports for weather
import pandas
import alinea.septo3d
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather
from alinea.alep.alep_weather import wetness_rapilly
from alinea.astk.TimeControl import *

# Imports for wheat
from alinea.alep.wheat import adel_one_leaf_element as leaf
from alinea.alep.architecture import set_properties, update_healthy_area

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septo3d_v2 import plugin_septoria
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.alep.disease_outputs import LeafInspector

def run_sim(year=1999, nb_dus=1, model='septoria_exchanging_rings', **kwds):
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Read weather and adapt it to septoria (add wetness)
    weather_file = 'meteo'+ str(year-1)[-2:] + '-' + str(year)[-2:] + '.txt'
    meteo_path = shared_data(alinea.septo3d, weather_file)
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    seq = pandas.date_range(start = str(year)+"-02-01 01:00:00",
                            end = str(year)+"-05-01 01:00:00", 
                            freq='H')

    # Generation of 1 leaf
    g = leaf()
    set_properties(g, label = 'LeafElement', area=20., green_area=20., position_senescence=1)
    source_leaf = g.node(10)
    blade_id = 8
    inspector = LeafInspector(g, blade_id=blade_id)

    # Initialize the models for septoria
    if 'alinea.alep.'+model in sys.modules:
        del(sys.modules['alinea.alep.'+model])
    septoria = plugin_septoria(model)
    DU = septoria.dispersal_unit(**kwds)
    growth_controler = NoPriorityGrowthControl()

    # Define the schedule of calls for each model
    every_h = time_filter(seq, delay=1)
    every_24h = time_filter(seq, delay=24)
    every_rain = rain_filter(seq, weather)
    weather_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    wheat_timing = IterWithDelays(*time_control(seq, every_24h, weather.data))
    septo_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))

    # Simulation ############################################################
    for i, controls in enumerate(zip(weather_timing, wheat_timing, septo_timing, rain_timing)):
        weather_eval, wheat_eval, septo_eval, rain_eval = controls
        
        # Get weather for date and add it as properties on leaves
        if weather_eval:
            # set_properties(g,label = 'LeafElement',
                           # temp = weather_eval.value.temperature_air[0],
                           # wetness = weather_eval.value.wetness[0],
                           # relative_humidity = weather_eval.value.relative_humidity[0],
                           # wind_speed = weather_eval.value.wind_speed[0])
            set_properties(g,label = 'LeafElement',
                           temp = 22.,
                           wetness = True,
                           relative_humidity = 90,
                           wind_speed = 0.)
        if rain_eval:
            # set_properties(g,label = 'LeafElement',
                           # rain_intensity = rain_eval.value.rain.mean(),
                           # rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)
            set_properties(g,label = 'LeafElement',
                           rain_intensity = 0.,
                           rain_duration = 0.) 
             
        # Develop disease
        if septo_eval:
            # Update g for the disease
            update_healthy_area(g, label = 'LeafElement')
            
            # Inoculation
            if i==1:
                # if i<=2.5*24:
                    # nb_dus=i
                # else:
                    # nb_dus= i - (i/2.5*24)
                dus = [DU(nb_spores=1., status='deposited') for i in range(int(nb_dus))]
                try:
                    source_leaf.dispersal_units += dus
                except:
                    source_leaf.dispersal_units = dus
            
            # Update dispersal units and lesions
            infect(g, septo_eval.dt, label='LeafElement')
            
            # if i>=753:
                # import pdb
                # pdb.set_trace()
            update(g, septo_eval.dt, growth_controler, label='LeafElement')
        
            # Save outputs   
            inspector.update_variables(g)
            inspector.update_du_variables(g)
    return inspector

class NoTransport:
    def disperse(self, g, dispersal_units, time_control = None):
        return {}
    
def run_sim_disp(year=1999, nb_dus=1, **kwds):
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Read weather and adapt it to septoria (add wetness)
    weather_file = 'meteo'+ str(year-1)[-2:] + '-' + str(year)[-2:] + '.txt'
    meteo_path = shared_data(alinea.septo3d, weather_file)
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    seq = pandas.date_range(start = str(year)+"-02-01 01:00:00",
                            end = str(year)+"-05-01 01:00:00", 
                            freq='H')

    # Generation of 1 leaf
    g = leaf()
    set_properties(g, label = 'LeafElement', area=20., green_area=20., position_senescence=1)
    source_leaf = g.node(10)
    blade_id = 8
    inspector = LeafInspector(g, blade_id=blade_id)

    # Initialize the models for septoria
    if 'alinea.alep.septoria_exchanging_rings' in sys.modules:
        del(sys.modules['alinea.alep.septoria_exchanging_rings'])
    septoria = plugin_septoria()
    DU = septoria.dispersal_unit(**kwds)
    growth_controler = NoPriorityGrowthControl()
    emitter = SeptoriaRainEmission(domain_area=1./2500)
    transporter = NoTransport()

    # Define the schedule of calls for each model
    every_h = time_filter(seq, delay=1)
    every_24h = time_filter(seq, delay=24)
    every_rain = rain_filter(seq, weather)
    weather_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    wheat_timing = IterWithDelays(*time_control(seq, every_24h, weather.data))
    septo_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))

    # Simulation ############################################################
    for i, controls in enumerate(zip(weather_timing, wheat_timing, septo_timing, rain_timing)):
        weather_eval, wheat_eval, septo_eval, rain_eval = controls
        
        # Get weather for date and add it as properties on leaves
        if weather_eval:
            # set_properties(g,label = 'LeafElement',
                           # temp = weather_eval.value.temperature_air[0],
                           # wetness = weather_eval.value.wetness[0],
                           # relative_humidity = weather_eval.value.relative_humidity[0],
                           # wind_speed = weather_eval.value.wind_speed[0])
            set_properties(g,label = 'LeafElement',
                           temp = 18.,
                           wetness = True,
                           relative_humidity = 90,
                           wind_speed = 0.,
                           rain_intensity = 0.,
                           rain_duration = 0.)
        if i==1000:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = 1.,
                           rain_duration = 1.) 
        # if rain_eval:
            # set_properties(g,label = 'LeafElement',
                           # rain_intensity = rain_eval.value.rain.mean(),
                           # rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)

             
        # Develop disease
        if septo_eval:
            # Update g for the disease
            update_healthy_area(g, label = 'LeafElement')
            
            # Inoculation
            if i==1:
                # if i<=2.5*24:
                    # nb_dus=i
                # else:
                    # nb_dus= i - (i/2.5*24)
                dus = [DU(nb_spores=1., status='deposited') for i in range(int(nb_dus))]
                try:
                    source_leaf.dispersal_units += dus
                except:
                    source_leaf.dispersal_units = dus
            
            # Update dispersal units and lesions
            infect(g, septo_eval.dt, label='LeafElement')
            update(g, septo_eval.dt, growth_controler, label='LeafElement')
            
        if i==1000:
            g, nb = disperse(g, emitter, transporter, "septoria", label='LeafElement')
        
        # Save outputs   
        inspector.update_variables(g)
        inspector.update_du_variables(g)
    return inspector  