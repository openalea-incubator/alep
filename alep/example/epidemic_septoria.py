""" Simulation #1 of the paper. Epidemics of septoria on wheat. """

# Imports ##########################################################
import random as rd
import numpy as np
import pandas

from alinea.alep.wheat import adel_mtg2
from alinea.alep.architecture import *
from alinea.alep.disease_operation import *
from alinea.alep.disease_outputs import *
from alinea.alep.inoculation import InoculationFirstLeaves
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.alep.protocol import *
from alinea.alep.alep_color import alep_colormap, green_yellow_red

from alinea.alep.alep_weather import add_wetness
from alinea.weather.global_weather import *
from datetime import datetime

from openalea.vpltk import plugin

from alinea.astk.TimeControl import *

# Simulation of an epidemics of septoria ###########################
def simulate_epi_septo(data_file='meteo01.csv',
                       start_date=None, end_date=None):
    """ Simulate an epidemics of septoria under weather given in data_file.
    
    Parameters
    ----------
    data_file: string
        Name of a file with weather data whose format complies with
        the recommendations in 'alinea.weather'
    start_date: format datetime
        Date of start of simulation. Must belong to data_file.
        Ex : start_date = datetime(2000, 10, 1, 1, 00, 00)
    end_date: format datetime
        Date of end of simulation. Must belong to data_file.
        Ex : start_date = datetime(2000, 12, 31, 1, 00, 00)
        
    Returns
    -------
    dates: numpy array
        List of dates of the simulation.
    avg_severity: float
        Average severity on the whole plant.
    ..TODO:: Define outputs
    """
    # Read weather data file and add wetness in weather data
    weather = Weather(data_file = 'meteo01.csv')
    weather = add_wetness(weather)
    
    # Find dates of start and end of simulation if not given as inputs
    if not start_date:
        start_date = weather.data.datetime[0]
    if not end_date:
        end_date = weather.data.datetime[-1]
        
    # Set schedule of calls for each model
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))

    wheat_timing = TimeControl(delay=1, steps = nb_steps)
    septo_timing = TimeControl(delay=1, steps = nb_steps)
    weather_timing = TimeControl(delay=1, steps = nb_steps)
    timer = TimeControler(wheat = wheat_timing, 
                          disease = septo_timing,
                          weather = weather_timing)
    
    # Generate a wheat canopy
    g = adel_mtg2()
    set_properties(g,label = 'LeafElement', surface=5.,
                    healthy_surface=5., position_senescence=None)
    
    # Create a pool of dispersal units for initiation and distribute it
    diseases=plugin.discover('alep.disease')
    septoria = diseases['septoria_exchanging_rings'].load()
    nb_du = 5
    dispersal_units = generate_stock_du(nb_du, disease=septoria)
    inoculator = InoculationFirstLeaves()
    initiate(g, dispersal_units, inoculator)
    
    # Call models that will be used in disease interface
    controler = NoPriorityGrowthControl()
    emitter = SeptoriaRainEmission()
    transporter = Septo3DSplash(reference_surface=1./200)
    
    # Initiate output list
    avg_severity = np.zeros(nb_steps)
    
    # Simulation loop
    for t in timer:
        # Get weather for date t
        if timer.numiter==1:
            date = start_date
        mgc, globalclimate = weather.get_weather(t['weather'].dt, date)
        
        # Set weather properties on every leaf element
        set_properties(g,label = 'LeafElement',
                        wetness=globalclimate.wetness.values[0],
                        temp=globalclimate.temperature_air.values[0],
                        rain_intensity=globalclimate.rain.values[0],
                        rain_duration=1.,
                        relative_humidity=globalclimate.relative_humidity.values[0],
                        wind_speed=globalclimate.wind_speed.values[0])
    
        # Compute fungal processes
        infect(g, t['disease'].dt, label='LeafElement')
        update(g, t['disease'].dt, controler, label='LeafElement')
        if globalclimate.rain.values[0]>0:
            disperse(g, emitter, transporter, "septoria", label='LeafElement')
        
        # Save output
        avg_severity[timer.numiter-1] = compute_total_severity(g)
        
        # /!\ TEMP /!\
        # Refill pool of initial inoculum to simulate differed availability of inoculum
        if timer.numiter%10 == 0 and timer.numiter < 100:
            initiate(g, dispersal_units, inoculator)
            
        # Advance date for next simulation step
        date = weather.next_date(t['weather'].dt, date)
    
    # Save dates of simulation as an output
    dates = weather.data.datetime[weather.data.datetime>=start_date]
    dates = dates[dates<=end_date]
    
    return dates, avg_severity