""" Make images of a simulation of vine/oidium epidemics for a movie.
"""
# Useful imports
import random as rd
import numpy as np
import pandas


# Imports for weather
import pandas
import alinea.septo3d
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather
from alinea.alep.alep_weather import wetness_rapilly
from alinea.astk.TimeControl import *

# Imports for grapevine
from alinea.alep.vine import Vine
from alinea.astk.plant_interface import new_canopy, grow_canopy
from alinea.alep.architecture import set_properties, update_healthy_area, add_area_topvine

# Imports for powdery mildew
from openalea.vpltk import plugin
from alinea.alep.protocol import *
from alinea.alep.powdery_mildew import *
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.dispersal_emission import PowderyMildewWindEmission
from alinea.alep.dispersal_transport import PowderyMildewWindDispersal
from alinea.alep.growth_control import GrowthControlVineLeaf
from alinea.alep.infection_control import BiotrophDUProbaModel

# Imports for display and saving
from alinea.alep.disease_outputs import plot_severity_vine, save_image

# Initialization ##############################################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Read weather and adapt it to powdery mildew (add wetness)
meteo_path = shared_data(alinea.septo3d, 'meteo05-06.txt')
weather = Weather(data_file=meteo_path)
weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
seq = pandas.date_range(start = "2006-03-01 01:00:00", end = "2006-04-10 01:00:00", freq='H')

# Initialize a vine canopy
vine = Vine()
g,_ = new_canopy(vine, age = 6)

# Initialize the models for powdery mildew
diseases=plugin.discover('alep.disease')
powdery_mildew = diseases['powdery_mildew'].load()
infection_controler = BiotrophDUProbaModel()
growth_controler = GrowthControlVineLeaf()
emitter = PowderyMildewWindEmission()
transporter = PowderyMildewWindDispersal()

# Define the schedule of calls for each model
every_h = time_filter(seq, delay=1)
every_24h = time_filter(seq, delay=24)
weather_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
vine_timing = IterWithDelays(*time_control(seq, every_24h, weather.data))
mildew_timing = IterWithDelays(*time_control(seq, every_h, weather.data))

# Simulation #########################################################
for ind, controls in enumerate(zip(weather_timing, vine_timing, mildew_timing)):
    weather_eval, vine_eval, mildew_eval = controls
    
    # Get weather for date and add it as properties on leaves
    if weather_eval:
        set_properties(g,label = 'lf',
                       temp = weather_eval.value.temperature_air[0],
                       rain_intensity = weather_eval.value.rain.mean(),
                       wetness = weather_eval.value.wetness[0],
                       relative_humidity = weather_eval.value.relative_humidity[0],
                       wind_speed = weather_eval.value.wind_speed[0],
                       wind_direction = (1,0,0))

    # Grow vine canopy
    if vine_eval:
        print(vine_eval.value.datetime[0])
        g,_ = grow_canopy(g, vine, vine_eval.value)
        add_area_topvine(g)
    
    # Develop disease
    if mildew_eval:
        # Update g for the disease
        update_healthy_area(g, label = 'lf')
        
        # Refill pool of initial inoculum to simulate differed availability
        labels = g.property('label')
        if ind < 500 and labels[188]=='lf' and ind%50==0:
            dispersal_units = generate_stock_du(nb_dus=100, disease=powdery_mildew)
            try:
                g.node(188).dispersal_units += dispersal_units
            except: 
                g.node(188).dispersal_units = dispersal_units
               
        # Update dispersal units and lesions
        infect(g, mildew_eval.dt, infection_controler, label='lf')
        update(g, mildew_eval.dt, growth_controler, senescence_model=None, label='lf')
        
        # Disperse
        disperse(g, emitter, transporter, "powdery_mildew", label='lf')

    if vine_eval:
        scene = plot_severity_vine(g, trunk=True, transparency=0.9)
        index = ind/24
        if index < 10 :
            image_name='./images_powdery_mildew/image0000%d.png' % index
        elif index < 100 :
            image_name='./images_powdery_mildew/image000%d.png' % index
        elif index < 1000 :
            image_name='./images_powdery_mildew/image00%d.png' % index
        elif index < 10000 :
            image_name='./images_powdery_mildew/image0%d.png' % index
        else :
            image_name='./images_powdery_mildew/image%d.png' % index
        save_image(scene, image_name=image_name)