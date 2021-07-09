""" Example generating images for a movie of a septoria epidemic.

Careful: 
    - This simulation is long
    - Set the angle of the viewer before saving images
"""
    
# General useful imports
import random as rd
import numpy as np

# Imports for weather
import pandas
import alinea.septo3d
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather
from alinea.alep.alep_weather import wetness_rapilly
from alinea.astk.TimeControl import *

# Imports for wheat
from alinea.alep.wheat import initialize_stand
from alinea.astk.plant_interface import grow_canopy
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septo3d_v2 import plugin_septoria
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.alep.dispersal_transport import SeptoriaRainDispersal
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence

# Imports for display and saving
from alinea.alep.disease_outputs import plot_severity_by_leaf, save_image

# Initiation of the simulation ##########################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Read weather and adapt it to septoria (add wetness)
meteo_path = shared_data(alinea.septo3d, 'meteo05-06.txt')
weather = Weather(data_file=meteo_path)
weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
seq = pandas.date_range(start = "2005-10-01 01:00:00", end = "2006-07-01 01:00:00", freq='H')
                        
# Initialize a wheat canopy
g, wheat, domain_area, domain = initialize_stand(age=0., length=1, width=1,
    sowing_density=250, plant_density=250, inter_row=0.12,  seed=3)

# Choose source leaf in canopy 
# (Here the value of the leaf is known but it changes with another initialize_stand)
source_leaf = g.node(21943)

# Initialize the models for septoria
septoria = plugin_septoria()
inoculator = RandomInoculation()
growth_controler = NoPriorityGrowthControl()
infection_controler = BiotrophDUProbaModel()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
emitter = SeptoriaRainEmission(domain_area=domain_area)
transporter = SeptoriaRainDispersal()
washor = RapillyWashing()

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
        set_properties(g,label = 'LeafElement',
                       temp = weather_eval.value.temperature_air[0],
                       wetness = weather_eval.value.wetness[0],
                       relative_humidity = weather_eval.value.relative_humidity[0],
                       wind_speed = weather_eval.value.wind_speed[0])
    if rain_eval:
        set_properties(g,label = 'LeafElement',
                       rain_intensity = rain_eval.value.rain.mean(),
                       rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)
    
    # Grow wheat canopy
    if wheat_eval:
        g,_ = grow_canopy(g, wheat, wheat_eval.value)
        # Note : The position of senescence goes back to its initial value after
        # a while for undetermined reason
        # --> temporary hack for keeping senescence position low when it is over
        positions = g.property('position_senescence')
        greens = g.property('is_green')
        areas = g.property('area')
        senesced_areas = g.property('senesced_area')
        leaves = get_leaves(g, label = 'LeafElement')
        vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
        positions.update({vid:(0 if (positions[vid]==1 and not greens[vid]) or
                                    (positions[vid]>0 and round(areas[vid],5)==round(senesced_areas[vid],5))
                                    else positions[vid]) for vid in vids})
                                    
    # Develop disease
    if septo_eval:
        # Update g for the disease
        sen_model.find_senescent_lesions(g, label = 'LeafElement')
        update_healthy_area(g, label = 'LeafElement')
        
        # Possibly refill pool of initial inoculum to simulate differed availability
        if rain_eval and i <= 700 and source_leaf.geometry!=None:
            dus = generate_stock_du(nb_dus=rd.randint(0,5), disease=septoria)
            try:
                source_leaf.dispersal_units += dus
            except:
                source_leaf.dispersal_units = dus
        
        # Update dispersal units and lesions
        infect(g, septo_eval.dt, infection_controler, label='LeafElement')
        update(g, septo_eval.dt, growth_controler, sen_model, label='LeafElement')
        
    # Disperse and wash
    if rain_eval:
        if rain_eval.value.rain.mean()>0:
            g, nb = disperse(g, emitter, transporter, "septoria", label='LeafElement')
            wash(g, washor, rain_eval.value.rain.mean(), label='LeafElement')
            
    # Display outputs and save image
    if wheat_eval:
        scene = plot_severity_by_leaf(g, senescence=False, transparency=0.9)
        index = i/24
        if index < 10 :
            image_name='./images_septo/image0000%d.png' % index
        elif index < 100 :
            image_name='./images_septo/image000%d.png' % index
        elif index < 1000 :
            image_name='./images_septo/image00%d.png' % index
        elif index < 10000 :
            image_name='./images_septo/image0%d.png' % index
        else :
            image_name='./images_septo/image%d.png' % index
        save_image(scene, image_name=image_name)
            