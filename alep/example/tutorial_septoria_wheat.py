""" Tutorial for septoria"""
import random as rd
import numpy as np
import pandas

from alinea.adel.astk_interface import AdelWheat
from alinea.astk.plant_interface import *
#from alinea.alep.wheat import adel_mtg2
from alinea.alep.architecture import *
from alinea.alep.disease_operation import *
from alinea.alep.disease_outputs import *
from alinea.alep.inoculation import InoculationFirstLeaves, RandomInoculation
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.protocol import *
from alinea.alep.alep_color import alep_colormap, green_yellow_red

from alinea.alep.alep_weather import add_wetness
from alinea.weather.global_weather import *
from datetime import datetime

from openalea.vpltk import plugin

from alinea.astk.TimeControl import *

# Useful functions ########################################################################
def update_plot(g):
    # Count lesion surfaces by id & add it as MTG property 
    surface_lesions_by_leaf = count_lesion_surfaces_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'surface_lesions', surface_lesions_by_leaf, label = 'LeafElement')

    # Compute severity by leaf
    severity_by_leaf = compute_severity_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'severity', severity_by_leaf, label = 'LeafElement')

    # Visualization
    g = alep_colormap(g, 'surface_lesions', cmap=green_yellow_red(levels=10), lognorm=False)
    # TODO : Normalize the colormap between 0 and 1 /!\
    scene = plot3d(g)
    Viewer.display(scene)
    return scene
    
# Initiation ##############################################################################
# Define a plant or canopy
wheat = AdelWheat()
g,_ = new_canopy(wheat,age=100)

# Add missing properties needed for the simulation
# The simulation requires the following properties on leaf elements:
#   - 'surface': total surface of the leaf element
#   - 'healthy_surface': surface of the leaf element without lesion or senescence
#   - 'position_senescence': position of the senescence on blade axis
#set_properties(g,label = 'LeafElement',
#               surface=5., healthy_surface=5., position_senescence=None)
               
# discover the disease implementation for septoriose
diseases=plugin.discover('alep.disease')
#septoria_classes = [kls for kls in diseases if kls.load().name == 'septoria']
# or 
try:
    septoria = diseases['septoria_exchanging_rings'].load()
except KeyError:
    from alinea.alep.septoria_exchanging_rings import Disease
    septoria=Disease()
# Create a pool of dispersal units (DU)
nb_du = 5
dispersal_units = generate_stock_du(nb_du, disease=septoria)

# Distribute the DU 
#inoculator = InoculationFirstLeaves()
inoculator = RandomInoculation()
initiate(g, dispersal_units, inoculator)
# g = color.colormap(g,'length',lognorm=False)
# scene = plot3d(g)
# Viewer.display(scene)

# Preparation of the simulation loop #####################################################
# Call models that will be used in disease interface
controler = NoPriorityGrowthControl()
dispersor = Septo3DSplash(reference_surface=1./200)

# Choose dates of simulation
start_date = datetime(2000, 10, 1, 1, 00, 00)
# end_date = datetime(2000, 10, 10, 00, 00)
end_date = datetime(2000, 12, 31, 00, 00)
date = start_date

# Read weather between date and add wetness
weather = Weather(data_file = 'meteo01.csv')
weather = add_wetness(weather)

# Set schedule of calls for each model
nsteps = len(pandas.date_range(start_date, end_date, freq='H'))

wheat_timing = TimeControl(delay=24, steps = nsteps, model = wheat, weather = weather, start_date = start_date)
septo_timing = TimeControl(delay=1, steps = nsteps)
weather_timing = TimeControl(delay=1, steps = nsteps)
plot_timing = TimeControl(delay=24, steps = nsteps)
timer = TimeControler(wheat = wheat_timing, disease = septo_timing,
                      weather = weather_timing, plotting = plot_timing)

for t in timer:
    print(timer.numiter)
    
    # Not needed if wheat model is complete
    #set_properties_on_new_leaves(g,label = 'LeafElement',
    #                         surface=5., healthy_surface=5.,
    #                         position_senescence=None)
    
    # Get weather for date
    mgc, globalclimate = weather.get_weather(t['weather'].dt, date)
    
    set_properties(g,label = 'LeafElement',
                    wetness=globalclimate.wetness.values[0],
                    temp=globalclimate.temperature_air.values[0],
                    rain_intensity=globalclimate.rain.values[0],
                    rain_duration=1.,
                    relative_humidity=globalclimate.relative_humidity.values[0],
                    wind_speed=globalclimate.wind_speed.values[0])
    
    #wheat
    grow_canopy(g,wheat,t['wheat'])
    # Disease
    infect(g, t['disease'].dt, label='LeafElement')
    update(g, t['disease'].dt, controler, label='LeafElement')

    if globalclimate.rain.values[0]>0:
        disperse(g, dispersor, "septoria", label='LeafElement')  

    if t['plotting'].dt > 0:
        print('plotting...')           
        update_plot(g)

    # Refill pool of initial inoculum to simulate differed availability of inoculum
    if timer.numiter%10 == 0 and timer.numiter < 100:
        initiate(g, dispersal_units, inoculator)

    # Advance date for next simulation step
        date = weather.next_date(t['weather'].dt, date)