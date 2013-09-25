""" Tutorial for septoria"""
import random as rd
import numpy as np
import pandas
from datetime import datetime

from alinea.astk.plantgl_utils import *

from alinea.astk.plant_interface import grow_canopy
from alinea.alep.wheat import initialize_stand

from alinea.alep.architecture import (set_properties, update_healthy_area,
                                      set_property_on_each_id, get_leaves,
                                      get_leaves)
from alinea.alep.septoria import plugin_septoria
from alinea.alep.disease_operation import *
from alinea.alep.disease_outputs import (compute_lesion_areas_by_leaf,
                                        compute_severity_by_leaf)
from alinea.alep.inoculation import InoculationFirstLeaves, RandomInoculation
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.protocol import *
from alinea.alep.alep_color import alep_colormap, green_yellow_red

from alinea.alep.alep_weather import get_septoria_weather
from alinea.weather.global_weather import *

from alinea.astk.TimeControl import *

from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer

# Useful functions ############################################################
def update_plot(g):
    # Count lesion surfaces by id & add it as MTG property
    surface_lesions_by_leaf = compute_lesion_areas_by_leaf(g, label = 'LeafElement')
    # print('Lesion Surface ', surface_lesions_by_leaf)
    # print('')get_septoria_weather
    # vids = get_leaves(g)
    # print('Position Senescence ', {vid:g.node(vid).position_senescence for vid in vids})
    # print('')
    # print('Healhty Surface ', g.property('healthy_area'))
    set_property_on_each_id(g, 'surface_lesions', surface_lesions_by_leaf, label = 'LeafElement')

    # Compute severity by leaf
    severity_by_leaf = compute_severity_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'severity', severity_by_leaf, label = 'LeafElement')

    # Visualization
    g = alep_colormap(g, 'severity', cmap=green_yellow_red(levels=100),
                      lognorm=False, zero_to_one=False, vmax=100)

    leaves = get_leaves(g, label='LeafElement')
    pos_sen = g.property('position_senescence')
    for leaf in leaves:
        if pos_sen[leaf]==0.:
            g.node(leaf).color = (157, 72, 7)

    scene = plot3d(g)
    Viewer.display(scene)
    return scene

# Initialization ##############################################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Choose dates of simulation and initialize the value of date
start_date = datetime(2000, 10, 1, 1, 00, 00)
# end_date = datetime(2001, 03, 01, 00, 00)
end_date = datetime(2000, 12, 31, 00, 00)
date = start_date

# Read weather and adapt it to septoria (add wetness)
# weather = get_septoria_weather(data_file='./../../example/meteo01.csv')
weather = get_septoria_weather(data_file='meteo01.csv')

# Initialize a wheat canopy
g, wheat, domain_area = initialize_stand(age=0., length=0.1,
                                        width=0.2, sowing_density=150,
                                        plant_density=150, inter_row=0.12)

# Note G. GARIN: 04/09/2013 : To move
domain_area_cm2 = domain_area * 10000 # hack (to be changed in adel)
geometries = g.property('geometry')
lai = get_lai(geometries, domain_area_cm2)

# Initialize the models for septoria
septoria = plugin_septoria()
inoculator = RandomInoculation()
controler = NoPriorityGrowthControl()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
dispersor = Septo3DSplash(reference_surface=domain_area)

# Define the schedule of calls for each model
nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
weather_timing = TimeControl(delay=1, steps=nb_steps)
wheat_timing = TimeControl(delay=24, steps=nb_steps, model=wheat, weather=weather, start_date=start_date)
septo_timing = TimeControl(delay=1, steps=nb_steps)
plot_timing = TimeControl(delay=24, steps=nb_steps)
timer = TimeControler(weather=weather_timing, wheat=wheat_timing, disease = septo_timing, plotting=plot_timing)

# Simulation ##################################################################

AA = []
GA = []
HA = []
LA = []

leaves = get_leaves(g, label = 'LeafElement')
green_leaves = {leaf:[] for leaf in leaves}
healthy_leaves = {leaf:[] for leaf in leaves}
sen_leaves = {leaf:[] for leaf in leaves}
len_leaves = {leaf:[] for leaf in leaves}
pos_sen = {leaf:[] for leaf in leaves}

h_l = {leaf:None for leaf in leaves}

for t in timer:
    print(timer.numiter)

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
    update_healthy_area(g, label = 'LeafElement')

    # Note : The position of senescence goes back to its initial value after
    # a while for undetermined reason
    # --> temporary hack for keeping senescence position low when it is over
    positions = g.property('position_senescence')
    are_green = g.property('is_green')
    vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
    positions.update({vid:(0 if positions[vid]==1 and not are_green[vid] else positions[vid])
                      for vid in vids})

    # Disease
    infect(g, t['disease'].dt, label='LeafElement')
    update(g, t['disease'].dt, controler, sen_model, label='LeafElement')

    if globalclimate.rain.values[0]>0:
        disperse(g, dispersor, "septoria", label='LeafElement')

    if t['plotting'].dt > 0:
        # print('plotting...')
        update_plot(g)

    # Refill pool of initial inoculum to simulate differed availability of inoculum
    if timer.numiter%10 == 0 and timer.numiter <= 500:
        nb_dus = 10
        dus = generate_stock_du(nb_dus, disease=septoria)
        initiate(g, dus, inoculator)

    # if timer.numiter==3193:
        # import pdb
        # pdb.set_trace()
        # To see what's happening on leaf 10...

    # Advance date for next simulation step
    date = weather.next_date(t['weather'].dt, date)

    AA.append(sum(aa for aa in g.property('area').itervalues()))
    GA.append(sum(ga for ga in g.property('green_area').itervalues()))
    HA.append(sum(ha for ha in g.property('healthy_area').itervalues()))
    les_area = compute_lesion_areas_by_leaf(g, label = 'LeafElement')
    LA.append(sum(la for la in les_area.itervalues()))
    for leaf in leaves:
        sen_leaves[leaf].append(g.node(leaf).position_senescence)
        green_leaves[leaf].append(g.node(leaf).green_area)
        healthy_leaves[leaf].append(g.node(leaf).healthy_area)
        len_leaves[leaf].append(g.node(leaf).length)
        pos_sen[leaf].append(g.node(leaf).position_senescence)
        
    # Temp
    # for leaf in leaves:
        # if g.node(leaf).healthy_area>0. and h_l[leaf]==None:
            # h_l[leaf] = 0.
        # if g.node(leaf).healthy_area==0. and h_l[leaf]==0.:
            # h_l[leaf] = 1.
        # if g.node(leaf).healthy_area>0. and h_l[leaf]==1.:
            # import pdb
            # pdb.set_trace()
    # print(h_l)
    
    # if False and t['plotting'].dt > 0:
        # print('area', g.property('area'))
        # print('')
        # print('green area', g.property('green_area'))
        # print('')
        # print('healthy area', g.property('healthy_area'))
        # print('')
        # print('Lesion area', les_area)
        # print('')
        # print('DUs', {k:len(v) for k,v in g.property('dispersal_units').iteritems()})
        # print('')
        # print('Lesions', {k:len(v) for k,v in g.property('lesions').iteritems()})

# Display results
from pylab import *
for k in green_leaves.iterkeys():
    plot(green_leaves[k])
    plot(healthy_leaves[k], '--')
    # plot(len_leaves[k])
show(False)
