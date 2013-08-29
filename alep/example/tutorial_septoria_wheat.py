""" Tutorial for septoria"""
import random as rd
import numpy as np
import pandas
from datetime import datetime

from alinea.astk.plantgl_utils import *

from alinea.adel.stand.stand import agronomicplot
from alinea.adel.astk_interface import AdelWheat
from alinea.astk.plant_interface import *
#from alinea.alep.wheat import adel_mtg2
from alinea.alep.architecture import (set_properties, set_healthy_area,
                                      set_property_on_each_id, get_leaves, 
                                      get_leaves)
from alinea.alep.disease_operation import *
from alinea.alep.disease_outputs import (compute_lesion_areas_by_leaf,
                                        compute_severity_by_leaf)
from alinea.alep.inoculation import InoculationFirstLeaves, RandomInoculation
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.protocol import *
from alinea.alep.alep_color import alep_colormap, green_yellow_red

from alinea.alep.alep_weather import add_wetness
from alinea.weather.global_weather import *

from openalea.vpltk import plugin

from alinea.astk.TimeControl import *

rd.seed(0)
np.random.seed(0)

# Useful functions ########################################################################
def update_plot(g):       
    # Count lesion surfaces by id & add it as MTG property 
    surface_lesions_by_leaf = compute_lesion_areas_by_leaf(g, label = 'LeafElement')
    # print('Lesion Surface ', surface_lesions_by_leaf)
    # print('')
    # vids = get_leaves(g)
    # print('Position Senescence ', {vid:g.node(vid).position_senescence for vid in vids})
    # print('')
    # print('Healhty Surface ', g.property('healthy_area'))
    set_property_on_each_id(g, 'surface_lesions', surface_lesions_by_leaf, label = 'LeafElement')

    # Compute severity by leaf
    severity_by_leaf = compute_severity_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'severity', severity_by_leaf, label = 'LeafElement')

    # Visualization
    g = alep_colormap(g, 'surface_lesions', cmap=green_yellow_red(levels=100), lognorm=False)
    # TODO : Normalize the colormap between 0 and 1 /!\
    scene = plot3d(g)
    Viewer.display(scene)
    return scene
    

def plugin_septoria():
    diseases=plugin.discover('alep.disease')
    #septoria_classes = [kls for kls in diseases if kls.load().name == 'septoria']

    try:
        septoria = diseases['septoria_exchanging_rings'].load()
        # septoria = diseases['septoria_continuous'].load()
        # septoria = diseases['septoria_with_rings'].load()
    except KeyError:
        from alinea.alep.septoria_exchanging_rings import Disease
        septoria=Disease()
    return septoria

# Initiation ##############################################################################
# Define a plant or canopy
nplants, positions, domain, domain_area = agronomicplot(0.1, 0.2, 150, 150, 0.12) 
#hack (to be changed in adel)
domain_area_cm2 = domain_area * 10000
wheat = AdelWheat(nplants=nplants, positions = positions)
g,_ = new_canopy(wheat,age=100)
geometries = g.property('geometry')
lai = get_lai(geometries, domain_area_cm2)
# Add the property 'healthy_area' on the leaves
set_healthy_area(g, label = 'LeafElement')

# from alinea.alep.wheat import adel_mtg2
# g = adel_mtg2()

# Add missing properties needed for the simulation
# The simulation requires the following properties on leaf elements:
#   - 'surface': total surface of the leaf element
#   - 'healthy_surface': surface of the leaf element without lesion or senescence
#   - 'position_senescence': position of the senescence on blade axis
# set_properties(g,label = 'LeafElement',
              # surface=5., healthy_surface=5., position_senescence=None)
               
# discover the disease implementation for septoriose
septoria = plugin_septoria()

# Create a pool of dispersal units (DU)
nb_dus = 10
# nb_dus = 100
dispersal_units = generate_stock_du(nb_dus, disease=septoria)

# Distribute the DU 
inoculator = RandomInoculation()
# initiate(g, dispersal_units, inoculator)

# Preparation of the simulation loop #####################################################
# Call models that will be used in disease interface
controler = NoPriorityGrowthControl()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
dispersor = Septo3DSplash(reference_surface=domain_area)

# Choose dates of simulation
start_date = datetime(2000, 10, 1, 1, 00, 00)
end_date = datetime(2001, 04, 01, 00, 00)
#end_date = datetime(2000, 12, 31, 00, 00)
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

AA = []
GA = []
HA = []
LA = []                      
leaves = get_leaves(g, label = 'LeafElement')
green_leaves = {leaf:[] for leaf in leaves}
healthy_leaves = {leaf:[] for leaf in leaves}
sen_leaves = {leaf:[] for leaf in leaves}

for t in timer:
    #print(timer.numiter)
    
    # Not needed if wheat model is complete
    #set_properties_on_new_leaves(g,label = 'LeafElement',
    #                         surface=5., healthy_surface=5.,
    #                         position_senescence=None)
    
    # position_senescence = g.property('position_senescence')
    # if len(position_senescence)==0:
        # g.remove_property('position_senescence')
    # set_properties(g,label = 'LeafElement', position_senescence=1., green_area=5.)
    
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
    if timer.numiter%10 == 0 and timer.numiter >= 500 and timer.numiter <= 1000:
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
        
    if False and t['plotting'].dt > 0:
        print('area', g.property('area'))
        print('')
        print('green area', g.property('green_area'))
        print('')
        print('healthy area', g.property('healthy_area'))
        print('')
        print('Lesion area', les_area)
        print('')
        print('DUs', {k:len(v) for k,v in g.property('dispersal_units').iteritems()})
        print('')
        print('Lesions', {k:len(v) for k,v in g.property('lesions').iteritems()})
        
# Display results
from pylab import *
for k in green_leaves.iterkeys():
    plot(green_leaves[k])
    plot(healthy_leaves[k], '--')
show(False)
