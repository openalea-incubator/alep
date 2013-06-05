""" Tutorial for septoria"""
import random as rd

from alinea.alep.wheat_examples import adel_mtg, adel_mtg2, adel_one_leaf
from alinea.alep.architecture import *
from alinea.alep.disease_operation import *
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.protocol import *
from alinea.alep.alep_color import alep_colormap, green_yellow_red

from openalea.vpltk import plugin

from alinea.astk.TimeControl import *

#######################################################################################

def update_plot(g):
    # Count lesion surfaces by id & add it as MTG property 
    surface_lesions_by_leaf = count_lesion_surfaces_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'surface_lesions', surface_lesions_by_leaf, label = 'LeafElement')
                       
    # Visualization
    g = alep_colormap(g, 'surface_lesions', cmap=green_yellow_red(levels=10), lognorm=False)
    scene = plot3d(g)
    Viewer.display(scene)
    return scene

#######################################################################################
# Define a plant or canopy
g = adel_mtg2()

# Add missing properties needed for the simulation
# The simulation requires the following properties on leaf elements:
#   - 'surface': total surface of the leaf element
#   - 'healthy_surface': surface of the leaf element without lesion or senescence
#   - 'position_senescence': position of the senescence on blade axis
set_properties(g,label = 'LeafElement',
               surface=5., healthy_surface=5., position_senescence=None)
               
# discover the disease implementation for septoriose
diseases=plugin.discover('alep.disease')
#septoria_classes = [kls for kls in diseases if kls.load().name == 'septoria']
# or 
septoria = diseases['septoria_exchanging_rings'].load()

# Create a pool of dispersal units (DU)
nb_du = 10
dispersal_units = generate_stock_du(nb_du, disease=septoria)

# Distribute the DU 
inoculator = RandomInoculation()
initiate(g, dispersal_units, inoculator)
# g = color.colormap(g,'length',lognorm=False)
scene = plot3d(g)
Viewer.display(scene)

# Preparation of the simulation loop
controler = NoPriorityGrowthControl()
dispersor = Septo3DSplash(reference_surface=1./200)
nsteps = 1000
wheat_timing = TimeControl(delay=1, steps = nsteps)
septo_timing = TimeControl(delay=1, steps = nsteps)
plot_timing = TimeControl(delay=24, steps = nsteps)
timer = TimeControler(wheat = wheat_timing, disease = septo_timing, ploting = plot_timing)


for t in timer:
    print(timer.numiter)
    set_properties_on_new_leaves(g,label = 'LeafElement',
                             surface=5., healthy_surface=5.,
                             position_senescence=None)
    
    # Rain setting:
    if timer.numiter>=400 and timer.numiter%100==0:
        rain_intensity=1.
        rain_duration=2.
    else:
        rain_intensity=0.
        rain_duration=0.
    
    set_properties(g,label = 'LeafElement',
                    wetness=True,
                    temp=22.,
                    rain_intensity=rain_intensity,
                    rain_duration=rain_duration,
                    relative_humidity=85.,
                    wind_speed=0.5)
                        
    infect(g, t['disease'].dt, label='LeafElement')
    update(g, t['disease'].dt, controler, label='LeafElement')

    if rain_intensity>0:
        disperse(g, scene, dispersor, "septoria", label='LeafElement')  

    if t['ploting'].dt > 0:
        print('ploting...')
        scene = update_plot(g)
