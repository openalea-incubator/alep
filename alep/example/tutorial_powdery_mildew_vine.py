"""
Demonstrate a simulation of vine/oidium epidemics
"""

from alinea.alep.vine import Vine
from alinea.astk.TimeControl import *

from alinea.alep.inoculation import RandomInoculation
from alinea.alep.growth_control import GrowthControlVineLeaf
from alinea.alep.dispersal import RandomDispersal, PowderyMildewWindDispersal
from alinea.alep.powdery_mildew import *
from alinea.alep.protocol import *
from alinea.adel.mtg_interpreter import plot3d

from alinea.alep.disease_operation import *
from alinea.alep.disease_outputs import *
from alinea.alep.architecture import (get_leaves, set_properties, set_property_on_each_id, 
                                      set_properties_on_new_leaves, add_area_topvine)
from alinea.alep.alep_color import alep_colormap, green_white

def update_plot(g):
    # Count lesions by id & add it as MTG property ####################################
    nb_lesions_by_leaf = count_lesions_by_leaf(g)
    set_property_on_each_id(g, 'nb_lesions', nb_lesions_by_leaf, label = 'lf')
                       
    # Visualization ###################################################################
    g = alep_colormap(g, 'nb_lesions', cmap=green_white(levels=10), lognorm=False, 
                        zero_to_one=False, vmax=100)
    labels = g.property('label')
    trunk_ids = [k for k,l in labels.iteritems() 
                if l.startswith('tronc') or l.startswith('en')]
    brown = (100,70,30)
    for id in trunk_ids:
        trunk = g.node(id)
        trunk.color = brown
    scene = plot3d(g)
    Viewer.display(scene)

vine = Vine()
g0 = vine.setup_canopy(age=6)
# vine.plot(g0)

# Add missing properties needed for the simulation
# The simulation requires the following properties on leaf elements:
#   - 'surface': total surface of the leaf element (in cm2)
#   - 'healthy_surface': surface of the leaf element without lesion or senescence (in cm2)
#   - 'age': age of the leaf (in decimal_days)
#   - 'position_senescence': position of the senescence on blade axis
add_area_topvine(g0, conversion_factor=1000., label = 'lf')
set_properties(g0, label = 'lf', position_senescence=None)
# set_properties(g0, label = 'lf', surface=5., healthy_surface=5., position_senescence=None)

fungus = powdery_mildew()
PowderyMildewDU.fungus=fungus
nb_dus = 110
dispersal_units = ([PowderyMildewDU(nb_spores=1, status="emitted") for i in range(nb_dus)])
inoculator=RandomInoculation()
initiate(g0,dispersal_units,inoculator, label='lf')

controler = GrowthControlVineLeaf()
# dispersor = RandomDispersal()
dispersor = PowderyMildewWindDispersal()
# nsteps = 49
nsteps = 1000

vine_timing = TimeControl(delay=24, steps = nsteps)
mildew_timing = TimeControl(delay =1, steps = nsteps)
plot_timing = TimeControl(delay=24, steps = nsteps)
timer = TimeControler(vine = vine_timing, disease = mildew_timing, plotting = plot_timing)

# TEMP
# initial_vids = get_leaves(g0, leaf_name='lf')

g=g0
scene = vine.generate_scene(g)

def step(t):
    global g
    global scene
    print(timer.numiter)
    
    set_properties(g,label = 'lf',
                    wetness=True,
                    temp=22.,
                    rain_intensity=0.,
                    rain_duration=0.,
                    relative_humidity=85.,
                    wind_speed=0.5,
                    wind_direction=(1.,0.,0.))
                        
    infect(g, t['disease'].dt, label='lf')
    update(g, t['disease'].dt, controler, label='lf')
    disperse(g, dispersor, "powdery_mildew", label='lf')  
    
    g = vine.grow(g,t['vine'])
    set_properties_on_new_leaves(g,label = 'lf', position_senescence=None)
    add_area_topvine(g, conversion_factor=1000., label = 'lf')
    # TEMP
    # new_vids = get_leaves(g, leaf_name='lf')
    # if len(new_vids)>len(initial_vids):
        # raise Exception('')
        
    scene = plot3d(g)
    if t['plotting'].dt > 0:
        print('plotting...')
        update_plot(g)


for i in range(nsteps-100):
    t = timer.next()
    step(t)
