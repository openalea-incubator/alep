"""
Demonstrate a simulation of vine/oidium epidemics
"""

from alinea.alep.vine import Vine
from alinea.astk.TimeControl import *

from alinea.alep.inoculation import RandomInoculation
from alinea.alep.growth_control import GrowthControlVineLeaf
from alinea.alep.dispersal import RandomDispersal
from alinea.alep.powdery_mildew import *
from alinea.alep.protocol import *
from alinea.adel.mtg_interpreter import plot3d

from alinea.alep.disease_operation import *
from alinea.alep.architecture import *
from alinea.alep.alep_color import alep_colormap, green_white

def update_plot(g):
    # Count lesions by id & add it as MTG property ####################################
    nb_lesions_by_leaf = count_lesions_by_leaf(g, label = 'lf')
    set_property_on_each_id(g, 'nb_lesions', nb_lesions_by_leaf, label = 'lf')
                       
    # Visualization ###################################################################
    g = alep_colormap(g, 'nb_lesions', cmap=green_white(levels=10), lognorm=False)
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
set_properties(g0,label = 'lf',
               surface=5., healthy_surface=5., position_senescence=None)

fungus = powdery_mildew()
PowderyMildewDU.fungus=fungus
nb_dus = 100
dispersal_units = ([PowderyMildewDU(nb_spores=1, status="emitted") for i in range(nb_dus)])
inoculator=RandomInoculation()
initiate(g0,dispersal_units,inoculator, label='lf')


controler = GrowthControlVineLeaf()
dispersor = RandomDispersal()
# nsteps = 49
nsteps = 100

vine_timing = TimeControl(delay=24, steps = nsteps)
mildew_timing = TimeControl(delay =1, steps = nsteps)
plot_timing = TimeControl(delay=10, steps = nsteps)
timer = TimeControler(vine = vine_timing, disease = mildew_timing, ploting = plot_timing)


g=g0
scene = vine.generate_scene(g)

for t in timer:
    print(timer.numiter)
    set_properties_on_new_leaves(g,label = 'lf',
                             surface=5., healthy_surface=5.,
                             position_senescence=None)
    set_properties(g,label = 'lf',
                    wetness=True,
                    temp=22.,
                    rain_intensity=0.,
                    rain_duration=0.,
                    relative_humidity=85.,
                    wind_speed=0.5)
                        
    infect(g, t['disease'].dt, label='lf')
    update(g, t['disease'].dt, controler, label='lf')
    disperse(g, scene, dispersor, "powdery_mildew", label='lf')  
    
    g = vine.grow(g,t['vine'])
    scene = vine.generate_scene(g)
    if t['ploting'].dt > 0:
        print('ploting...')
        update_plot(g)
