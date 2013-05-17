"""
Demonstrate a simulation of vine/oidium epidemics
"""

from alinea.alep.vine import Vine

from alinea.alep.inoculation import RandomInoculation
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.dispersal import RandomDispersal
from alinea.alep.architecture import *
from alinea.alep.powdery_mildew import *
from alinea.alep.protocol import *
from alinea.adel.mtg_interpreter import plot3d

vine = Vine()
g0 = vine.setup_canopy(age=6)
vine.plot(g0)

# Add missing properties needed for the simulation
# The simulation requires the following properties on leaf elements:
#   - 'surface': total surface of the leaf element (in cm2)
#   - 'healthy_surface': surface of the leaf element without lesion or senescence (in cm2)
#   - 'age': age of the leaf (in hours)
#   - 'position_senescence': position of the senescence on blade axis
set_properties(g0,label = 'lf',
               surface=5, healthy_surface=5, position_senescence=None)

fungus = powdery_mildew()
PowderyMildewDU.fungus=fungus
nb_dus = 100
dispersal_units = ([PowderyMildewDU(nb_spores=1, status="emitted") for i in range(nb_dus)])
inoculator=RandomInoculation()
initiate(g0,dispersal_units,inoculator, label='lf')

dt=1
controler = NoPriorityGrowthControl()
dispersor = RandomDispersal()
nsteps = 2

g=g0

for i in range(nsteps):
    g = vine.grow(g,dt)
    set_properties(g,label = 'lf',
                    wetness=True,
                    temp=22.,
                    rain_intensity=0.,
                    rain_duration=0.,
                    relative_humidity=85.,
                    wind_speed=0.2)
    infect(g,dt)
    update(g,dt,controler)
    
    scene = plot3d(g)
    disperse(g, scene, dispersor, "PowderyMildew")
