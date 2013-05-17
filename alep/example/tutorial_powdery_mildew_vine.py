"""
Demonstrate a simulation of vine/oidium epidemics
"""

from alinea.alep.vine import Vine

from alinea.alep.inoculation import RandomInoculation
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.powdery_mildew import *
from alinea.alep.protocol import *

vine = Vine()
g0 = vine.setup_canopy(age=6)
vine.plot(g0)

fungus = powdery_mildew()
PowderyMildewDU.fungus=fungus
nb_dus = 100
dispersal_units = ([PowderyMildewDU(nb_spores=1, status="emitted") for i in range(nb_dus)])
inoculator=RandomInoculation()
initiate(g0,dispersal_units,inoculator, label='lf')

dt=1
controler = NoPriorityGrowthControl()
nsteps = 2

g=g0

for i in range(nsteps):
    g = vine.grow(g,dt)
    infect(g,dt)
    update(g,dt,controler)
