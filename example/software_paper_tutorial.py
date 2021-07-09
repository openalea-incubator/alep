#Tutorial by Marc LABADIE for software paper rebuild echap loop with adel protocol

# 1. Get echap reconstruction

from pathlib import Path

from alinea.echap.architectural_reconstructions import EchapReconstructions
from alinea.adel.newmtg import move_properties

from alinea.alep.architecture import update_healthy_area

filename= 'echap_reconstruction.pckl'

#if not os.path.exists(filename):
if not Path(filename).exists():
    echap = EchapReconstructions()
    echap.save(filename=filename)
else:
    echap= EchapReconstructions.load(filename=filename)

adel = echap.get_reconstruction(name="Mercia", nplants=2,seed=0)

# 2. run simulation 

#init sumulation
timestep = 30 #Day degree
steps = 30



# 2.1 Grow wheat canopy and vizualized development

## init canopy
canopy_age=300
g = adel.setup_canopy(age=canopy_age)

## init alep
### Add the property 'healthy_area' on the leaves
update_healthy_area(g, label = 'LeafElement')

adel.plot(g)
for i in range(steps):
    canopy_age+=timestep
    
    # update canopy
    newg = adel.setup_canopy(age=canopy_age)
    adel.canopy_age=canopy_age
    move_properties(g, newg)
    g= newg
    update_healthy_area(g, label = 'LeafElement')
    adel.plot(g)
    
# 3. Disease development (septo3D)
import pandas

from openalea.deploy.shared_data import shared_data

from alinea.astk.Weather import Weather

import alinea.alep
from alinea.alep.septo3d_v2 import plugin_septoria
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.alep.dispersal_transport import SeptoriaRainDispersal
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUProbaModel

from alinea.alep.senescence import WheatSeptoriaPositionedSenescence

## Initialize the models for septoria
septoria = plugin_septoria()
inoculator = RandomInoculation()
growth_controler = NoPriorityGrowthControl()
infection_controler = BiotrophDUProbaModel()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
emitter = SeptoriaRainEmission(domain_area=domain_area)
transporter = SeptoriaRainDispersal()
washor = RapillyWashing()

# Read weather and adapt it to septoria (add wetness)
meteo_path = shared_data(alinea.alep, 'meteo05-06.txt')
weather = Weather(data_file=meteo_path)
weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
seq = pandas.date_range(start = "2005-10-01 01:00:00", end = "2006-07-01 01:00:00", freq='H')
