#Tutorial by Marc LABADIE for software paper rebuild echap loop with adel protocol

##########################################
# Warning : Hack DU add call !!!dispersal unit=nb_spore & 
#            position senecesence = senecesence_length
######################################

# 1. Get echap reconstruction

from pathlib import Path

from alinea.echap.architectural_reconstructions import EchapReconstructions
from alinea.adel.newmtg import move_properties

from alinea.alep.architecture import update_healthy_area, set_properties, get_leaves

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
#update_healthy_area(g, label = 'LeafElement')

#adel.plot(g)
#for i in range(steps):
#    canopy_age+=timestep
    
    # update canopy
#    newg = adel.setup_canopy(age=canopy_age)
#    adel.canopy_age=canopy_age
#    move_properties(g, newg)
#    g= newg
#    update_healthy_area(g, label = 'LeafElement')
#    adel.plot(g)
    
# 3. Disease development (septo3D)
import pandas
import random as rd

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
from alinea.alep.mini_models import leaf_wetness_rapilly
from alinea.alep.disease_operation import generate_stock_du

from alinea.astk.TimeControl import time_filter , rain_filter, IterWithDelays,time_control
from alinea.alep.protocol import infect , update, disperse

from alinea.alep.disease_outputs import plot_severity_by_leaf

#from alinea.alep.senescence import WheatSeptoriaPositionedSenescence

## Initialize the models for septoria
septoria = plugin_septoria()
inoculator = RandomInoculation()
growth_controler = NoPriorityGrowthControl()
infection_controler = BiotrophDUProbaModel()
emitter = SeptoriaRainEmission(domain_area=adel.domain_area)
transporter = SeptoriaRainDispersal()
washor = RapillyWashing()

# Read weather and adapt it to septoria (add wetness)
meteo_path = shared_data(alinea.alep, 'meteo05-06.txt')
weather = Weather(data_file=meteo_path)
weather.check(varnames=['wetness'], models={'wetness':lambda data: leaf_wetness_rapilly(data.rain, data.relative_humidity, data.PPFD)})
seq = pandas.date_range(start = "2005-10-01 01:00:00", end = "2006-07-01 01:00:00", freq='H')

every_h = time_filter(seq, delay=1)
every_24h = time_filter(seq, delay=24)
every_rain = rain_filter(seq, weather)

weather_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
wheat_timing = IterWithDelays(*time_control(seq, every_24h, weather.data))
septo_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))

# Choose source leaf in canopy 
# (Here the value of the leaf is known but it changes with another initialize_stand)
source_leaf = g.node(12)

for i, controls in enumerate(zip(weather_timing, wheat_timing, septo_timing, rain_timing)):
    weather_eval, wheat_eval, septo_eval, rain_eval = controls

    # Get weather for date and add it as properties on leaves
    if weather_eval.eval:
        #HACK MARC 2021: Considered that position scencescene is senescence_lenght
        if "position_senescence" not in g.property_names():
            g.add_property("position_senescence")

        g.property("position_senescence").update({k:v for k,v in g.property("senesced_length").items()}) 

        set_properties(g,label = 'LeafElement',
                       temp = weather_eval.value.temperature_air[0],
                       wetness = weather_eval.value.wetness[0],
                       relative_humidity = weather_eval.value.relative_humidity[0],
                       wind_speed = weather_eval.value.wind_speed[0])
                    
        
    if rain_eval.eval:
        set_properties(g,label = 'LeafElement',
                       rain_intensity = rain_eval.value.rain.mean(),
                       rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)
    
    # Grow wheat canopy
    if wheat_eval.eval:
        g = adel.grow(g, wheat_eval.value)

        # Note : The position of senescence goes back to its initial value after
        # a while for undetermined reason
        # --> temporary hack for keeping senescence position low when it is over
        positions = g.property('position_senescence')
        greens = g.property('is_green')
        areas = g.property('area')
        senesced_areas = g.property('senesced_area')
        leaves = get_leaves(g, label = 'LeafElement')
        vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
#        positions.update({vid:(0 if (positions[vid]==1 and not greens[vid]) or
#                                    (positions[vid]>0 and round(areas[vid],5)==round(senesced_areas[vid],5))
#                                    else positions[vid]) for vid in vids})
                                    
    # Develop disease
    if septo_eval.eval:
        # Update g for the disease
        #sen_model.find_senescent_lesions(g, label = 'LeafElement')
        update_healthy_area(g, label = 'LeafElement')
        
        # Possibly refill pool of initial inoculum to simulate differed availability
        if rain_eval.eval and i <= 700 and source_leaf.geometry!=None:
            dus = generate_stock_du(nb_dus=rd.randint(0,5), disease=septoria)
            try:
                source_leaf.dispersal_units += dus
            except:
                source_leaf.dispersal_units = dus
        
        # Update dispersal units and lesions
        infect(g, septo_eval.dt, infection_controler, label='LeafElement')
        update(g, septo_eval.dt, growth_controler, label='LeafElement')
        
    # Disperse and wash
    if rain_eval.eval:
        if rain_eval.value.rain.mean()>0:
            g = disperse(g, emitter, transporter, "septoria", label='LeafElement')
            # wash(g, washor, rain_eval.value.rain.mean(), label='LeafElement')
    
    if wheat_eval:
        scene = plot_severity_by_leaf(g, senescence=False, transparency=0.9)
 