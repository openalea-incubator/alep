#Tutorial by Marc LABADIE for software paper rebuild echap loop with adel protocol

##########################################
# Warning : Hack DU add call !!!dispersal unit=nb_spore & 
#            position senecesence = senecesence_length
######################################

# 1. Get echap reconstruction

from pathlib import Path
import pandas
import random as rd

from alinea.echap.architectural_reconstructions import EchapReconstructions

import alinea.alep
from alinea.alep.architecture import set_properties
from alinea.alep.septo3d_v2 import plugin_septoria
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.alep.dispersal_transport import SeptoriaRainDispersal
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.mini_models import leaf_wetness_rapilly, linear_degree_days
from alinea.alep.disease_operation import DU_Generator
from alinea.alep.protocol import infect , update, disperse,external_contamination


from alinea.echap.weather_data import get_weather
from alinea.echap.interception_leaf import pesticide_applications, InterceptModel
from alinea.echap.interfaces import record as do_record, pesticide_interception, pesticide_surfacic_decay, pesticide_penetrated_decay, pesticide_efficacy#to avoid confusion with numpy record
from alinea.echap.microclimate_leaf import microclimate_leaf
from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.milne_leaf import PenetratedDecayModel
from alinea.echap.pesticide_efficacy import PesticideEfficacyModel
from alinea.echap.contamination import SimpleContamination,SimpleSoilInoculum
from alinea.echap.tests_nodes import plot_pesticide

from alinea.astk.TimeControl import time_filter , rain_filter, IterWithDelays,time_control,date_filter,DegreeDayModel, thermal_time_filter
from alinea.astk.Weather import Weather, linear_degree_days




#1. Wheat architecture costruction from Echap data
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
canopy_age=150
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
weather=get_weather(variety='Mercia')
weather.check(["degree_days",'temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'],models={'degree_days':linear_degree_days})
    
#weather.check(varnames=['wetness'], models={'wetness':lambda data: leaf_wetness_rapilly(data.rain, data.relative_humidity, data.PPFD)})
periods = 500 # 5000 pour le cycle complet
seq = pandas.date_range(start = "2010-11-02", periods=periods, freq='H')  
seq_bid = pandas.date_range(start = "2010-10-15", periods=periods+7000, freq='H')
applications = """date,dose, product_name
2010-11-02 00:00:00, 0, bug
2010-11-02 1:00:00, 1, Opus
2010-11-02 2:00:00, 0, bug
2010-04-29 10:00:00, 1, Opus
2010-04-29 11:00:00, 0, bug
2010-05-22 10:00:00, 1, Opus
2010-05-22 11:00:00, 0, bug
"""

pest_calendar = pesticide_applications(applications)
productsDB={'Opus': {'Epoxiconazole': 125}, 'Banko 500': {'Chlorothalonil': 500}}
interception = InterceptModel(productsDB)
Milne_efficacy = PesticideEfficacyModel()
#pgen, adel, domain, domain_area, convUnit, nplants = reconst(name=name, dTT_stop=dTT_stop, original=original, n=nplants)
#tx = pgen['dynT_user'].a_cohort[0]
SspoSol = 0.01
recorder = LeafElementRecorder()
Milne_leaf = PenetratedDecayModel()
TTmodel = DegreeDayModel(Tbase = 0)
every_rain = rain_filter(seq, weather)
every_h = time_filter(seq, delay = 6)
every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 15)
every_pest = date_filter(seq, pest_calendar)

rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
canopy_timing = IterWithDelays(*time_control(seq, every_dd, weather.data))
doses_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
pest_timing = IterWithDelays(*time_control(seq, every_pest, pest_calendar))

def update_pesticides(g, weather_data):
    #g = pesticide_surfacic_decay(g, PearlLeaf, weather_data)
    g = pesticide_penetrated_decay(g, Milne_leaf, weather_data)
    g = pesticide_efficacy(g, Milne_efficacy, weather_data)
    return g

def simple_contamination(g, weather_data, level, domain, domain_area, convUnit):
    inoc = SimpleSoilInoculum(DU_generator=DU_Generator(disease=septoria), sporulating_fraction =level, domain_area = domain_area, convUnit=convUnit)
    contaminator = SimpleContamination(domain = domain, domain_area =domain_area, convUnit=convUnit)
    g = external_contamination(g, inoc, contaminator, weather_data)
    return g  

# Choose source leaf in canopy 
# (Here the value of the leaf is known but it changes with another initialize_stand)
#source_leaf = g.node(12)

for i,controls in enumerate(zip(canopy_timing, doses_timing, rain_timing, pest_timing)):
    canopy_iter, doses_iter, rain_iter, pest_iter = controls
    
    if canopy_iter.eval:
        print('--', canopy_iter.value.index[-1], '--')
        print('-- update microclimate / canopy --')
        g = adel.grow(g, canopy_iter.value)
        #_=rain_and_light_star(g, light_sectors = '1', domain=adel.domain, convUnit=adel.convUnit)
        _=update(g, canopy_iter.dt, growth_controler)
        #_=do_record(g, canopy_iter.value, recorder)
    if pest_iter.eval:
        print('--', pest_iter.value.index[-1], '--')
        print('-- update microclimate / pesticide --')
        
        #_=pesticide_intercept(g, pest_iter.value)
        _=pesticide_interception(g, interception, pest_iter.value, label='LeafElement')
        #_=plot_pesticide(g)
    if doses_iter.eval:
        print('--', doses_iter.value.index[-1], '--')
        print('-- update microclimate / doses --')
        #_=microclimate_leaf(g, doses_iter.value, domain = adel.domain, convUnit = adel.convUnit)
        _=update_pesticides(g, doses_iter.value)
        # plot_pesticide(g)
        #_=do_record(g, doses_iter.value, recorder, header={'iter':i, 'TT':adel.canopy_age})
    if rain_iter.eval:
        print('--', rain_iter.value.index[-1], '--')
        print('-- rain --')

        wdata = rain_iter.value
        
        # Get weather for date and add it as properties on leaves
        set_properties(g,label = 'LeafElement',
            temperature_sequence = wdata.temperature_air.tolist(),
            relative_humidity_sequence= wdata.relative_humidity.tolist(),
            wetness_sequence = leaf_wetness_rapilly(rain=wdata.rain,relative_humidity=wdata.relative_humidity,PPFD=wdata.PPFD),
            #dd_sequence = wdata.degree_days.tolist()            
            )

#        set_properties(g,label = 'LeafElement',
#            rain_intensity = wdata.rain.mean(),
#            rain_duration = len(wdata.rain) if wdata.rain.sum() > 0 else 0.)

        g = disperse(g, emitter, transporter, "septoria", label='LeafElement')
        #dispersion(g, wdata, domain, domain_area, convUnit)
        _=simple_contamination(g,wdata, SspoSol, adel.domain, adel.domain_area, adel.convUnit)
        _=infect(g, rain_iter.dt)


    # # Get weather for date and add it as properties on leaves
    # if weather_eval.eval:
    #     #HACK MARC 2021: Considered that position scencescene is senescence_lenght
    #     if "position_senescence" not in g.property_names():
    #         g.add_property("position_senescence")

    #     g.property("position_senescence").update({k:v for k,v in g.property("senesced_length").items()}) 

    #     set_properties(g,label = 'LeafElement',
    #                    temp = weather_eval.value.temperature_air[0],
    #                    wetness = weather_eval.value.wetness[0],
    #                    relative_humidity = weather_eval.value.relative_humidity[0],
    #                    wind_speed = weather_eval.value.wind_speed[0])
                    
        
    # if rain_eval.eval:
    #     set_properties(g,label = 'LeafElement',
    #                    rain_intensity = rain_eval.value.rain.mean(),
    #                    rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)
    
    # Grow wheat canopy
#     if wheat_eval.eval:
#         g = adel.grow(g, wheat_eval.value)

#         # Note : The position of senescence goes back to its initial value after
#         # a while for undetermined reason
#         # --> temporary hack for keeping senescence position low when it is over
#         positions = g.property('position_senescence')
#         greens = g.property('is_green')
#         areas = g.property('area')
#         senesced_areas = g.property('senesced_area')
#         leaves = get_leaves(g, label = 'LeafElement')
#         vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
# #        positions.update({vid:(0 if (positions[vid]==1 and not greens[vid]) or
# #                                    (positions[vid]>0 and round(areas[vid],5)==round(senesced_areas[vid],5))
# #                                    else positions[vid]) for vid in vids})
                                    
#     # Develop disease
#     if septo_eval.eval:
#         # Update g for the disease
#         #sen_model.find_senescent_lesions(g, label = 'LeafElement')
#         update_healthy_area(g, label = 'LeafElement')
        
#         # Possibly refill pool of initial inoculum to simulate differed availability
#         if rain_eval.eval and i <= 700 and source_leaf.geometry!=None:
#             dus = generate_stock_du(nb_dus=rd.randint(0,5), disease=septoria)
#             try:
#                 source_leaf.dispersal_units += dus
#             except:
#                 source_leaf.dispersal_units = dus
        
#         # Update dispersal units and lesions
#         infect(g, septo_eval.dt, infection_controler, label='LeafElement')
#         update(g, septo_eval.dt, growth_controler, label='LeafElement')
        
#     # Disperse and wash
#     if rain_eval.eval:
#         if rain_eval.value.rain.mean()>0:
#             g = disperse(g, emitter, transporter, "septoria", label='LeafElement')
#             # wash(g, washor, rain_eval.value.rain.mean(), label='LeafElement')
    
#     if wheat_eval:
#         scene = plot_severity_by_leaf(g, senescence=False, transparency=0.9)
 