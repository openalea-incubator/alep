""" Execute the first simulation of the paper : simulate a septoria epidemics on wheat. """

# Useful imports
import random as rd
import numpy as np
import pandas
from alinea.astk.TimeControl import *
from alinea.alep.alep_color import alep_colormap, green_yellow_red

# Imports for weather
from datetime import datetime
from alinea.alep.alep_weather import get_septoria_weather

# Imports for wheat
from alinea.alep.wheat import initialize_stand, find_blade_id, find_leaf_ids
from alinea.astk.plant_interface import grow_canopy
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.disease_outputs import LeafInspector
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence

# Initialization #####################################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Choose dates of simulation and initialize the value of date
start_date = datetime(2000, 10, 1, 1, 00, 00)
# end_date = datetime(2000, 12, 31, 00, 00)
end_date = datetime(2001, 07, 01, 00, 00)
date = None

# Read weather and adapt it to septoria (add wetness)
weather = get_septoria_weather(data_file='meteo01.csv')

# Initialize a wheat canopy
g, wheat, domain_area = initialize_stand(age=0., length=0.1, 
                                        width=0.2, sowing_density=150,
                                        plant_density=150, inter_row=0.12)

# Initialize the models for septoria
septoria = plugin_septoria()
inoculator = RandomInoculation()
growth_controler = NoPriorityGrowthControl()
infection_controler = BiotrophDUPositionModel()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
emitter = SeptoriaRainEmission()
transporter = Septo3DSplash(reference_surface=domain_area)

# Define the schedule of calls for each model
nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
weather_timing = TimeControl(delay=1, steps=nb_steps)
wheat_timing = TimeControl(delay=24, steps=nb_steps, model=wheat, weather=weather, start_date=start_date)
septo_timing = TimeControl(delay=1, steps=nb_steps)
timer = TimeControler(weather=weather_timing, wheat=wheat_timing, disease = septo_timing)

# Call leaf inspectors for target blades (top 3)
inspectors = {}
for rank in range(1,4):
    inspectors[rank] = LeafInspector(g, blade_id=find_blade_id(g, leaf_rank = rank, only_visible=False))

# Simulation #########################################################
for t in timer:
    # print(timer.numiter)
    # Update date
    date = (weather.next_date(t['weather'].dt, date) if date!=None else start_date)
    print(date)
    
    # Get weather for date and add it as properties on leaves
    _, data = weather.get_weather(t['weather'].dt, date)
    set_properties(g,label = 'LeafElement',
                    wetness = data.wetness.values[0],
                    temp = data.temperature_air.values[0],
                    rain_intensity = data.rain.values[0],
                    rain_duration = data.rain_duration.values[0],
                    relative_humidity = data.relative_humidity.values[0],
                    wind_speed = data.wind_speed.values[0])

    # Grow wheat canopy
    grow_canopy(g,wheat,t['wheat'])
    update_healthy_area(g, label = 'LeafElement')
    # Note : The position of senescence goes back to its initial value after
    # a while for undetermined reason
    # --> temporary hack for keeping senescence position low when it is over
    positions = g.property('position_senescence')
    are_green = g.property('is_green')
    leaves = get_leaves(g, label = 'LeafElement')
    vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
    positions.update({vid:(0 if positions[vid]==1 and not are_green[vid] else positions[vid]) 
                      for vid in vids})
    
    # Develop disease
    if timer.numiter%10 == 0 and timer.numiter <= 500:
        # Refill pool of initial inoculum to simulate differed availability
        dispersal_units = generate_stock_du(nb_dus=10, disease=septoria)
        initiate(g, dispersal_units, inoculator)
    infect(g, t['disease'].dt, infection_controler, label='LeafElement')
    update(g, t['disease'].dt, growth_controler, sen_model, label='LeafElement')
    if data.dispersal_event.values[0]==True:
        disperse(g, emitter, transporter, "septoria", label='LeafElement')
        
    # Save outputs
    for inspector in inspectors.itervalues():
        inspector.update_variables(g)