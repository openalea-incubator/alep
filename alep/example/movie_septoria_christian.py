
# Imports for weather
from datetime import datetime
from alinea.alep.alep_weather import get_septoria_weather
import alinea.septo3d
from openalea.deploy.shared_data import shared_data

# Imports for wheat
from alinea.alep.wheat import initialize_stand, find_blade_id, find_leaf_ids
from alinea.astk.plant_interface import grow_canopy
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.disease_operation import generate_stock_du
from alinea.aled.epidemics import SeptoriaSepto3D as septo




meteo_path = shared_data(alinea.septo3d, 'meteo00-01.txt')
# Read weather and adapt it to septoria (add wetness)
weather = get_septoria_weather(data_file=meteo_path, sep='\t')

# Choose dates of simulation and initialize the value of date
start_date = datetime(2000, 10, 1, 1, 00, 00)
end_date = datetime(2001, 07, 01, 00, 00)




# Initialize a wheat canopy
# g, wheat, domain_area, domain = initialize_stand(age=0., length=0.5, 
                                                 # width=0.5, sowing_density=250,
                                                 # plant_density=250, inter_row=0.12, seed=1)

g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, 
                                                width=0.2, sowing_density=150,
                                                plant_density=150, inter_row=0.12)
                                                
# Define the schedule of calls for each model
nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
weather_timing = TimeControl(delay=1, steps=nb_steps)
wheat_timing = TimeControl(delay=24, steps=nb_steps, model=wheat, weather=weather, start_date=start_date)
septo_timing = TimeControl(delay=1, steps=nb_steps)
timer = TimeControler(weather=weather_timing, wheat=wheat_timing, disease = septo_timing)

for t in timer:
    date = (weather.next_date(t['weather'].dt, date) if date!=None else start_date)
    print(date)
    grow_canopy(g,wheat,t['wheat'])