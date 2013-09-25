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
from alinea.alep.inoculation import RandomInoculation
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence
from alinea.alep.disease_outputs import LeafInspector
from alinea.alep.disease_outputs import compute_total_necrotic_area

# Modified senescence model ##########################################
def generate_lesions_with_senescence_treshold(nb_lesions=1, senescence_treshold=250.):
    """ Generate a pool of lesion with a modified senescence reponse.
    
    Senescence responds to an age treshold independent of lesion status.
    
    Parameters
    ----------
    nb_lesions: int
        Number of lesions in stock
    senescence_treshold: float
        Age at which the lesion is affected by senescence
        
    Returns
    -------
    lesions: list
        List of lesions with senescence response according to a treshold
    """
    
    # New model of senescence response
    def senescence_response_new(self):
        """ Find surfaces killed by senescence.
        
        Senescence affects surfaces whose age is below the age treshold.
        """
        from math import floor
        treshold = self.senescence_treshold
        f = self.fungus
        rings = self.surface_rings
        surface = self.surface
        delta_ring = self.delta_age_ring
        width_ring = f.delta_age_ring
        time_to_chlo = f.degree_days_to_chlorosis
        if treshold < time_to_chlo:
            raise Exception('Do not use a treshold below time to chlorosis')
        ddday = self.ddday
        ddday_sen = self.ddday_before_senescence
        age_dday = self.age_dday
        
        # Stop growth
        self.disable_growth()
        
        if age_dday <= treshold:
            # Everything is affected
            self.disable()
            self.first_ring.disable()
            surface_dead = self.surface
            surface_non_senescent = 0.
            if len(rings)>0:
                rings[:]=0.
        else:
            treshold -= time_to_chlo
            diff = delta_ring - treshold
            nb_full_rings = floor(diff/width_ring)
            if len(rings)>0:
                if nb_full_rings == 0:
                    # The treshold is inside the current ring, which is not full
                    portion_non_senescent = (age_dday-treshold)/(age_dday+width_ring-delta_ring)
                    surface_non_senescent = portion_non_senescent*rings[-1]
                    rings[-1] *= portion_non_senescent
                    rings[:-1] = 0.
                else:
                    # The treshold is inside an already full ring
                    portion_non_senescent = (diff%width_ring)/width_ring
                    if portion_non_senescent > 0.:
                        surface_non_senescent = (sum(rings[-nb_full_rings:]) + 
                                                 portion_non_senescent*rings[-(nb_full_rings+1)])
                        rings[-(nb_full_rings+1)] *= portion_non_senescent
                        rings[:-(nb_full_rings+1)] = 0.
                    else:
                        surface_non_senescent = sum(rings[-nb_full_rings:])
                        rings[:-nb_full_rings] = 0.
                        # update surfaces rings ...
                surface_non_senescent += self.first_ring.surface
                surface_dead = self.surface - surface_non_senescent

        self.surface_alive = surface_non_senescent
        self.surface_dead = surface_dead
        
        # Complete the age of the lesion up to the end of time step
        self.ddday = ddday - ddday_sen
        self.age_dday += ddday - ddday_sen
        
        # Manage first ring
        self.first_ring.update(lesion=self)
        # Manage the other rings
        if len(self.surface_rings)>0:
            if self.age_dday-f.degree_days_to_chlorosis > self.delta_age_ring:
                self.surface_rings = np.append(self.surface_rings, 0.)
                self.delta_age_ring += f.delta_age_ring
            self.exchange_surfaces()

        # Update stock of spores
        if self.status==f.SPORULATING:
            self.update_stock()
    
    # Generation of lesion pool
    septoria = plugin_septoria()
    LesionKlass = septoria.lesion()
    LesionKlass.senescence_response = senescence_response_new
    lesions = [LesionKlass(nb_spores = 1.) for i in range(nb_lesions)]
    for les in lesions:
        les.senescence_treshold = senescence_treshold
    return lesions

# Initialization #####################################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Choose dates of simulation and initialize the value of date
start_date = datetime(2000, 10, 1, 1, 00, 00)
end_date = datetime(2000, 12, 31, 00, 00)
# end_date = datetime(2001, 07, 01, 00, 00)
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
controler = NoPriorityGrowthControl()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
dispersor = Septo3DSplash(reference_surface=domain_area)

# Define the schedule of calls for each model
nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
weather_timing = TimeControl(delay=1, steps=nb_steps)
wheat_timing = TimeControl(delay=24, steps=nb_steps, model=wheat, weather=weather, start_date=start_date)
septo_timing = TimeControl(delay=1, steps=nb_steps)
timer = TimeControler(weather=weather_timing, wheat=wheat_timing, disease = septo_timing)

# Call leaf inspectors for target blades (top 3)
necrosis = []
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
    infect(g, t['disease'].dt, label='LeafElement')
    update(g, t['disease'].dt, controler, sen_model, label='LeafElement')
    if data.dispersal_event.values[0]==True:
        disperse(g, dispersor, "septoria", label='LeafElement')
        
    # Save outputs
    necrosis.append(compute_total_necrotic_area(g, label='LeafElement'))
    for inspector in inspectors.itervalues():
        inspector.update_variables(g)