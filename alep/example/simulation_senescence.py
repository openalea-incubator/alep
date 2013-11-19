""" Simulate a septoria epidemics on wheatwith a given response to senescence. """

# Useful imports
import random as rd
import numpy as np
import pandas

from alinea.alep.alep_color import alep_colormap, green_yellow_red

# Imports for weather
from alinea.astk.TimeControl import *
import alinea.septo3d
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather
from alinea.alep.alep_weather import wetness_rapilly, is_raining

# Imports for wheat
from alinea.alep.wheat import initialize_stand, find_blade_id, find_leaf_ids
from alinea.astk.plant_interface import grow_canopy
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septoria_exchanging_rings import SeptoriaExchangingRings
from alinea.alep.septoria import Disease as _Disease, SeptoriaParameters as _SeptoriaParameters
from alinea.alep.septoria import plugin_septoria
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence
from alinea.alep.disease_outputs import LeafInspector
from alinea.alep.disease_outputs import compute_total_necrotic_area

# Plot ########################################################################
from alinea.alep.architecture import set_property_on_each_id
from alinea.alep.disease_outputs import compute_severity_by_leaf
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer
def update_plot(g):
    # Compute severity by leaf
    severity_by_leaf = compute_severity_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'severity', severity_by_leaf, label = 'LeafElement')

    # Visualization
    g = alep_colormap(g, 'severity', cmap=green_yellow_red(levels=100),
                      lognorm=False, zero_to_one=False, vmax=100)

    leaves = get_leaves(g, label='LeafElement')
    pos_sen = g.property('position_senescence')
    for leaf in leaves:
        if pos_sen[leaf]==0.:
            g.node(leaf).color = (157, 72, 7)

    scene = plot3d(g)
    Viewer.display(scene)
    return scene
    
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
    NewSeptoria = new_septoria_lesion()
    lesions = [NewSeptoria(nb_spores = 1.) for i in range(nb_lesions)]
    return lesions

def new_septoria(senescence_treshold=330.):
    
    class NewSeptoria(SeptoriaExchangingRings):
        def __init__(self, nb_spores=None, position=None):
            super(NewSeptoria, self).__init__(nb_spores=nb_spores, position=position)
            self.fungus.senescence_treshold = senescence_treshold
        
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
                
            # todo : modify condition
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
        
    class Parameters(_SeptoriaParameters):
        def __init__(self,**kwds):
            _SeptoriaParameters.__init__(self, model=NewSeptoria, **kwds)
            
    class Disease(_Disease):
        @classmethod
        def parameters(cls, **kwds):
            return Parameters(**kwds)
        
        @classmethod
        def lesion(cls, **kwds):
            NewSeptoria.fungus=cls.parameters(**kwds)
            return NewSeptoria
    
    return Disease
    
# Initialization #####################################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Read weather and adapt it to septoria (add wetness)
meteo_path = shared_data(alinea.septo3d, 'meteo00-01.txt')
weather = Weather(data_file=meteo_path)
weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
seq = pandas.date_range(start = "2000-10-01 01:00:00",
                        end = "2001-07-01 01:00:00", 
                        freq='H')


# Initialize a wheat canopy
g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, 
                                                width=0.2, sowing_density=150,
                                                plant_density=150, inter_row=0.12)

# Initialize the models for septoria
# septoria = new_septoria()
septoria = plugin_septoria()
inoculator = RandomInoculation()
controler = NoPriorityGrowthControl()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
emitter = SeptoriaRainEmission(domain_area=domain_area)
transporter = Septo3DSplash(reference_surface=domain_area)

# Define the schedule of calls for each model
every_h = time_filter(seq, delay=1)
every_24h = time_filter(seq, delay=24)
every_rain = rain_filter(seq, weather)
weather_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
wheat_timing = IterWithDelays(*time_control(seq, every_24h, weather.data))
septo_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))

# Call leaf inspectors for target blades (top 3)
# inspector = LeafInspector(g, blade_id=48)

# Simulation #########################################################
for i,controls in enumerate(zip(weather_timing, wheat_timing, 
                                septo_timing, rain_timing)):
    weather_eval, wheat_eval, septo_eval, rain_eval = controls
    
    # Update date
    date = weather_eval.value.index[0]
    print(date)
    
    # Get weather for date and add it as properties on leaves
    set_properties(g,label = 'LeafElement',
                   temp = weather_eval.value.temperature_air[0],
                   wetness = weather_eval.value.wetness[0],
                   rain_intensity = rain_eval.value.rain.mean(),
                   rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.,
                   relative_humidity = weather_eval.value.relative_humidity[0],
                   wind_speed = weather_eval.value.wind_speed[0])


    # Grow wheat canopy
    if wheat_eval:
        grow_canopy(g, wheat, wheat_eval.value)
    # update_healthy_area(g, label = 'LeafElement')
    
    # # Note : The position of senescence goes back to its initial value after
    # # a while for undetermined reason
    # # --> temporary hack for keeping senescence position low when it is over
    # positions = g.property('position_senescence')
    # are_green = g.property('is_green')
    # leaves = get_leaves(g, label = 'LeafElement')
    # vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
    # positions.update({vid:(0 if positions[vid]==1 and not are_green[vid] else positions[vid]) 
                      # for vid in vids})
    
    # # Develop disease
    # if data.dispersal_event.values[0]==True and timer.numiter <= 1500:
        # # Refill pool of initial inoculum to simulate differed availability
        # dispersal_units = generate_stock_du(nb_dus=10, disease=septoria)
        # initiate(g, dispersal_units, inoculator)
    # infect(g, t['disease'].dt, label='LeafElement')
    # update(g, t['disease'].dt, controler, sen_model, label='LeafElement')
    # if data.dispersal_event.values[0]==True:
        # disperse(g, emitter, transporter, "septoria", label='LeafElement')
        
    # Save outputs
    # inspector.update_variables(g)
    
    if wheat_eval:
        update_plot(g)