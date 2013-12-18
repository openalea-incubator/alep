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
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
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
def generate_lesions_with_senescence_threshold(nb_lesions=1, senescence_threshold=250.):
    """ Generate a pool of lesion with a modified senescence reponse.
    
    Senescence responds to an age threshold independent of lesion status.
    
    Parameters
    ----------
    nb_lesions: int
        Number of lesions in stock
    senescence_threshold: float
        Age at which the lesion is affected by senescence
        
    Returns
    -------
    lesions: list
        List of lesions with senescence response according to a threshold
    """
    NewSeptoria = new_septoria_lesion()
    lesions = [NewSeptoria(nb_spores = 1.) for i in range(nb_lesions)]
    return lesions

def new_septoria(senescence_threshold=330.):
    
    class NewSeptoria(SeptoriaExchangingRings):
        def __init__(self, nb_spores=None, position=None):
            super(NewSeptoria, self).__init__(nb_spores=nb_spores, position=position)
            self.fungus.senescence_threshold = senescence_threshold
        
        def senescence_response(self):
            """ Find surfaces killed by senescence.
            
            Senescence affects surfaces whose age is below the age threshold.
            """
            from math import floor
            f = self.fungus
            threshold = f.senescence_threshold
            rings = np.copy(self.surface_rings)
            surface = self.surface
            delta_ring = self.delta_age_ring
            width_ring = f.delta_age_ring
            time_to_chlo = f.degree_days_to_chlorosis
            if threshold < time_to_chlo:
                raise Exception('Do not use a threshold below time to chlorosis')
                
            # todo : modify condition
            ddday = self.ddday
            ddday_sen = self.ddday_before_senescence
            age_dday = self.age_dday
            
            # Stop growth
            self.disable_growth()
            
            if age_dday <= threshold:
                # Everything is affected
                self.disable()
                self.first_ring.disable()
                self.surface_dead = self.surface_alive
                self.surface_alive = 0.
                self.surface_rings = np.array([])
            else:
                if len(rings)>0 and self.surface_dead==0.:
                    threshold -= time_to_chlo
                    diff = delta_ring - threshold
                    nb_full_rings = floor(diff/width_ring)
                    portion_non_senescent = (diff%width_ring)/width_ring
                    if len(rings)>nb_full_rings:
                        surface_dead = sum(rings[:-(nb_full_rings+1)])+(1-portion_non_senescent)*rings[-(nb_full_rings+1)]
                        rings[:-(nb_full_rings+1)] = 0.
                        rings[-(nb_full_rings+1)] *= portion_non_senescent
                    else:
                        surface_dead = 0.
                    self.surface_rings = rings
                    self.surface_alive = sum(self.surface_rings) + self.first_ring.surface
                    self.surface_dead = surface_dead
                
                # if len(rings)>0:
                    # if nb_full_rings == 0:
                        # # The threshold is inside the current ring, which is not full
                        # portion_non_senescent = (age_dday-threshold)/(age_dday+width_ring-delta_ring)
                        # surface_non_senescent = portion_non_senescent*rings[-1]
                        # rings[-1] *= portion_non_senescent
                        # rings[:-1] = 0.
                    # else:
                        # # The threshold is inside an already full ring
                        # portion_non_senescent = (diff%width_ring)/width_ring
                        # if portion_non_senescent > 0.:
                            # surface_non_senescent = (sum(rings[-nb_full_rings:]) + 
                                                     # portion_non_senescent*rings[-(nb_full_rings+1)])
                            # rings[-(nb_full_rings+1)] *= portion_non_senescent
                            # rings[:-(nb_full_rings+1)] = 0.
                        # else:
                            # surface_non_senescent = sum(rings[-nb_full_rings:])
                            # rings[:-nb_full_rings] = 0.
                            # # update surfaces rings ...
                    # surface_non_senescent += self.first_ring.surface
                    # surface_dead = self.surface - surface_non_senescent
                # else:
                    # surface_non_senescent = self.first_ring.surface
                    # surface_dead = 0.
                
            # self.surface_alive = surface_non_senescent
            # self.surface_dead = surface_dead
            
            # Complete the age of the lesion up to the end of time step
            # self.ddday = ddday - ddday_sen
            self.ddday = ddday_sen
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
    
def run_simulation(senescence_threshold=330.):
    # Initialization #####################################################
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Read weather and adapt it to septoria (add wetness)
    meteo_path = shared_data(alinea.septo3d, 'meteo98-99.txt')
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    # seq = pandas.date_range(start = "2000-10-01 01:00:00",
                            # end = "2001-03-01 01:00:00", 
                            # freq='H')
    seq = pandas.date_range(start = "1998-10-01 01:00:00",
                            end = "1999-07-01 01:00:00", 
                            freq='H')

    # Initialize a wheat canopy
    g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, 
                                                    width=0.2, sowing_density=150,
                                                    plant_density=150, inter_row=0.12, 
                                                    seed=3)
    
    # temp
    g.node(66).id = 66
    
    
    # Initialize the models for septoria
    septoria = new_septoria(senescence_threshold=senescence_threshold)
    # septoria = plugin_septoria()
    inoculator = RandomInoculation()
    growth_controler = NoPriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
    emitter = SeptoriaRainEmission(domain_area=domain_area)
    transporter = Septo3DSplash(reference_surface=domain_area)
    washor = RapillyWashing()

    # Define the schedule of calls for each model
    every_h = time_filter(seq, delay=1)
    every_24h = time_filter(seq, delay=24)
    every_rain = rain_filter(seq, weather)
    weather_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    wheat_timing = IterWithDelays(*time_control(seq, every_24h, weather.data))
    septo_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))

    # Call leaf inspectors for target blades (top 3)
    inspectors = {}
    first_blade = 8
    ind = 0.
    for blade in range(8,104,8):
        ind += 1
        inspectors[ind] = LeafInspector(g, blade_id=blade)
    nb_lesions = []
    nb_dus = []
    from alinea.alep.disease_outputs import count_lesions

    # Simulation #########################################################
    for i,controls in enumerate(zip(weather_timing, wheat_timing, 
                                    septo_timing, rain_timing)):
        weather_eval, wheat_eval, septo_eval, rain_eval = controls
        
        # Update date
        date = weather_eval.value.index[0]

        # Get weather for date and add it as properties on leaves
        if weather_eval:
            set_properties(g,label = 'LeafElement',
                           temp = weather_eval.value.temperature_air[0],
                           wetness = weather_eval.value.wetness[0],
                           relative_humidity = weather_eval.value.relative_humidity[0],
                           wind_speed = weather_eval.value.wind_speed[0])
        if rain_eval:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = rain_eval.value.rain.mean(),
                           rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)

        # Grow wheat canopy
        if wheat_eval:
            print(date)
            g,_ = grow_canopy(g, wheat, wheat_eval.value)
            # Note : The position of senescence goes back to its initial value after
            # a while for undetermined reason
            # --> temporary hack for keeping senescence position low when it is over
            positions = g.property('position_senescence')
            are_green = g.property('is_green')
            leaves = get_leaves(g, label = 'LeafElement')
            positions.update({leaf:(0 if positions[leaf]==1 and not are_green[leaf] else positions[leaf]) 
                              for leaf in leaves})
            
        # Develop disease
        if septo_eval:
            sen_model.find_senescent_lesions(g, label = 'LeafElement')
            update_healthy_area(g, label = 'LeafElement')
            if rain_eval and i <= 500:
                # Refill pool of initial inoculum to simulate differed availability
                dispersal_units = generate_stock_du(nb_dus=rd.randint(0,3), disease=septoria)
                initiate(g, dispersal_units, inoculator)
            infect(g, septo_eval.dt, infection_controler, label='LeafElement')
            update(g, septo_eval.dt, growth_controler, sen_model, label='LeafElement')          
        if rain_eval:
            g, nb = disperse(g, emitter, transporter, "septoria", label='LeafElement')
            wash(g, washor, rain_eval.value.rain.mean(), label='LeafElement')
        else:
            nb = 0.
        nb_dus.append(nb)
        
        nb_lesions.append(count_lesions(g))
        # Save outputs
        if wheat_eval:
            for inspector in inspectors.itervalues():
                inspector.update_variables(g)
                inspector.update_du_variables(g)
        
        if wheat_eval:
            update_plot(g)
    
    return inspectors, g, nb_dus, nb_lesions

def test_update(senescence_threshold=330.):
    """ 
    """
    from alinea.alep.wheat import adel_one_leaf_element
    
    # Read weather and adapt it to septoria (add wetness)
    meteo_path = shared_data(alinea.septo3d, 'meteo00-01.txt')
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    seq = pandas.date_range(start = "2000-10-01 01:00:00",
                            end = "2001-10-31 01:00:00", 
                            freq='H')
    
    # Generate a wheat MTG
    g = adel_one_leaf_element()
    set_properties(g, label = 'LeafElement', 
                   area=5., green_area=5., healthy_area=5., position_senescence=1,
                   wetness=True, temp=22., rain_intensity=0., relative_humidity=90.)
    
    # Generate one lesion of septoria and distribute it on g
    septoria = new_septoria(senescence_threshold=senescence_threshold)
    # septoria = plugin_septoria()
    Lesion = septoria.lesion()
    leaf = g.node(10)
    leaf.lesions = [Lesion(nb_spores=1, position=[0.5, 0])]
    
    # Call model of growth control
    growth_controler = NoPriorityGrowthControl()
    sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
    
    #Temp
    emitter = SeptoriaRainEmission(domain_area=0.0004)
    transporter = Septo3DSplash(reference_surface=0.0004)
    
    # Loop of simulation
    every_h = time_filter(seq, delay=1)
    every_rain = rain_filter(seq, weather)
    septo_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    surface = []
    surface_alive = []
    surface_empty = []
    nb_dus = []
    stock_spores = []
    for i,controls in enumerate(zip(septo_timing, rain_timing)):
        septo_eval, rain_eval = controls
        
        if i==500:
            set_properties(g, label = 'LeafElement', position_senescence=0) 
        if rain_eval:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = rain_eval.value.rain.mean(),
                           rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)
            
        # Update healthy area
        sen_model.find_senescent_lesions(g, label = 'LeafElement')
        update_healthy_area(g, label = 'LeafElement')
    
        # Update
        update(g, septo_eval.dt, growth_controler, senescence_model=sen_model, label='LeafElement')
        
        if rain_eval:
            g, nb = disperse(g, emitter, transporter, "septoria", label='LeafElement')
        else:
            nb=0.
        nb_dus.append(nb)
        
        # Check that the lesion is in the right status and has the right surface
        lesion = g.property('lesions')
        if lesion:
            assert sum(len(l) for l in lesion.itervalues()) == 1
            l = lesion.values()[0][0]
            surface.append(l.surface)
            surface_alive.append(l.surface_alive)
            surface_empty.append(l.surface_empty)
            stock_spores.append(l.stock_spores)
            # if i==299:
                # import pdb
                # pdb.set_trace()
            l.compute_all_surfaces()
            f = l.fungus
            print('lesion surface: %f' % round(l.surface, 6))
            
    return g, surface, surface_alive, surface_empty, nb_dus, stock_spores