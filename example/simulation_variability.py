""" Simulate septoria epidemics on wheat with various thresholds of latency. """

# Useful imports
import random as rd
rnd = rd.Random(1)
import numpy as np
import pandas
import sys

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
from alinea.alep.septoria import plugin_septoria
from alinea.alep.septoria import Disease as _Disease, SeptoriaParameters as _SeptoriaParameters
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence
from alinea.alep.disease_outputs import LeafInspector, plot_severity_by_leaf

# Model of septoria with various latency thresholds ##########################################
def variable_septoria(mu=None, sigma=None):
    from alinea.alep.septoria_age_physio import SeptoriaAgePhysio
    
    class VariableSeptoria(SeptoriaAgePhysio):
        def __init__(self, nb_spores=None, position=None, **kwds):
            super(VariableSeptoria, self).__init__(nb_spores=nb_spores, position=position)
            self.fungus = _SeptoriaParameters(model=self.fungus.model,**self.fungus.__dict__)
            self.fungus.degree_days_to_chlorosis = rnd.gauss(mu=mu, sigma=sigma)

    class Parameters(_SeptoriaParameters):
        def __init__(self,**kwds):
            _SeptoriaParameters.__init__(self, model=VariableSeptoria, **kwds)
            
    class Disease(_Disease):
        @classmethod
        def parameters(cls, **kwds):
            return Parameters(**kwds)
        
        @classmethod
        def lesion(cls, **kwds):
            VariableSeptoria.fungus=cls.parameters(**kwds)
            return VariableSeptoria
    
    return Disease
    
# Simulation for a situation ################################################################
def run_simulation(start_year, variability=True, **kwds):
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Read weather and adapt it to septoria (add wetness)
    weather_file = 'meteo'+ str(start_year)[-2:] + '-' + str(start_year+1)[-2:] + '.txt'
    meteo_path = shared_data(alinea.septo3d, weather_file)
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    seq = pandas.date_range(start = str(start_year)+"-10-01 01:00:00",
                            end = str(start_year+1)+"-07-01 01:00:00", 
                            freq='H')

    # Initialize a wheat canopy
    g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, 
                                                    width=0.2, sowing_density=150,
                                                    plant_density=150, inter_row=0.12, 
                                                    seed=3)
    
    # Initialize the models for septoria
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
    if variability==True:
        septoria = variable_septoria(**kwds)
    else:
        septoria = plugin_septoria(model="septoria_age_physio")
    DU = septoria.dispersal_unit()
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
    first_blade = 80
    ind = 4.
    for blade in range(first_blade,104,8):
        ind -= 1
        inspectors[ind] = LeafInspector(g, blade_id=blade)
    
    # Simulation loop
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
                if rd.random()<0.4:
                    dispersal_units = [DU(nb_spores=rd.randint(1,100), status='emitted') for i in range(rd.randint(0,3))]
                    initiate(g, dispersal_units, inoculator)
            infect(g, septo_eval.dt, infection_controler, label='LeafElement')
            update(g, septo_eval.dt, growth_controler, sen_model, label='LeafElement')                
        
        les = g.property('lesions')
        lesions = sum([l for l in list(les.values())], [])
        
        print([l.fungus.degree_days_to_chlorosis for l in lesions])
        
        # if len(lesions)>10:
            # import pdb
            # pdb.set_trace()
        
        
        if rain_eval:
            g, nb = disperse(g, emitter, transporter, "septoria", label='LeafElement')
            wash(g, washor, rain_eval.value.rain.mean(), label='LeafElement')
        
        # Save outputs
        for inspector in inspectors.values():
            inspector.update_variables(g)
            inspector.update_du_variables(g)
        
        if wheat_eval:
            plot_severity_by_leaf(g)
    
    return inspectors
    
def run_all_simulations(mu=220., sigma=30.):
    # Set parameters of simulations
    start_years = [1998, 2001, 2002]

    # Run simulation for each year with variable and fixed chlorosis threshold
    for year in start_years:
        # Compute and store results for simulation with variability
        inspector = run_simulation(start_year=year, variability=True, mu=220., sigma=30.)
        stored_insp = '.\simulation_variability\inspector_'+str(year)+'_var.pckl'
        f_insp = open(stored_insp, 'w')
        pickle.dump(inspector, f_insp)
        f_insp.close()
        del inspector
        
        # Compute and store results for simulation without variability
        inspector = run_simulation(start_year=year, variability=False)
        stored_insp = '.\simulation_variability\inspector_'+str(year)+'_fix.pckl'
        f_insp = open(stored_insp, 'w')
        pickle.dump(inspector, f_insp)
        f_insp.close()
        del inspector
        
def read_outputs():
    start_years = [1998, 2001, 2002]
    out = {}
    for year in start_years:
        stored_insp =  '.\simulation_variability\inspector_'+str(year)+'_var.pckl'
        f_insp = open(stored_insp)
        out[str(year)+'_var'] = pickle.load(f_insp)
        f_insp.close()
        
        stored_insp =  '.\simulation_variability\inspector_'+str(year)+'_fix.pckl'
        f_insp = open(stored_insp)
        out[str(year)+'_fix'] = pickle.load(f_insp)
        f_insp.close()
    return out