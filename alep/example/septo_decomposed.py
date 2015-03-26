"""

Run annual loop for the model of septoria in alep on the basis of 'annual_loop_decomposed' in echap

"""
import pandas
import sys
import numpy as np

# Imports for wheat and rain/light interception
from alinea.adel.newmtg import move_properties, adel_ids
from alinea.echap.architectural_reconstructions import (EchapReconstructions,
                                                        get_EchapReconstructions)
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves
from alinea.caribu.caribu_star import rain_and_light_star

# Imports for weather
from alinea.alep.alep_weather import wetness_rapilly, linear_degree_days
from alinea.alep.alep_time_control import *
from alinea.astk.TimeControl import *
from alinea.echap.weather_data import *
import alinea.septo3d
import alinea.alep
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather

# Imports for alep septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum, Septo3DSoilContamination, Septo3DEmission
from alinea.popdrops.alep_interface import PopDropsSoilContamination, PopDropsEmission, PopDropsTransport
from alinea.alep.growth_control import PriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.disease_outputs import plot_severity_by_leaf, save_image, AdelSeptoRecorder
from variable_septoria import *

def get_weather(start_date="2010-10-15 12:00:00", end_date="2011-06-20 01:00:00"):
    """ Get weather data for simulation. """
    start = datetime.datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    if start.year >= 2010:
        filename = 'Boigneville_0109'+str(start.year)+'_3108'+str(start.year+1)+'_h.csv'
        meteo_path = shared_data(alinea.echap, filename)
        weather = Weather(meteo_path, reader = arvalis_reader)
        weather.check(['temperature_air', 'PPFD', 'relative_humidity',
                       'wind_speed', 'rain', 'global_radiation', 'vapor_pressure'])
        notation_dates_file = shared_data(alinea.alep, 'notation_dates/notation_dates_'+str(start.year+1)+'.csv')
        weather.check(varnames=['notation_dates'], models={'notation_dates':add_notation_dates}, notation_dates_file = notation_dates_file)
    else:
        start_yr = str(start.year)[2:4]
        end_yr = str(start.year+1)[2:4]
        filename = 'meteo'+ start_yr + '-' + end_yr + '.txt'
        meteo_path = shared_data(alinea.septo3d, filename)
        weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    weather.check(varnames=['degree_days'], models={'degree_days':linear_degree_days}, start_date=start_date, base_temp=0., max_temp=30.)
    # weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':linear_degree_days}, start_date=start_date, base_temp=0., max_temp=25.)
    weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':linear_degree_days}, start_date=start_date, base_temp=0., max_temp=30.)
    return weather

def setup(start_date="2010-10-15 12:00:00", end_date="2011-06-20 01:00:00", variety='Mercia',
            nplants = 30, nsect = 7, disc_level = 5, Tmin = 10., Tmax = 25., WDmin = 10., 
            rain_min = 0.2, reset_reconst = False):
    """ Get plant model, weather data and set scheduler for simulation. """
    # Initialize wheat plant
    if reset_reconst == False:
        reconst = get_EchapReconstructions()
    else:
        reconst = EchapReconstructions()
    # ULTRA TEMP ET MOCHE 
    adel = None
    while adel is None:
        try:
            adel = reconst.get_reconstruction(name=variety, nplants = nplants,
                nsect = nsect, disc_level = disc_level, aspect = 'line', 
                seed = np.random.random_integers(100))
        except:
            pass
            
    # Manage weather
    weather = get_weather(start_date = start_date, end_date = end_date)
    
    # Define the schedule of calls for each model
    seq = pandas.date_range(start = start_date, end = end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase = 0)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 20.)
    every_rain = rain_filter(seq, weather, rain_min = rain_min)
    # every_dd_or_rain = filter_or([every_dd, every_rain])
    # canopy_timing = IterWithDelays(*time_control(seq, every_dd_or_rain, weather.data))
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')
    septo_filter = septoria_filter(seq, weather, degree_days = 10., Tmin = Tmin, Tmax = Tmax, WDmin = WDmin, rain_min = rain_min)
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    septo_timing = CustomIterWithDelays(*time_control(seq, septo_filter, weather.data), eval_time='end')
    return adel, weather, seq, rain_timing, canopy_timing, septo_timing

def septo_disease(adel, sporulating_fraction, layer_thickness, distri_chlorosis = None, **kwds):
    """ Choose models to assemble the disease model. """
    domain = adel.domain
    domain_area = adel.domain_area
    convUnit = adel.convUnit
    if distri_chlorosis is not None:
        fungus = variable_septoria(distri_chlorosis = distri_chlorosis)
        mutable = True
    else:
        fungus = plugin_septoria()
        mutable = False
    fungus.parameters(group_dus=True, nb_rings_by_state=1, **kwds)
    inoculum = SoilInoculum(fungus, sporulating_fraction=sporulating_fraction,
                            domain_area=domain_area, mutable = mutable)
    contaminator = PopDropsSoilContamination(domain=domain, domain_area=domain_area)
    growth_controler = PriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    emitter = PopDropsEmission(domain=domain, compute_star = False)
    transporter = PopDropsTransport(domain = domain, domain_area = domain_area, dh = layer_thickness, convUnit = convUnit)
    return inoculum, contaminator, infection_controler, growth_controler, emitter, transporter
   
def make_canopy(start_date = "2010-10-15 12:00:00", end_date = "2011-06-20 01:00:00", variety = 'Mercia',
                nplants = 30, nsect = 7, disc_level = 5, dir = './adel/mercia_2011_30pl_7sect', 
                reset_reconst = False):
    """ Simulate and save canopy (prior to simulation). """
    adel, weather, seq, rain_timing, canopy_timing, septo_timing = setup(start_date = start_date,
        end_date = end_date, variety = variety, nplants = nplants, nsect = nsect, disc_level = disc_level,
        reset_reconst = reset_reconst)
    
    domain = adel.domain
    convUnit = adel.convUnit
    g = adel.setup_canopy(age=0.)
    rain_and_light_star(g, light_sectors = '1', domain = domain, convUnit = convUnit)
    it_wheat = 0
    adel.save(g, it_wheat, dir=dir)
    for i, controls in enumerate(zip(canopy_timing, septo_timing)):
        canopy_iter, septo_iter = controls
        if canopy_iter:
            it_wheat += 1
            g = adel.grow(g, canopy_iter.value)
            rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            adel.save(g, it_wheat, dir=dir)

def run_disease(start_date = "2010-10-15 12:00:00", end_date = "2011-06-20 01:00:00", 
                variety = 'Mercia', nplants = 30, nsect = 7,
                disc_level = 5, dir = './adel/mercia_2011_30pl_7sect', 
                sporulating_fraction = 1e-4, layer_thickness = 0.01, record = True, 
                save_images = False,
                adel = None, weather = None, seq = None, rain_timing = None, 
                canopy_timing = None, septo_timing = None, distri_chlorosis = None, **kwds):
    """ Simulate epidemics with canopy saved before simulation """
    if any(x==None for x in [adel, weather, seq, rain_timing, canopy_timing, septo_timing]):
        if 'temp_min' in kwds:
            Tmin = kwds['temp_min']
        else:
            Tmin = 0.
        adel, weather, seq, rain_timing, canopy_timing, septo_timing = setup(start_date = start_date,
                            end_date = end_date, variety = variety, nplants = nplants, nsect = nsect, 
                            disc_level = disc_level, Tmin = Tmin)
            
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
        
    inoculum, contaminator, infection_controler, growth_controler, emitter, transporter = septo_disease(adel, 
                                    sporulating_fraction, layer_thickness, distri_chlorosis, **kwds)
    it_wheat = 0
    g,TT = adel.load(it_wheat, dir=dir)
    leaf_ids = adel_ids(g)
        
    # Prepare saving of outputs
    if record == True:
        recorder = AdelSeptoRecorder()
    else:
        recorder = None
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, septo_timing)):
        canopy_iter, rain_iter, septo_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            # print canopy_iter.value.index[-1]
            it_wheat += 1
            newg,TT = adel.load(it_wheat, dir=dir)
            move_properties(g, newg)
            g = newg
            leaf_ids = adel_ids(g)
        
        # Get weather for date and add it as properties on leaves
        if septo_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = septo_iter.value.temperature_air,
                           wetness_sequence = septo_iter.value.wetness,
                           dd_sequence = septo_iter.value.degree_days)
        if rain_iter:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = rain_iter.value.rain.mean(),
                           rain_duration = len(rain_iter.value.rain) if rain_iter.value.rain.sum() > 0 else 0.)

        # External contamination
        geom = g.property('geometry')
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.:
            g = external_contamination(g, inoculum, contaminator, rain_iter.value)

        # Develop disease (infect for dispersal units and update for lesions)
        if septo_iter:
            infect(g, septo_iter.dt, infection_controler, label='LeafElement')
            update(g, septo_iter.dt, growth_controler, senescence_model=None, label='LeafElement')
            
        # Disperse and wash
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.:
            g = disperse(g, emitter, transporter, "septoria", label='LeafElement', weather_data=rain_iter.value)

        # if save_images == True:
            # if canopy_iter:
                # scene = plot_severity_by_leaf(g)
                # if it_wheat < 10 :
                    # image_name='./images_septo/image0000%d.png' % it_wheat
                # elif it_wheat < 100 :
                    # image_name='./images_septo/image000%d.png' % it_wheat
                # elif it_wheat < 1000 :
                    # image_name='./images_septo/image00%d.png' % it_wheat
                # elif it_wheat < 10000 :
                    # image_name='./images_septo/image0%d.png' % it_wheat
                # else :
                    # image_name='./images_septo/image%d.png' % it_wheat
                # save_image(scene, image_name=image_name)
            
        # Save outputs
        if septo_iter and record == True:     
            date = septo_iter.value.index[-1]
            recorder.record(g, date, degree_days = septo_iter.value.degree_days[-1])
                    
    if record == True:
        recorder.post_treatment(variety = 'tremie')
    
    return g, recorder
    
def run_disease_and_canopy(start_date = "2010-10-15 12:00:00", end_date = "2011-06-20 01:00:00", 
                            variety = 'Mercia', nplants = 30, nsect = 7, disc_level = 5, 
                            sporulating_fraction = 1e-4, layer_thickness = 0.01, record = True, 
                            save_images = False, adel = None, weather = None, seq = None, 
                            rain_timing = None, canopy_timing = None, septo_timing = None, 
                            distri_chlorosis = None, **kwds):
    """ Simulate epidemics with canopy simulated during simulation """
    if any(x==None for x in [adel, weather, seq, rain_timing, canopy_timing, septo_timing]):
        if 'temp_min' in kwds:
            Tmin = kwds['temp_min']
        else:
            Tmin = 0.
        adel, weather, seq, rain_timing, canopy_timing, septo_timing = setup(start_date = start_date,
                            end_date = end_date, variety = variety, nplants = nplants, nsect = nsect, 
                            disc_level = disc_level, Tmin = Tmin)
            
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
        
    inoculum, contaminator, infection_controler, growth_controler, emitter, transporter = septo_disease(adel, 
                                    sporulating_fraction, layer_thickness, distri_chlorosis, **kwds)
    domain = adel.domain
    convUnit = adel.convUnit
    g = adel.setup_canopy(age=0.)
    rain_and_light_star(g, light_sectors = '1', domain = domain, convUnit = convUnit)
    leaf_ids = adel_ids(g)
        
    # Prepare saving of outputs
    if record == True:
        recorder = AdelSeptoRecorder()
    else:
        recorder = None
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, septo_timing)):
        canopy_iter, rain_iter, septo_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            # print canopy_iter.value.index[-1]
            g = adel.grow(g, canopy_iter.value)
            rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            leaf_ids = adel_ids(g)
        
        # Get weather for date and add it as properties on leaves
        if septo_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = septo_iter.value.temperature_air,
                           wetness_sequence = septo_iter.value.wetness,
                           dd_sequence = septo_iter.value.degree_days)
        if rain_iter:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = rain_iter.value.rain.mean(),
                           rain_duration = len(rain_iter.value.rain) if rain_iter.value.rain.sum() > 0 else 0.)

        # External contamination
        geom = g.property('geometry')
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.:
            g = external_contamination(g, inoculum, contaminator, rain_iter.value)

        # Develop disease (infect for dispersal units and update for lesions)
        if septo_iter:
            infect(g, septo_iter.dt, infection_controler, label='LeafElement')
            update(g, septo_iter.dt, growth_controler, senescence_model=None, label='LeafElement')
            
        # Disperse and wash
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.:
            g = disperse(g, emitter, transporter, "septoria", label='LeafElement', weather_data=rain_iter.value)

        # if save_images == True:
            # if canopy_iter:
                # scene = plot_severity_by_leaf(g)
                # if it_wheat < 10 :
                    # image_name='./images_septo/image0000%d.png' % it_wheat
                # elif it_wheat < 100 :
                    # image_name='./images_septo/image000%d.png' % it_wheat
                # elif it_wheat < 1000 :
                    # image_name='./images_septo/image00%d.png' % it_wheat
                # elif it_wheat < 10000 :
                    # image_name='./images_septo/image0%d.png' % it_wheat
                # else :
                    # image_name='./images_septo/image%d.png' % it_wheat
                # save_image(scene, image_name=image_name)
            
        # Save outputs
        if septo_iter and record == True:     
            date = septo_iter.value.index[-1]
            recorder.record(g, date, degree_days = septo_iter.value.degree_days[-1])
                    
    if record == True:
        recorder.post_treatment(variety = 'tremie')
    
    return g, recorder
    
def stat_profiler(call='run_disease()'):
    import cProfile
    import pstats
    cProfile.run(call, 'restats')
    return pstats.Stats('restats')