"""

Run annual loop for the model of septoria in alep on the basis of 'annual_loop_decomposed' in echap

"""
import pandas
import sys

# Imports for wheat
from simulation_tools import wheat_path, init_canopy, grow_canopy
from alinea.echap.architectural_reconstructions import echap_reconstructions
from alinea.alep.architecture import set_properties
from alinea.adel.newmtg import adel_labels

# Imports for weather
from simulation_tools import get_weather
from alinea.alep.alep_time_control import CustomIterWithDelays, septoria_filter
from alinea.astk.TimeControl import (IterWithDelays, rain_filter, time_filter,
                                     thermal_time_filter, DegreeDayModel,
                                     time_control)

# Imports for alep septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum
from alinea.popdrops.alep_interface import PopDropsSoilContamination, PopDropsEmission, PopDropsTransport
from alinea.alep.growth_control import PriorityGrowthControl,SeptoRustCompetition, GeometricPoissonCompetition
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.disease_outputs import plot_severity_by_leaf, save_image, AdelSeptoRecorder
from variable_septoria import *

def setup(sowing_date="2010-10-15 12:00:00", start_date = None,
          end_date="2011-06-20 01:00:00", variety='Mercia',
          nplants = 30, nsect = 7, disc_level = 5, Tmin = 10., Tmax = 25., WDmin = 10., 
          rain_min = 0.2, recording_delay = 24., reset_reconst = True):
    """ Get plant model, weather data and set scheduler for simulation. """
    # Set canopy
    it_wheat = 0
    reconst = echap_reconstructions(reset=True, reset_data=True)
    adel = reconst.get_reconstruction(name=variety, nplants=nplants, nsect=nsect)
    year = int(end_date[:4])    
    wheat_dir = wheat_path((year, variety, nplants, nsect))
    g, wheat_is_loaded = init_canopy(adel, wheat_dir, rain_and_light=True) 
            
    # Manage weather
    weather = get_weather(start_date = sowing_date, end_date = end_date)
    
    # Define the schedule of calls for each model
    if start_date is None:
        start_date = sowing_date
    seq = pandas.date_range(start = start_date, end = end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase = 0)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 20.)
    every_rain = rain_filter(seq, weather, rain_min = rain_min)
    every_recording = time_filter(seq, delay=recording_delay)
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')
    septo_filter = septoria_filter(seq, weather, degree_days = 10., Tmin = Tmin, Tmax = Tmax, WDmin = WDmin, rain_min = rain_min)
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    septo_timing = CustomIterWithDelays(*time_control(seq, septo_filter, weather.data), eval_time='end')
    recorder_timing = IterWithDelays(*time_control(seq, every_recording, weather.data))
    return (g, adel, weather, seq, rain_timing, canopy_timing, septo_timing, 
            recorder_timing, it_wheat, wheat_dir, wheat_is_loaded)

def septo_disease(adel, sporulating_fraction, layer_thickness, distri_chlorosis = None, **kwds):
    """ Choose models to assemble the disease model. """
               
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
            
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
    inoculum = SoilInoculum(fungus=fungus, 
                            sporulating_fraction=sporulating_fraction,
                            domain_area=domain_area, mutable=mutable)
    contaminator = PopDropsSoilContamination(fungus=fungus,
                                             group_dus=True,
                                             domain=domain, domain_area=domain_area, 
                                             compute_star = False)
#    growth_controler = PriorityGrowthControl()
#    growth_controler = GeometricPoissonCompetition()
    growth_controler = SeptoRustCompetition()
    infection_controler = BiotrophDUProbaModel()
    emitter = PopDropsEmission(domain=domain, compute_star = False)
    transporter = PopDropsTransport(fungus=fungus, group_dus=True,
                                    domain = domain, domain_area = domain_area,
                                    dh = layer_thickness, convUnit = convUnit,
                                    compute_star = False)
    return inoculum, contaminator, infection_controler, growth_controler, emitter, transporter

def annual_loop_septo(year = 2012, variety = 'Tremie13', sowing_date = '10-15',
                      nplants = 30, nsect = 7,
                      sporulating_fraction = 1e-4, layer_thickness = 0.01, 
                      record = True, output_file = None,
                      save_images = False, reset_reconst = True, 
                      distri_chlorosis = None, **kwds):
    """ Simulate epidemics with canopy saved before simulation """
    if 'temp_min' in kwds:
        Tmin = kwds['temp_min']
    else:
        Tmin = 0.
    (g, adel, weather, seq, rain_timing, 
     canopy_timing, septo_timing, recorder_timing, it_wheat, wheat_dir,
     wheat_is_loaded) = setup(sowing_date=str(year-1)+"-"+sowing_date+" 12:00:00", 
                              end_date=str(year)+"-07-01 00:00:00",
                              variety = variety, nplants = nplants,
                              nsect = nsect, Tmin = Tmin, 
                              reset_reconst = reset_reconst)

    (inoculum, contaminator, infection_controler, growth_controler, emitter, 
     transporter) = septo_disease(adel, sporulating_fraction, layer_thickness, 
                                    distri_chlorosis, **kwds)

    a_labels = adel_labels(g)
    leaves = [k for k,v in a_labels.iteritems() if v.startswith('plant1_MS_metamer1_blade_LeafElement')]
    for lf in leaves:
        g.node(lf).tag = 'tag'

    # Prepare saving of outputs
    if record == True:
        recorder = AdelSeptoRecorder()
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, 
                                     septo_timing, recorder_timing)):
        canopy_iter, rain_iter, septo_iter, record_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            it_wheat += 1
            g = grow_canopy(g, adel, canopy_iter, it_wheat,
                        wheat_dir, wheat_is_loaded,rain_and_light=True)
                        
            a_labels = adel_labels(g)
            leaves = [k for k,v in a_labels.iteritems() if v.startswith('plant1_MS_metamer1_blade_LeafElement')]
            for lf in leaves:
                g.node(lf).tag = 'tag'
                
                
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
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.2:
            g = external_contamination(g, inoculum, contaminator, rain_iter.value,
                                       domain=adel.domain, 
                                       domain_area=adel.domain_area)
        # Develop disease (infect for dispersal units and update for lesions)
        if septo_iter:
            infect(g, septo_iter.dt, infection_controler, label='LeafElement')
            update(g, septo_iter.dt, growth_controler, senescence_model=None, label='LeafElement')            
        # Disperse and wash
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.2:
            g = disperse(g, emitter, transporter, "septoria",
                         label='LeafElement', weather_data=rain_iter.value,
                         domain=adel.domain, domain_area=adel.domain_area)
        # Save images
        if save_images == True:
            if canopy_iter:
                scene = plot_severity_by_leaf(g)
                if it_wheat < 10 :
                    image_name='./images_septo/image0000%d.png' % it_wheat
                elif it_wheat < 100 :
                    image_name='./images_septo/image000%d.png' % it_wheat
                elif it_wheat < 1000 :
                    image_name='./images_septo/image00%d.png' % it_wheat
                elif it_wheat < 10000 :
                    image_name='./images_septo/image0%d.png' % it_wheat
                else :
                    image_name='./images_segpto/image%d.png' % it_wheat
                save_image(scene, image_name=image_name)
        # Save outputs
        if record_iter and record == True:
            date = record_iter.value.index[-1]
            print date
            recorder.record(g, date, degree_days = record_iter.value.degree_days[-1])
                    
    if record == True:
        recorder.post_treatment(variety = 'tremie')
        if output_file is not None:
            recorder.save(output_file)
        else:
            return g, recorder
    else:
        return g
    
def stat_profiler(call='run_disease()'):
    import cProfile
    import pstats
    cProfile.run(call, 'septo_stats')
    return pstats.Stats('septo_stats')