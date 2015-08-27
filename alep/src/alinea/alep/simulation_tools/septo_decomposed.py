"""

Run annual loop for the model of septoria in alep on the basis of 'annual_loop_decomposed' in echap

"""
import pandas as pd
import sys

# Imports for wheat
from alinea.alep.simulation_tools.simulation_tools import (wheat_path, 
                                                           init_canopy, 
                                                           grow_canopy,
                                                           alep_echap_reconstructions,
                                                           get_iter_rep_wheats,
                                                           get_filename)
from alinea.alep.architecture import set_properties

# Imports for weather
from simulation_tools import get_weather
from alinea.alep.alep_time_control import CustomIterWithDelays, septoria_filter_ddays
from alinea.astk.TimeControl import (IterWithDelays, rain_filter, time_filter,
                                     thermal_time_filter, DegreeDayModel,
                                     time_control)

# Imports for alep septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.alep.simulation_tools.simulation_tools import group_duplicates_in_cohort
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum
from alinea.popdrops.alep_interface import PopDropsSoilContamination, PopDropsEmission, PopDropsTransport
from alinea.alep.growth_control import PriorityGrowthControl,SeptoRustCompetition, GeometricPoissonCompetition
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.disease_outputs import save_image, AdelSeptoRecorder
from variable_septoria import *

# Temp
import alinea.alep
from openalea.deploy.shared_data import shared_data

from alinea.alep.disease_outputs import plot_by_leaf
from alinea.alep.simulation_tools.simulation_tools import add_leaf_dates_to_data

def setup(sowing_date="2010-10-15 12:00:00", start_date = None,
          end_date="2011-06-20 01:00:00", variety='Mercia',
          nplants = 30, nsect = 7, disc_level = 5, septo_delay_dday = 10.,
          rain_min = 0.2, recording_delay = 24., rep_wheat = None):
    """ Get plant model, weather data and set scheduler for simulation. """
    # Set canopy
    it_wheat = 0
    reconst = alep_echap_reconstructions()
    adel = reconst.get_reconstruction(name=variety, nplants=nplants, nsect=nsect)
    year = int(end_date[:4])
    wheat_dir = wheat_path(year, variety, nplants, nsect, rep_wheat)
    g, wheat_is_loaded = init_canopy(adel, wheat_dir, rain_and_light=True) 
            
    # Manage weather
    weather = get_weather(start_date = sowing_date, end_date = end_date)
    
    # Define the schedule of calls for each model
    if start_date is None:
        start_date = sowing_date
    seq = pd.date_range(start = start_date, end = end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase = 0)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 20.)
    every_rain = rain_filter(seq, weather, rain_min = rain_min)
    every_recording = time_filter(seq, delay=recording_delay)
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')
    septo_filter = septoria_filter_ddays(seq, weather, delay = septo_delay_dday, rain_min = rain_min)
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    septo_timing = CustomIterWithDelays(*time_control(seq, septo_filter, weather.data), eval_time='end')
    recorder_timing = IterWithDelays(*time_control(seq, every_recording, weather.data))
    return (g, adel, weather, seq, rain_timing, canopy_timing, septo_timing, 
            recorder_timing, it_wheat, wheat_dir, wheat_is_loaded)

def septo_disease(adel, sporulating_fraction, layer_thickness,
                  distri_chlorosis = None, **kwds):
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
                            domain_area=domain_area)
    contaminator = PopDropsSoilContamination(fungus=fungus,
                                             group_dus=True,
                                             domain=domain, domain_area=domain_area, 
                                             compute_star = False,
                                             mutable = mutable)
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

def annual_loop_septo(year = 2013, variety = 'Tremie13', sowing_date = '10-29',
                      nplants = 15, nsect = 7, septo_delay_dday = 10.,
                      sporulating_fraction = 1e-4, layer_thickness = 0.01, 
                      record = True, output_file = None,
                      save_images = False, reset_reconst = True, 
                      distri_chlorosis = None, rep_wheat = None, **kwds):
    """ Simulate epidemics with canopy saved before simulation """
    (g, adel, weather, seq, rain_timing, 
     canopy_timing, septo_timing, recorder_timing, it_wheat, wheat_dir,
     wheat_is_loaded) = setup(sowing_date=str(year-1)+"-"+sowing_date+" 12:00:00", 
                              end_date=str(year)+"-08-01 00:00:00",
                              variety = variety, nplants = nplants,
                              nsect=nsect, rep_wheat=rep_wheat, 
                              septo_delay_dday=septo_delay_dday)

    (inoculum, contaminator, infection_controler, growth_controler, emitter, 
     transporter) = septo_disease(adel, sporulating_fraction, layer_thickness, 
                                    distri_chlorosis, **kwds)
    
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
                        wheat_dir, wheat_is_loaded, rain_and_light=True)               
                
        # Get weather for date and add it as properties on leaves
        if septo_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = septo_iter.value.temperature_air,
                           wetness_sequence = septo_iter.value.wetness,
                           relative_humidity_sequence = septo_iter.value.relative_humidity,
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
            group_duplicates_in_cohort(g) # Additional optimisation (group identical cohorts)
            update(g, septo_iter.dt, growth_controler, senescence_model=None, label='LeafElement')            
        # Disperse and wash
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.2:
            g = disperse(g, emitter, transporter, "septoria",
                         label='LeafElement', weather_data=rain_iter.value,
                         domain=adel.domain, domain_area=adel.domain_area)
        # Save images
        if save_images == True:
            if canopy_iter:
                scene = plot_severity_septo_by_leaf(g)
                if it_wheat < 10 :
                    image_name=variety+'_image0000%d.png' % it_wheat
                elif it_wheat < 100 :
                    image_name=variety+'_image000%d.png' % it_wheat
                elif it_wheat < 1000 :
                    image_name=variety+'_image00%d.png' % it_wheat
                elif it_wheat < 10000 :
                    image_name=variety+'_image0%d.png' % it_wheat
                else :
                    image_name='image%d.png' % it_wheat
                image_name = str(shared_data(alinea.alep)/'images_septo'/image_name)
                save_image(scene, image_name=image_name)
        # Save outputs
        if record_iter and record == True:
            date = record_iter.value.index[-1]
            print date
            recorder.record(g, date, degree_days = record_iter.value.degree_days[-1])
                    
    if record == True:
        recorder.post_treatment(variety = variety)
        if output_file is not None:
            recorder.save(output_file)
        else:
            return g, recorder
    else:
        return g
        
def run_reps_septo(year = 2013, variety = 'Tremie13', 
                   nplants = 15, nsect = 7, sowing_date = '10-15',
                   sporulating_fraction = 5e-3, layer_thickness=0.01, 
                   nreps = 10, suffix = None, **kwds):
    df = pd.DataFrame()
    rep_wheats = get_iter_rep_wheats(year, variety, nplants, nsect, nreps)
    for rep in range(nreps):
        g, recorder = annual_loop_septo(year=year, variety=variety,
                                            sowing_date=sowing_date,
                                            nplants=nplants, nsect=nsect,
                                            sporulating_fraction=sporulating_fraction,
                                            layer_thickness=layer_thickness, 
                                            rep_wheat=next(rep_wheats), **kwds)
        df_ = recorder.data
        df_['rep'] = rep
        df = pd.concat([df, df_])
    output_file = get_filename(fungus='septoria', year=year, variety=variety,
                               nplants=nplants, inoc=sporulating_fraction,
                               suffix=suffix)
    df.to_csv(output_file)
    
def stat_profiler(call='annual_loop_septo()'):
    import cProfile
    import pstats
    cProfile.run(call, 'septo_stats')
    return pstats.Stats('septo_stats')
    
def plot_severity_septo_by_leaf(g, senescence=True,
                                transparency=None, 
                                label='LeafElement'):
    """ Display the MTG with colored leaves according to disease severity 
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    senescence: bool
        True if senescence must be displayed, False otherwise
    transparency: float[0:1]
        Transparency of the part of the MTG without lesion
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    scene:
        Scene containing the MTG attacked by the disease
    
    """
    from alinea.alep.architecture import set_property_on_each_id, get_leaves
    from alinea.alep.alep_color import alep_colormap, green_yellow_red
    from alinea.adel.mtg_interpreter import plot3d
    from alinea.alep.disease_outputs import plot3d_transparency
    from openalea.plantgl.all import Viewer
    # Visualization
    lesions = g.property('lesions')
    leaves = get_leaves(g, label=label)
    severity_by_leaf = {}
    for lf in leaves:
        if lf in lesions:
            leaf = g.node(lf)
            severity_by_leaf[lf] = sum([l.surface_spo+l.surface_empty for l in leaf.lesions])*100./leaf.area if leaf.area>0 else 0.
        else:
            severity_by_leaf[lf] = 0.
    set_property_on_each_id(g, 'severity', severity_by_leaf, label=label)
    g = alep_colormap(g, 'severity', cmap=green_yellow_red(levels=100),
                      lognorm=False, zero_to_one=False, vmax=100)

    if senescence==True:
        leaves = get_leaves(g, label=label)
        sen_lengths = g.property('senesced_length')
        green_lengths = g.property('green_length')
        for leaf in leaves:
            if sen_lengths[leaf]>0. and round(green_lengths[leaf],15)==0.:
                g.node(leaf).color = (157, 72, 7)
    
    if transparency!=None:
        for id in g:           
            if not id in severity_by_leaf:
                g.node(id).color = (255,255,255)
                g.node(id).transparency = 0.9
            elif severity_by_leaf[id]==0.:
                g.node(id).color = (255,255,255)
                g.node(id).transparency = transparency
            else:
                g.node(id).transparency = 0.
        scene = plot3d_transparency(g)
    else:
        scene = plot3d_transparency(g)
    Viewer.display(scene)
    return scene
    
# TEMP
def temp1():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=1.,
                   growth_rate=0.0006, density_dus_emitted_ref=1.79e3,
                   nreps=5, suffix='ref')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=1e-2,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=1.,
                   growth_rate=0.0006, density_dus_emitted_ref=1.79e3, 
                   nreps=5, suffix='inoc')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.03, proba_inf=1.,
                   growth_rate=0.0006, density_dus_emitted_ref=1.79e3, 
                   nreps=5, suffix='smin')
                   
def temp2():
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=0.5,
                   growth_rate=0.0006, density_dus_emitted_ref=1.79e3,
                   nreps=5, suffix='proba')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=1.,
                   growth_rate=0.001, density_dus_emitted_ref=1.79e3, 
                   nreps=5, suffix='rate')
    run_reps_septo(year=2013, variety='Tremie13', sowing_date='10-29',
                   nplants=15, sporulating_fraction=5e-3,
                   degree_days_to_chlorosis=120., Smin=0.01, proba_inf=1,
                   growth_rate=0.0006, density_dus_emitted_ref=3e3, 
                   nreps=5, suffix='emission')
                   