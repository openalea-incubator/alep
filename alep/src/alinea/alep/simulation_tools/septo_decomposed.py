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
                                                           alep_custom_reconstructions,
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
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum, Septo3DEmission
from alinea.popdrops.alep_interface import PopDropsSoilContamination, PopDropsEmission, PopDropsTransport
from alinea.alep.growth_control import PriorityGrowthControl,SeptoRustCompetition, GeometricPoissonCompetition
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.disease_outputs import save_image, AdelSeptoRecorder
from variable_septoria import *

# Temp
import alinea.alep
from openalea.deploy.shared_data import shared_data
from alinea.echap.disease.alep_septo_evaluation import *
from alinea.alep.disease_outputs import plot_by_leaf
from alinea.alep.simulation_tools.simulation_tools import add_leaf_dates_to_data
from alinea.adel.newmtg import adel_labels

def setup(sowing_date="2010-10-15 12:00:00", start_date = None,
          end_date="2011-06-20 01:00:00", variety='Mercia',
          nplants = 30, nsect = 7, disc_level = 5, septo_delay_dday = 10.,
          rain_min = 0.2, recording_delay = 24., rep_wheat = None,
          save_images=False, keep_leaves=False, leaf_duration=2., **kwds):
    """ Get plant model, weather data and set scheduler for simulation. """
    # Set canopy
    it_wheat = 0
    if variety!='Custom':
        reconst = alep_echap_reconstructions(keep_leaves=keep_leaves, leaf_duration=leaf_duration)
        adel = reconst.get_reconstruction(name=variety, nplants=nplants, nsect=nsect)
        if save_images:
            adel.stand.density_curve=None
    else:
        adel = alep_custom_reconstructions(variety='Tremie13', nplants=nplants, nsect=nsect, **kwds)
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
                  distri_chlorosis = None, competition='poisson',
                  age_infection=False, compute_star=False, **kwds):
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
                                             compute_star = compute_star,
                                             mutable = mutable)
#    growth_controler = PriorityGrowthControl()
#    growth_controler = GeometricPoissonCompetition()
    if competition=='poisson':
        growth_controler = SeptoRustCompetition()
    elif competition=='simple':
        growth_controler = PriorityGrowthControl()
    else:
        raise ValueError('Unknown competition model')
    infection_controler = BiotrophDUProbaModel(age_infection=age_infection)
#    emitter = Septo3DEmission(domain=domain)
    emitter = PopDropsEmission(domain=domain, compute_star = compute_star)
    transporter = PopDropsTransport(fungus=fungus, group_dus=True,
                                    domain = domain, domain_area = domain_area,
                                    dh = layer_thickness, convUnit = convUnit,
                                    compute_star = compute_star, wash=True)
    return inoculum, contaminator, infection_controler, growth_controler, emitter, transporter

def annual_loop_septo(year = 2013, variety = 'Tremie13', sowing_date = '10-29',
                      nplants = 15, nsect = 7, septo_delay_dday = 10.,
                      sporulating_fraction = 5e-3, layer_thickness = 0.01, 
                      record = True, output_file = None, 
                      competition = 'poisson', save_images = False, 
                      reset_reconst = True, distri_chlorosis = None, 
                      rep_wheat = None, age_infection=False, keep_leaves=False,
                      leaf_duration = 2., compute_star = False, **kwds):
    """ Simulate epidemics with canopy saved before simulation """
    (g, adel, weather, seq, rain_timing, 
     canopy_timing, septo_timing, recorder_timing, it_wheat, wheat_dir,
     wheat_is_loaded) = setup(sowing_date=str(year-1)+"-"+sowing_date+" 12:00:00", 
                              end_date=str(year)+"-07-30 00:00:00",
                              variety = variety, nplants = nplants,
                              nsect=nsect, rep_wheat=rep_wheat, 
                              septo_delay_dday=septo_delay_dday,
                              save_images=save_images, 
                              keep_leaves=keep_leaves, 
                              leaf_duration=leaf_duration,**kwds)
#    (g, adel, weather, seq, rain_timing, 
#     canopy_timing, septo_timing, recorder_timing, it_wheat, wheat_dir,
#     wheat_is_loaded) = setup(sowing_date=str(year-1)+"-"+sowing_date+" 12:00:00", 
#                              end_date=str(year)+"-07-30 00:00:00",
#                              variety = variety, nplants = nplants,
#                              nsect=nsect, rep_wheat=rep_wheat, 
#                              septo_delay_dday=septo_delay_dday,
#                              save_images=save_images, 
#                              keep_leaves=keep_leaves, 
#                              leaf_duration=leaf_duration, **kwds)

    (inoculum, contaminator, infection_controler, growth_controler, emitter, 
     transporter) = septo_disease(adel, sporulating_fraction, layer_thickness, 
                                    distri_chlorosis, competition=competition,
                                    age_infection=age_infection, 
                                    compute_star=compute_star, 
                                    **kwds)
    
    # Prepare saving of outputs
    if record == True:
        recorder = AdelSeptoRecorder(add_height=True)
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, 
                                     septo_timing, recorder_timing)):
        canopy_iter, rain_iter, septo_iter, record_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            it_wheat += 1
            g = grow_canopy(g, adel, canopy_iter, it_wheat,
                        wheat_dir, wheat_is_loaded, rain_and_light=True)
                
        # TEMP: debug
        g.add_property('a_label')
        a_labels = g.property('a_label')
        a_labels.update( adel_labels(g))
                        
                
        # Get weather for date and add it as properties on leaves
        if septo_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = septo_iter.value.temperature_air.tolist(),
                           wetness_sequence = septo_iter.value.wetness.tolist(),
                           relative_humidity_sequence = septo_iter.value.relative_humidity.tolist(),
                           dd_sequence = septo_iter.value.degree_days.tolist())
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
                scene = plot_severity_septo_by_leaf(g, senescence=False)
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
    
                   
def set_canopy_visu(year=2013, variety='Tremie13', sowing_date='10-29', nplants=15):
    (g, adel, weather, seq, rain_timing, 
     canopy_timing, septo_timing, recorder_timing, it_wheat, wheat_dir,
     wheat_is_loaded) = setup(sowing_date=str(year-1)+"-"+sowing_date+" 12:00:00", 
                              end_date=str(year)+"-08-01 00:00:00",
                              variety = variety, nplants = nplants)
    g = adel.setup_canopy(1400)
    adel.plot(g)
        
def temp_films():
    g, recorder = annual_loop_septo(year=2003, variety='Rht3', 
                                    sowing_date='10-15', nplants=30,
                                    degree_days_to_chlorosis=150., 
                                    save_images=True)
    g, recorder = annual_loop_septo(year=2003, variety='Mercia', 
                                    sowing_date='10-15', nplants=30,
                                    degree_days_to_chlorosis=150., 
                                    save_images=True)

def temp_plot_simu(df, multiply_sev = True, xaxis='degree_days',
                   leaves=None, only_severity=False, xlims = [800, 2200],
                   display_lesions=False, title=None, correct_fnl=False):
    from alinea.echap.disease.alep_septo_evaluation import data_reader, plot_confidence_and_boxplot, plot_one_sim
    from alinea.alep.alep_weather import plot_rain_and_temp
    import cPickle as pickle
    import matplotlib.pyplot as plt
#    from math import ceil
    df_sim = df.copy()
    if title is not None:
        ttle = False
    else:
        ttle = True
    if df_sim.iloc[-1]['date'].year == 2012:
        try:
            f = open('data_obs_2012.pckl')
            data_obs = pickle.load(f)
            f.close()
            f = open('weather_2012.pckl')
            weather = pickle.load(f)
            f.close()
        except:
            data_obs, weather = data_reader(year = 2012,
                                              variety = 'Tremie12',
                                              from_file = 'control')
            f = open('data_obs_2012.pckl', 'w')
            pickle.dump(data_obs, f)
            f.close()
            f = open('weather_2012.pckl', 'w')
            pickle.dump(weather, f)
            f.close()
        if leaves is None:
            leaves = range(1,7)
        (df_mean_obs, df_low, df_high, fig, 
         axs) = plot_confidence_and_boxplot(data_obs, weather, 
                                            leaves=leaves, variable='severity', 
                                            xaxis=xaxis, xlims=xlims,
                                            title=ttle, return_fig=True)
    elif df_sim.iloc[-1]['date'].year == 2013:
        try:
            f = open('data_obs_2013.pckl')
            data_obs = pickle.load(f)
            f.close()
            f = open('weather_2013.pckl')
            weather = pickle.load(f)
            f.close()
        except:
            data_obs, weather = data_reader(year = 2013,
                                              variety = 'Tremie13',
                                              from_file = 'control')
            f = open('data_obs_2013.pckl', 'w')
            pickle.dump(data_obs, f)
            f.close()
            f = open('weather_2013.pckl', 'w')
            pickle.dump(weather, f)
            f.close()
        if leaves is None:
            leaves = range(1,6)
        (df_mean_obs, df_low, df_high, fig, 
         axs) = plot_confidence_and_boxplot(data_obs, weather, 
                                            leaves=leaves, variable='severity', 
                                            xaxis=xaxis, xlims=xlims,
                                            title=ttle, return_fig=True)
    else:
        raise ValueError('Unavailable year')
        
    if multiply_sev:
        df_sim['severity']*=100
    if correct_fnl:
        df_sim = resample_fnl(data_obs, df_sim, weather, variable='severity')
    if only_severity:
        plot_one_sim(df_sim, 'severity', xaxis, axs, leaves, 'r')
    else:
        df_sim['sev_tot'] = df_sim['leaf_disease_area'] *100./ df_sim['leaf_area']
#        df_sim['nec_spo'] = (df_sim['ratio_nec']+df_sim['ratio_spo']+df_sim['ratio_empty']) *100
        plot_one_sim(df_sim, 'severity', xaxis, axs, leaves, 'm')
        plot_one_sim(df_sim, 'sev_tot', xaxis, axs, leaves, 'r')
        plot_one_sim(df_sim, 'leaf_green_area', xaxis, axs, leaves, 'g')
#        plot_one_sim(df_sim, 'nec_spo', xaxis, axs, leaves, 'c')
        
    if title is not None:
        plt.text(0.5, 0.98, title, fontsize=18,
                 transform=fig.transFigure, horizontalalignment='center')
                 
    if display_lesions:
        fig, axs = plt.subplots(int(ceil(len(leaves)/2.))+1, 2, figsize=(15,10))
        axs_les = axs[:-1]
        axs_weather = axs[-1]
        plot_one_sim(df_sim, 'nb_lesions', xaxis, axs_les, leaves, 'b', xlims)
        for ax in axs_weather.flat:
            plot_rain_and_temp(weather, xaxis=xaxis, ax=ax, xlims=xlims, title='')


def temp_plot_comparison(data_obs_2012, data_sim_2012, weather_2012, 
                         data_obs_2013, data_sim_2013, weather_2013,
                         multiply_sev = False, xaxis='age_leaf',
                         xlims = [0, 1200], title=None, correct_fnl=False):
    from alinea.echap.disease.alep_septo_evaluation import plot_confidence_and_boxplot, plot_one_sim
    import matplotlib.pyplot as plt
    leaves = [5,4,3,2,1]
    fig, axs = plt.subplots(1, len(leaves), figsize=(30., 5.))
    (df_mean_obs_2012, df_low_2012, df_high_2012, fig, 
     axs) = plot_confidence_and_boxplot(data_obs_2012, weather_2012, 
                                        leaves=leaves, variable='severity', 
                                        xaxis=xaxis, xlims=xlims,
                                        minimum_sample_size=15,
                                        linestyle='',
                                        title=False, return_fig=True, 
                                        fig=fig, axs=axs,
                                        display_legend=False,
                                        tight_layout=False,
                                        fixed_color='b')
    (df_mean_obs_2013, df_low_2013, df_high_2013, fig, 
     axs) = plot_confidence_and_boxplot(data_obs_2013, weather_2013, 
                                        leaves=leaves, variable='severity', 
                                        xaxis=xaxis, xlims=xlims,
                                        minimum_sample_size=15,
                                        linestyle='',
                                        title=False, return_fig=True, 
                                        fig=fig, axs=axs,
                                        display_legend=False,
                                        tight_layout=False,
                                        fixed_color='r')
    if multiply_sev:
        data_sim_2012['severity']*=100
        data_sim_2013['severity']*=100
    if correct_fnl:
        data_sim_2012 = resample_fnl(data_obs_2012, data_sim_2012, weather_2012, variable='severity')
        data_sim_2013 = resample_fnl(data_obs_2013, data_sim_2013, weather_2013, variable='severity')
    plot_one_sim(data_sim_2012, 'severity', xaxis, axs, leaves, 'b', linewidth=2)
    plot_one_sim(data_sim_2013, 'severity', xaxis, axs, leaves, 'r', linewidth=2)
    
    ticks = np.arange(xlims[0], xlims[1], 200)
    for ax in axs.flat:
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks)
        ax.set_ylabel('')
        ax.set_xlabel('Age of leaf ($^\circ$Cd)', fontsize=16)
    axs[0].set_ylabel('Severity (in %)', fontsize=18)

    proxy = [plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='-'),
             plt.Rectangle((0,0), 0,0, facecolor='k', alpha=0.3)]
    labels = ['Mean', 'Confidence\n interval']
    lgd = axs[-1].legend(proxy, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    locbox = (1.1, 0.7)
    locim = (1.12, 0.7)
    txtbox = TextArea('        Distribution', minimumdescent=False)
    ab1 = AnnotationBbox(txtbox, xy=locbox, xycoords="axes fraction", 
                         frameon=True, fontsize=32, box_alignment=(0., 0.5))
    axs[-1].add_artist(ab1)
    try:
        img = read_png('box.png')
        imagebox = OffsetImage(img, zoom=.25)
        ab2 = AnnotationBbox(imagebox, xy=locim, 
                             xycoords="axes fraction", frameon=False)         
        axs[-1].add_artist(ab2)
    except:
        pass

    proxy2 = [plt.Rectangle((0,0), 0,0, facecolor='b', alpha=0.3),
              plt.Rectangle((0,0), 0,0, facecolor='r', alpha=0.3)]
    labels2 = ['2012', '2013']
    lgd2 = axs[-2].legend(proxy2, labels2, bbox_to_anchor=(2.6, 0.4), loc='lower right', borderaxespad=0.)

    if title is not None:
        plt.text(0.5, 0.98, title, fontsize=18,
                 transform=fig.transFigure, horizontalalignment='center')
    fig.savefig('comparison', bbox_extra_artists=(lgd,ab1,ab2,lgd2), bbox_inches='tight')
#    plt.tight_layout()

def temp_plot_comparison_complete(data_obs_2012, data_sim_2012, weather_2012, 
                                 data_obs_2013, data_sim_2013, weather_2013,
                                 multiply_sev = True, xaxis='degree_days',
                                 xlims = [800, 2200], correct_fnl=False):
    from alinea.echap.disease.alep_septo_evaluation import plot_confidence_and_boxplot, plot_one_sim
    import matplotlib.pyplot as plt
    leaves = range(1,6)
    fig, axs = plt.subplots(len(leaves), 3, figsize=(16.5, 20))
    axs_left = np.array([[ax[0]] for ax in axs])
    axs_center = np.array([[ax[1]] for ax in axs])
    axs_right = np.array([[ax[2]] for ax in axs])
    
    # Left column 2012 vs 2013
    (df_mean_obs_2012, df_low_2012, df_high_2012, fig, 
     axs_left) = plot_confidence_and_boxplot(data_obs_2012, weather_2012, 
                                        leaves=leaves, variable='severity', 
                                        xaxis=xaxis, xlims=xlims,
                                        minimum_sample_size=15,
                                        title=False, return_fig=True, 
                                        fig=fig, axs=axs_left,
                                        linestyle='',
                                        fixed_color='b', display_legend=False,
                                        tight_layout=False)
    (df_mean_obs_2013, df_low_2013, df_high_2013, fig, 
     axs_left) = plot_confidence_and_boxplot(data_obs_2013, weather_2013, 
                                        leaves=leaves, variable='severity', 
                                        xaxis=xaxis, xlims=xlims,
                                        minimum_sample_size=15,
                                        title=False, return_fig=True, 
                                        fig=fig, axs=axs_left,
                                        linestyle='',
                                        fixed_color='r', display_legend=False,
                                        tight_layout=False)
    if multiply_sev:
        data_sim_2012['severity']*=100
        data_sim_2013['severity']*=100
    if correct_fnl:
        data_sim_2012_resampled = resample_fnl(data_obs_2012, data_sim_2012, weather_2012, variable='severity')
        data_sim_2013_resampled = resample_fnl(data_obs_2013, data_sim_2013, weather_2013, variable='severity')
        plot_one_sim(data_sim_2012_resampled, 'severity', xaxis, axs_left, leaves, 'b')
        plot_one_sim(data_sim_2013_resampled, 'severity', xaxis, axs_left, leaves, 'r')
    else:
        plot_one_sim(data_sim_2012, 'severity', xaxis, axs_left, leaves, 'b')
        plot_one_sim(data_sim_2013, 'severity', xaxis, axs_left, leaves, 'r')
    
    # Center column fnls for 2012
    colors = iter([(51/255.,204/255.,51/255.), (0,0,0.5)])
    axs_center = plot_confidence_and_boxplot_by_fnl(data_obs_2012, weather_2012, 
                                            leaves = leaves, variable = 'severity', 
                                            xaxis = xaxis, xlims=xlims,
                                            fig=fig, axs=axs_center,
                                            colors=colors, linestyles=iter(['','']),
                                            return_ax = True, 
                                            display_legend=False,
                                            tight_layout=False)
#    colors = iter([(0.5,0.5,1), (0,0,0.5)])
    colors = iter([(51/255.,204/255.,51/255.), (0,0,0.5)])
    for fnl in np.unique(data_sim_2012['fnl']):
        df = data_sim_2012[data_sim_2012['fnl']==fnl]
        plot_one_sim(df, 'severity', xaxis, axs_center, leaves, next(colors))
    
    # Right column fnls for 2013
    colors = iter([(1,102/255.,0.), (0.5,0,0)])
    axs_right = plot_confidence_and_boxplot_by_fnl(data_obs_2013, weather_2013, 
                                            leaves = leaves, variable = 'severity', 
                                            xaxis = xaxis, xlims=xlims,
                                            fig=fig, axs=axs_right,
                                            colors=colors,
                                            return_ax = True, 
                                            display_legend=False,
                                            tight_layout=False)
#    colors = iter([(1,0.5,0.5), (0.5,0,0)])
    colors =iter([(1,102/255.,0.), (0.5,0,0)])
    for fnl in np.unique(data_sim_2013['fnl']):
        df = data_sim_2013[data_sim_2013['fnl']==fnl]
        plot_one_sim(df, 'severity', xaxis, axs_right, leaves, next(colors))
        
    # Custom
    ticks = np.arange(xlims[0], xlims[1], 200)
    for ax in axs.flat:
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks)
        
    proxy = [plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='-'),
             plt.Rectangle((0,0), 0,0, facecolor='k', alpha=0.3)]
    labels = ['Mean', 'Confidence\n interval']
    lgd = axs[0][-1].legend(proxy, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    locbox = (1.1, 0.7)
    locim = (1.12, 0.7)
    txtbox = TextArea('        Distribution', minimumdescent=False)
    ab1 = AnnotationBbox(txtbox, xy=locbox, xycoords="axes fraction", 
                         frameon=True, fontsize=32, box_alignment=(0., 0.5))
    axs[0][-1].add_artist(ab1)
    try:
        img = read_png('box.png')
        imagebox = OffsetImage(img, zoom=.25)
        ab2 = AnnotationBbox(imagebox, xy=locim, 
                             xycoords="axes fraction", frameon=False)         
        axs[0][-1].add_artist(ab2)
    except:
        pass
    proxy2 = [plt.Rectangle((0,0), 0,0, facecolor='b', alpha=0.3),
              plt.Rectangle((0,0), 0,0, facecolor='r', alpha=0.3),
              plt.Rectangle((0,0), 0,0, facecolor=(51/255.,204/255.,51/255.), alpha=0.3),
              plt.Rectangle((0,0), 0,0, facecolor=(0,0,0.5), alpha=0.3),
              plt.Rectangle((0,0), 0,0, facecolor=(1,102/255.,0.), alpha=0.3),
              plt.Rectangle((0,0), 0,0, facecolor=(0.5,0,0), alpha=0.3)]
    labels2 = ['total 2012', 'total 2013', 'FNL12 2012', 'FNL13 2012', 'FNL11 2013', 'FNL12 2013']
    lgd2 = axs[1][-1].legend(proxy2, labels2, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.savefig('complete', bbox_extra_artists=(lgd,ab1,ab2,lgd2), bbox_inches='tight')
  
def temp_plot_by_leaf(df, multiply_sev = True, xaxis='degree_days',
                      leaves=None, xlims = [800, 2200], correct_fnl=False):
    from alinea.echap.disease.alep_septo_evaluation import data_reader
    from alinea.echap.disease.septo_data_treatment import plot_by_leaf as plt_lf
    import cPickle as pickle
    df_sim = df.copy()
    fig, ax = plt.subplots()
    if multiply_sev:
        df_sim['severity']*=100
        df_sim['severity_on_green']*=100
    if df_sim.iloc[-1]['date'].year == 2011:
        variety = df_sim.iloc[-1]['variety']
        try:
            f = open('data_obs_2011_'+variety.lower()+'.pckl')
            data_obs_2011 = pickle.load(f)
            f.close()
            f = open('weather_2011.pckl')
            weather_2011 = pickle.load(f)
            f.close()
        except:
            data_obs_2011, weather_2011 = data_reader(year = 2011,
                                                      variety = variety,
                                                      from_file = 'control')
            f = open('data_obs_2011_'+variety.lower()+'.pckl', 'w')
            pickle.dump(data_obs_2011, f)
            f.close()
            f = open('weather_2011.pckl', 'w')
            pickle.dump(weather_2011, f)
            f.close()
        if leaves is None:
            leaves = range(1,10)
        plt_lf(df_sim, weather_2011, variable='severity_on_green',
               xaxis=xaxis, leaves=leaves,
               xlims=xlims, linestyle='-', marker=None,
               ax=ax, minimum_sample_size=0)
        plt_lf(data_obs_2011, weather_2011, variable='septo_green', 
               xaxis=xaxis, leaves=leaves, error_bars=True,
               xlims=xlims,linestyle='', ax=ax, minimum_sample_size=15)
    elif df_sim.iloc[-1]['date'].year == 2012:
        try:
            f = open('data_obs_2012.pckl')
            data_obs_2012 = pickle.load(f)
            f.close()
            f = open('weather_2012.pckl')
            weather_2012 = pickle.load(f)
            f.close()
        except:
            data_obs_2012, weather_2012 = data_reader(year = 2012,
                                                      variety = 'Tremie12',
                                                      from_file = 'control')
            f = open('data_obs_2012.pckl', 'w')
            pickle.dump(data_obs_2012, f)
            f.close()
            f = open('weather_2012.pckl', 'w')
            pickle.dump(weather_2012, f)
            f.close()
        if leaves is None:
            leaves = range(1,7)
        if correct_fnl:
            df_sim = resample_fnl(data_obs_2012, df_sim, weather_2012, variable='severity')
        plt_lf(data_obs_2012, weather_2012, xaxis=xaxis, leaves=leaves,
               xlims=xlims,linestyle='--', ax=ax, minimum_sample_size=15)
#        leaves = range(1,10)
        plt_lf(df_sim, weather_2012, xaxis=xaxis, leaves=leaves,
               xlims=xlims, linestyle='-', marker=None,
               ax=ax, minimum_sample_size=0)
    elif df_sim.iloc[-1]['date'].year == 2013:
        try:
            f = open('data_obs_2013.pckl')
            data_obs_2013 = pickle.load(f)
            f.close()
            f = open('weather_2013.pckl')
            weather_2013 = pickle.load(f)
            f.close()
        except:
            data_obs_2013, weather_2013 = data_reader(year = 2013,
                                                      variety = 'Tremie13',
                                                      from_file = 'control')
            f = open('data_obs_2013.pckl', 'w')
            pickle.dump(data_obs_2013, f)
            f.close()
            f = open('weather_2013.pckl', 'w')
            pickle.dump(weather_2013, f)
            f.close()
        if leaves is None:
            leaves = range(1,6)
        if correct_fnl:
            df_sim = resample_fnl(data_obs_2013, df_sim, weather_2013, variable='severity')
        plt_lf(data_obs_2013, weather_2013, xaxis=xaxis, leaves=leaves,
               xlims=xlims,linestyle='--', ax=ax, minimum_sample_size=15)
        plt_lf(df_sim, weather_2013, xaxis=xaxis, leaves=leaves,
               xlims=xlims, linestyle='-', marker=None,
               ax=ax, minimum_sample_size=0)
    else:
        raise ValueError('Unavailable year')

def temp_plot_comparison_2011(data_obs_mercia, data_sim_mercia,
                             data_obs_rht3, data_sim_rht3, weather_2011,
                             multiply_sev = True, xaxis='age_leaf_vs_flag_emg',
                             xlims = [200, 1200]):
    from alinea.echap.disease.septo_data_treatment import plot_by_leaf as plt_lf
    fig, axs = plt.subplots(1,2, figsize = (10,5))
    leaves = [1,2,3]
    plt_lf(data_sim_mercia, weather_2011, variable='severity_on_green',
           xaxis=xaxis, leaves=leaves,
           xlims=xlims, linestyle='-', marker=None, title='', legend=False,
           ax=axs[0], minimum_sample_size=0, with_brewer=True)
    plt_lf(data_obs_mercia, weather_2011, variable='septo_green', 
           xaxis=xaxis, leaves=leaves, error_bars=True, title='Mercia',
           xlims=xlims,linestyle='', ax=axs[0], marker='d', legend=False, 
           minimum_sample_size=15, with_brewer=True)
    plt_lf(data_sim_rht3, weather_2011, variable='severity_on_green',
           xaxis=xaxis, leaves=leaves, title='',
           xlims=xlims, linestyle='-', marker=None,
           ax=axs[1], minimum_sample_size=0, with_brewer=True)
    plt_lf(data_obs_rht3, weather_2011, variable='septo_green', 
           xaxis=xaxis, leaves=leaves, error_bars=True, title='Rht3',
           xlims=xlims,linestyle='', ax=axs[1], marker='d', 
           minimum_sample_size=15, with_brewer=True)
    for ax in axs:
        ax.set_xlabel('Thermal time since\nflag leaf emergence ($^\circ$Cd)', fontsize=16)
    axs[0].set_ylabel('Severity (% of green leaf area)', fontsize=18)
    axs[1].set_ylabel('')
    fig.savefig('Mercia_Rht3', bbox_inches='tight')
    
def plot_states_leaf(df, leaf=2):
    from alinea.echap.disease.alep_septo_evaluation import plot_one_sim
    import matplotlib.pyplot as plt
    df_sim = df.copy()
    fig, ax = plt.subplots()
    ax = np.array([ax])
    inc = (0, 222/256., 0)
    nec = (185./256, 93./256, 37./256)
    xaxis = 'degree_days'
    leaves = [leaf]
    plot_one_sim(df_sim, 'leaf_green_area', xaxis, ax, leaves, 'g')
    plot_one_sim(df_sim, 'leaf_disease_area', xaxis, ax, leaves, 'r')
    plot_one_sim(df_sim, 'surface_inc', xaxis, ax, leaves, inc)
    plot_one_sim(df_sim, 'surface_chlo', xaxis, ax, leaves, 'y')
    plot_one_sim(df_sim, 'surface_nec', xaxis, ax, leaves, nec)
    df_sim['tot_spo'] = df_sim['surface_spo'] + df_sim['surface_empty']
    plot_one_sim(df_sim, 'tot_spo', xaxis, ax, leaves, 'k')
#    df_sim['non_prio'] = df_sim['surface_inc'] + df_sim['surface_chlo']
#    df_sim['prio'] = df_sim['surface_nec'] + df_sim['surface_spo'] + df_sim['surface_empty']
#    plot_one_sim(df_sim, 'non_prio', xaxis, ax, leaves, 'm')
#    plot_one_sim(df_sim, 'prio', xaxis, ax, leaves, 'c')
#    plot_one_sim(df_sim, 'surface_spo', xaxis, ax, leaves, 'm')
#    plot_one_sim(df_sim, 'surface_non_spo', xaxis, ax, leaves, 'k')
#    plot_one_sim(df_sim, 'surface_empty', xaxis, ax, leaves, 'c')
    
def temp_plot_simu_by_fnl(df, multiply_sev = True, xaxis='degree_days',
                   leaves=None, only_severity=False, xlims = [800, 2200],
                   display_lesions=False):
    from alinea.echap.disease.alep_septo_evaluation import (data_reader,
                                                            plot_confidence_and_boxplot_by_fnl,
                                                            plot_one_sim)
    import cPickle as pickle
    import matplotlib.pyplot as plt
#    from math import ceil
    df_sim = df.copy()
    if df_sim.iloc[-1]['date'].year == 2012:
        try:
            f = open('data_obs_2012.pckl')
            data_obs_2012 = pickle.load(f)
            f.close()
            f = open('weather_2012.pckl')
            weather = pickle.load(f)
            f.close()
        except:
            data_obs_2012, weather = data_reader(year = 2012,
                                                  variety = 'Tremie12',
                                                  from_file = 'control')
            f = open('data_obs_2012.pckl', 'w')
            pickle.dump(data_obs_2012, f)
            f.close()
            f = open('weather_2012.pckl', 'w')
            pickle.dump(weather, f)
            f.close()
        if leaves is None:
            leaves = range(1,7)
        axs = plot_confidence_and_boxplot_by_fnl(data_obs_2012, weather, 
                                            leaves = leaves, variable = 'severity', 
                                            xaxis = xaxis, xlims=xlims,
                                            return_ax = True)                        
    elif df_sim.iloc[-1]['date'].year == 2013:
        try:
            f = open('data_obs_2013.pckl')
            data_obs_2013 = pickle.load(f)
            f.close()
            f = open('weather_2013.pckl')
            weather = pickle.load(f)
            f.close()
        except:
            data_obs_2013, weather = data_reader(year = 2013,
                                                  variety = 'Tremie13',
                                                  from_file = 'control')
            f = open('data_obs_2013.pckl', 'w')
            pickle.dump(data_obs_2013, f)
            f.close()
            f = open('weather_2013.pckl', 'w')
            pickle.dump(weather, f)
            f.close()
        if leaves is None:
            leaves = range(1,6)
        axs = plot_confidence_and_boxplot_by_fnl(data_obs_2013, weather, 
                                            leaves = leaves, variable = 'severity', 
                                            xaxis = xaxis, xlims=xlims,
                                            return_ax = True)
    else:
        raise ValueError('Unavailable year')
        
    if multiply_sev:
        df_sim['severity']*=100
    colors = iter(['b', 'r'])
    for fnl in np.unique(df_sim['fnl']):
        df = df_sim[df_sim['fnl']==fnl]
        plot_one_sim(df, 'severity', xaxis, axs, leaves, next(colors))

def resample_fnl(df_obs, df_sim, weather, variable='severity'):
    from alinea.echap.disease.septo_data_reader import get_ratio_fnl
    
    def hack_convert(x):
        if np.isreal(x):
            return x
        else:
            return float(''.join(ch for ch in x if ch.isalnum()))
    
    def cross_product(date, num_lf, fnl):
        if fnl == main_fnl:
            if num_lf in df.columns:
                ratio_obs = 1-df.loc[np.datetime64(date, 'ns'), num_lf]
            else:
                ratio_obs = 1-ratio_obs_default
            return ratio_obs/(1-ratio_sim)
        else:
            if num_lf in df.columns:
                ratio_obs = df.loc[np.datetime64(date, 'ns'), num_lf]
            else:
                ratio_obs = ratio_obs_default
            return ratio_obs/ratio_sim
    
    df = get_ratio_fnl(df_obs, weather, 'severity')
    if df_sim.iloc[-1]['date'].year == 2012:
        ratio_sim = 0.33333333333333331
        ratio_obs_default = 0.21
        main_fnl = 13
        df.iloc[[-3,-2,-1], 4] = 0.
    elif df_sim.iloc[-1]['date'].year == 2013:
        ratio_sim = 0.28888888888888886
        ratio_obs_default = 20./43.
        main_fnl = 11

    df = df.reset_index()
    for date in set(df_sim['date']):
        row = df_sim[df_sim['date']==date].iloc[0]
        line = {'Date':row['date'], 'Degree days':row['degree_days']}
        for i in range(1,8):
            line[i] = np.nan
        df = df.append(pd.Series(line), ignore_index=True)
    df = df.sort('Date')
    for i in range(1, len(df.columns)-1):
        df.iloc[0,i+1] = df[i][~np.isnan(df[i])].values[0] if len(df[i][~np.isnan(df[i])])>0 else 0.
    for i in range(1, len(df.columns)-1):
        df.iloc[-1,i+1] = df[i][~np.isnan(df[i])].values[-1]
    df = df.interpolate()
    if df_sim.iloc[-1]['date'].year == 2012:
        df.loc[:,8][np.isnan(df[8])] = 0.
    df = df.set_index('Date')
    
    fun_resample = np.frompyfunc(cross_product, 3, 1)
    df_sim['cross_product'] = fun_resample(df_sim['date'], df_sim['num_leaf_top'], df_sim['fnl'])
    df_sim['cross_product'] = df_sim['cross_product'].astype(float)
    for col in ['leaf_area','leaf_green_area', 'leaf_length', 'leaf_senesced_length', 
               'nb_dispersal_units', 'nb_lesions', 'nb_lesions_on_green', 
               'surface_inc', 'surface_chlo', 'surface_nec', 'surface_nec_on_green', 
               'surface_spo', 'surface_spo_on_green', 'surface_empty', 
               'surface_empty_on_green', 'surface_dead', variable]:
        try:
            df_sim[col] = df_sim[col].astype(float)
        except:
            df_sim[col] = df_sim[col].map(hack_convert)
        df_sim[col] *= df_sim['cross_product']
    return df_sim
#for suffix in ['new_calib_age', 'new_calib_smin', 'new_calib_rate', 'new_calib_states', 'new_calib_states_2']:
#    data_sim_new_calib = get_aggregated_data_sim(variety = 'Tremie12', nplants = 15,
#                                                 sporulating_fraction=7e-2, suffix=suffix)
#    temp_plot_simu(data_sim_new_calib, multiply_sev = False)
#    plt.savefig('C:/Users/ggarin/Desktop/20150409_results/'+suffix[10:]+'.png')
#    plot_by_leaf(data_sim_new_calib, variable='severity')
#    plt.savefig('C:/Users/ggarin/Desktop/20150409_results/'+suffix[10:]+'_all.png')
#    plot_by_leaf(data_sim_new_calib, variable='nb_lesions')
#    plt.savefig('C:/Users/ggarin/Desktop/20150409_results/'+suffix[10:]+'_lesions.png')

def plot_variables_leaf(df, leaf=2, add_leaf_dates=False,
                        xaxis='age_leaf', xlims=[0,1400]):
    if add_leaf_dates==True:
        df = add_leaf_dates_to_data(df)
    fig, axs = plt.subplots(1,2, figsize=(15,6))
    leaves = [leaf]
    df = df[df['num_leaf_top']==leaf]
    plot_one_sim(df, 'leaf_area', xaxis, np.array(axs[0]), leaves, 'g', linestyle='-')
    plot_one_sim(df, 'leaf_green_area', xaxis, np.array(axs[0]), leaves, 'g', linestyle='--')
    plot_one_sim(df, 'leaf_necrotic_area', xaxis, np.array(axs[0]), leaves, 'r', linestyle='-')
    plot_one_sim(df, 'leaf_necrotic_area_on_green', xaxis, np.array(axs[0]), leaves, 'r', linestyle='--')

    df['ratio_green'] = 0.
    df['ratio_green'][df['leaf_area']>0] = df['leaf_green_area'][df['leaf_area']>0]*100./df['leaf_area'][df['leaf_area']>0]
    plot_one_sim(df, 'ratio_green', xaxis, np.array(axs[1]), leaves, 'g', linestyle='-')
    plot_one_sim(df, 'severity', xaxis, np.array(axs[1]), leaves, 'r', linestyle='-')
    df_fill1 = get_mean(df, column='severity', xaxis=xaxis)
    df_green = get_mean(df, column='leaf_green_area', xaxis=xaxis)
    date_end_sen = df_green[df_green[leaf]>0].index[-1]
    df_fill1[df_fill1.index>date_end_sen+210.] = 0. 
    axs[1].fill_between(df_fill1.index, df_fill1.values.flat, alpha=0.3,
                        color='None', edgecolor='r', hatch='//')
    
    axs[0].set_xlim(xlims)
    axs[1].set_xlim(xlims)
    axs[0].set_ylim([0,max(df_green[leaf])*1.05])
    axs[1].set_ylim([0,105])
    axs[0].set_ylabel('Leaf area (cm2)', fontsize=18)
    axs[1].set_ylabel('Severity (%)', fontsize=18)
    axs[0].set_xlabel('Age of leaf (degree days)', fontsize=18)
    axs[1].set_xlabel('Age of leaf (degree days)', fontsize=18)
    axs[0].legend(['Total leaf', 'Green leaf', 'Sporulating', 'Sporulating\nin green'], loc='best')
    axs[1].legend(['% Green leaf', '% Sporulating'], loc='best')
    plt.tight_layout()

def plot_3_weathers():
    import matplotlib.pyplot as plt
    f = open('weather_2011.pckl')
    weather_2011 = pickle.load(f)
    f.close()
    f = open('weather_2012.pckl')
    weather_2012 = pickle.load(f)
    f.close()
    f = open('weather_2013.pckl')
    weather_2013 = pickle.load(f)
    f.close()

    fig, axs = plt.subplots(2,3, figsize=(16,4))
    years = [2011,2012,2013]
    for i, w in enumerate([weather_2011, weather_2012, weather_2013]):
        yr = years[i]
        if i==0:
            ylabel = True
        else:
            ylabel = False
        plot_relative_humidity(w, ax = axs[0][i], title=yr, ylims=[0,100],
                               xlims=[0, 2500],  xlabel=False, ylabel=ylabel)
        if i==2:
            ylabel2 = True
        else:
            ylabel2 = False
        plot_rain_and_temp(w, ax = axs[1][i], title='', xlims=[0, 2500],
                           ylims_rain = [0,10], ylims_temp=[-10, 40],
                             xlabel=True, ylabel1=ylabel, ylabel2=ylabel2)
    fig.tight_layout()

def get_rmse(data_obs, data_sim):
    def get_mean_add_doy(df):
        df = get_mean(df, column='severity', xaxis='date')
        df['doy'] = [d.dayofyear for d in df.index]
        return df
    df_mean_obs, df_mean_sim = map(lambda df: get_mean_add_doy(df), 
                                   (data_obs, data_sim))
    df_mean_sim = df_mean_sim[df_mean_obs.columns][df_mean_sim['doy'].isin(df_mean_obs['doy'])]
    df_mean_sim.index = df_mean_obs.index
    df_mean_obs, df_mean_sim = map(lambda df: df.drop('doy',1), (df_mean_obs, df_mean_sim))
    return np.sqrt((df_mean_sim - df_mean_obs) ** 2).mean().mean()