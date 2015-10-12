"""
Main steps to run a simulation of a complex epidemics of septoria & brown rust
"""
# General imports
import pandas as pd
import numpy as np
import sys

# Imports for weather and scheduling of simulation
from alinea.alep.simulation_tools.simulation_tools import get_weather
from alinea.alep.alep_time_control import CustomIterWithDelays, septoria_filter_ddays
from alinea.astk.TimeControl import (time_filter, rain_filter, IterWithDelays,
                                     thermal_time_filter, DegreeDayModel,
                                     time_control)
                                     
# Imports for septoria
from alinea.alep.septoria import plugin_septoria
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum
from alinea.popdrops.alep_interface import (PopDropsSoilContamination, 
                                            PopDropsEmission,
                                            PopDropsTransport)

# Imports for rust
from alinea.alep.brown_rust import BrownRustFungus
from alinea.alep.inoculation import AirborneContamination
from alinea.alep.dispersal_transport import BrownRustDispersal

# Imports for both diseases
from alinea.alep.simulation_tools.simulation_tools import group_duplicates_in_cohort
from alinea.alep.protocol import infect, update, disperse, external_contamination
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.growth_control import SeptoRustCompetition
from alinea.alep.disease_outputs import SeptoRustRecorder

# Imports for wheat
from alinea.alep.simulation_tools.simulation_tools import (wheat_path, 
                                                           init_canopy, 
                                                           grow_canopy,
                                                           alep_echap_reconstructions,
                                                           alep_custom_reconstructions,
                                                           get_iter_rep_wheats,
                                                           get_filename,
                                                           get_data_sim)
from alinea.alep.architecture import set_properties
from alinea.alep.disease_outputs import plot_by_leaf

# Temp
from alinea.echap.disease.alep_septo_evaluation import *
from alinea.alep.disease_outputs import plot_by_leaf
from alinea.alep.simulation_tools.simulation_tools import add_leaf_dates_to_data
from alinea.adel.newmtg import adel_labels

def setup_simu(sowing_date="2000-10-15 12:00:00", start_date = None,
               end_date="2001-05-25 01:00:00", 
               variety = 'Mercia', nplants = 15, nsect = 7, 
               sporulating_fraction = 1e-4, Tmin = 0., Tmax = 25., WDmin = 10.,
               rain_min = 0.2, TT_delay = 20., septo_delay_dday = 10.,
               rust_dispersal_delay = 24, recording_delay = 24,
               record=True, layer_thickness_septo = 0.01,
               layer_thickness_rust = 1., rep_wheat = None, group_dus = True,
               leaf_duration = 2., **kwds):
    """ Setup the simulation 
    
    Note : kwds are used to modify disease parameters, if same name of parameter for septoria and 
    rust, need to specify in name with '_septoria' or '_rust', otherwise they are given the same 
    value.
    
    """
    def modif_params(fungus, suffix='_septoria'):
        pars = fungus.parameters()
        for par in pars.iterkeys():
            if par in kwds:
                pars[par] = kwds[par]
            elif par+'_septoria' in kwds:
                pars[par] = kwds[par+'_septoria']
        fungus.parameters(**pars)
        return fungus
    
    # Get weather
    weather = get_weather(start_date=sowing_date, end_date=end_date)
    
    # Set canopy
    it_wheat = 0
    if variety!='Custom':
        reconst = alep_echap_reconstructions(leaf_duration=leaf_duration)
        adel = reconst.get_reconstruction(name=variety, nplants=nplants, nsect=nsect)
    else:
        adel = alep_custom_reconstructions(variety='Tremie13', nplants=nplants, nsect=nsect, **kwds)
    year = int(end_date[:4])    
    wheat_dir = wheat_path(year, variety, nplants, nsect, rep_wheat)
    g, wheat_is_loaded = init_canopy(adel, wheat_dir, rain_and_light=True)
    domain = adel.domain
    domain_area = adel.domain_area
    convUnit = adel.convUnit
    
    # Manage temporal sequence  
    if start_date is None:
        start_date = sowing_date
    seq = pd.date_range(start=start_date, end=end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase=0.)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay=TT_delay)
    every_rain = rain_filter(seq, weather, rain_min = rain_min)
    every_rust_dispersal = time_filter(seq, delay=rust_dispersal_delay)
    every_recording = time_filter(seq, delay=recording_delay)
    septo_rust_filter = septoria_filter_ddays(seq, weather, delay = septo_delay_dday, rain_min = rain_min)
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')
    septo_dispersal_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    rust_dispersal_timing = IterWithDelays(*time_control(seq, every_rust_dispersal, weather.data))
    septo_rust_timing = CustomIterWithDelays(*time_control(seq, septo_rust_filter, weather.data), eval_time='end')
    recorder_timing = IterWithDelays(*time_control(seq, every_recording, weather.data))
    
    # Set up models
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del (sys.modules['alinea.alep.septoria_age_physio'])
    if 'alinea.alep.brown_rust' in sys.modules:
        del (sys.modules['alinea.alep.brown_rust'])
    kwds['nb_rings_by_state'] = 1.
    kwds['group_dus'] = group_dus
    septoria = plugin_septoria()
    septoria = modif_params(septoria, suffix='_septoria')
    brown_rust = BrownRustFungus()
    brown_rust = modif_params(brown_rust, suffix='_rust')
    if record==True:
        recorder = SeptoRustRecorder()
    else:
        recorder = None
    growth_controler = SeptoRustCompetition()
    infection_controler = BiotrophDUProbaModel()
    septo_inoculum = SoilInoculum(fungus=septoria, 
                                  sporulating_fraction=sporulating_fraction,
                                  domain_area=domain_area)
    septo_contaminator = PopDropsSoilContamination(fungus=septoria,
                                                   group_dus=group_dus,
                                                   domain=domain,
                                                   domain_area=domain_area, 
                                                   compute_star=False,
                                                   mutable = False)
    rust_contaminator = AirborneContamination(fungus=brown_rust,
                                              group_dus=group_dus,
                                              domain_area=domain_area,
                                              layer_thickness=layer_thickness_septo)
    septo_emitter = PopDropsEmission(domain=domain, compute_star = False)
    septo_transporter = PopDropsTransport(fungus=septoria, group_dus=group_dus,
                                          domain = domain, domain_area = domain_area,
                                          dh = layer_thickness_septo, convUnit = convUnit,
                                          compute_star = False)
    rust_dispersor = BrownRustDispersal(fungus = brown_rust,
                                        group_dus = group_dus,
                                        domain_area = adel.domain_area,
                                        layer_thickness=layer_thickness_rust)
    return (g, adel, brown_rust, septoria, canopy_timing, septo_dispersal_timing, 
            rust_dispersal_timing, septo_rust_timing, recorder_timing,
            recorder, growth_controler, infection_controler, septo_inoculum,
            septo_contaminator, rust_contaminator, septo_emitter, 
            septo_transporter, rust_dispersor, it_wheat,
            wheat_dir, wheat_is_loaded)
            
def annual_loop_septo_rust(year = 2013, variety = 'Tremie13', sowing_date = '10-29',
                           nplants = 15, nsect = 7, sporulating_fraction = 5e-3, 
                           density_dispersal_units = 150, 
                           layer_thickness_septo = 0.01, layer_thickness_rust = 1.,
                           record = True, output_file = None,
                           reset_reconst = True, rep_wheat = None, 
                           leaf_duration = 2., date_inoc_rust=1000., 
                           length_inoc_rust=300., **kwds):
    """ Simulate epidemics with canopy saved before simulation """
    if 'temp_min' in kwds:
        Tmin = kwds['temp_min']
    else:
        Tmin = 0.
    (g, adel, brown_rust, septoria, canopy_timing, septo_dispersal_timing, 
    rust_dispersal_timing, septo_rust_timing, recorder_timing,
    recorder, growth_controler, infection_controler, septo_inoculum,
    septo_contaminator, rust_contaminator, septo_emitter, 
    septo_transporter, rust_dispersor, it_wheat,
    wheat_dir, wheat_is_loaded) = setup_simu(sowing_date=str(year-1)+"-"+sowing_date+" 12:00:00", 
                                              end_date=str(year)+"-07-30 00:00:00",
                                              variety = variety, nplants = nplants, nsect=nsect,
                                              sporulating_fraction=sporulating_fraction,
                                              Tmin=Tmin, rep_wheat=rep_wheat, record=record,
                                              layer_thickness_septo=layer_thickness_septo,
                                              layer_thickness_rust=layer_thickness_rust,
                                              leaf_duration=leaf_duration, **kwds)

    for i, controls in enumerate(zip(canopy_timing, septo_dispersal_timing, rust_dispersal_timing,
                                     septo_rust_timing, recorder_timing)):
        (canopy_iter, septo_dispersal_iter, rust_dispersal_iter,
        septo_rust_iter, record_iter) = controls
        
        # Grow wheat canopy
        if canopy_iter:
            it_wheat += 1
            g = grow_canopy(g, adel, canopy_iter, it_wheat,
                        wheat_dir, wheat_is_loaded,rain_and_light=True)               
                
        # Get weather for date and add it as properties on leaves
        if septo_rust_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = septo_rust_iter.value.temperature_air.tolist(),
                           wetness_sequence = septo_rust_iter.value.wetness.tolist(),
                           relative_humidity_sequence = septo_rust_iter.value.relative_humidity.tolist(),
                           dd_sequence = septo_rust_iter.value.degree_days.tolist())
        if septo_dispersal_iter:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = septo_dispersal_iter.value.rain.mean(),
                           rain_duration = len(septo_dispersal_iter.value.rain) if septo_dispersal_iter.value.rain.sum() > 0 else 0.)
        # External contamination
        geom = g.property('geometry')
        if septo_dispersal_iter and len(geom)>0 and septo_dispersal_iter.value.rain.mean()>0.2:
            g = external_contamination(g, septo_inoculum, septo_contaminator, septo_dispersal_iter.value,
                                       domain=adel.domain, 
                                       domain_area=adel.domain_area)
#        if (rust_dispersal_iter and len(geom)>0 and
#            rust_dispersal_iter.value.index[0] > pd.to_datetime(str(year)+'-03-01') and
#            rust_dispersal_iter.value.index[-1] < pd.to_datetime(str(year)+'-03-15')):
#        if (rust_dispersal_iter and len(geom)>0 and
#            rust_dispersal_iter.value.degree_days.tolist()[-1] > date_inoc_rust and
#            rust_dispersal_iter.value.degree_days.tolist()[-1] < date_inoc_rust+length_inoc_rust):
#            print 'RUST INOC'
        if (rust_dispersal_iter and len(geom)>0):
            g = external_contamination(g, rust_contaminator, rust_contaminator, 
                                       density_dispersal_units=density_dispersal_units,
                                       domain_area=adel.domain_area)
        # Develop disease (infect for dispersal units and update for lesions)
        if septo_rust_iter:
            infect(g, septo_rust_iter.dt, infection_controler, label='LeafElement')
#            group_duplicates_in_cohort(g) # Additional optimisation (group identical cohorts)
            update(g, septo_rust_iter.dt, growth_controler, senescence_model=None, label='LeafElement')            
        # Disperse and wash
        if septo_dispersal_iter and len(geom)>0 and septo_dispersal_iter.value.rain.mean()>0.2:
            g = disperse(g, septo_emitter, septo_transporter, "septoria",
                         label='LeafElement', weather_data=septo_dispersal_iter.value,
                         domain=adel.domain, domain_area=adel.domain_area)
        if rust_dispersal_iter and len(geom)>0:
            g = disperse(g, rust_dispersor, rust_dispersor,
                         fungus_name = "brown_rust",
                         label='LeafElement', 
                         weather_data=rust_dispersal_iter.value,
                         domain_area=adel.domain_area)
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
        
def concat_inocs(inoc_septo=5e-3, inoc_rust=150.):
    inoc = str(inoc_septo)+'_'+str(inoc_rust)
    return inoc.replace('.', '_')
    
def run_reps_septo_rust(year = 2013, variety = 'Tremie13', 
                   nplants = 15, nsect = 7, sowing_date = '10-15',
                   sporulating_fraction = 5e-3, density_dispersal_units=150.,
                   layer_thickness_septo = 0.01, layer_thickness_rust = 1.,
                   nreps = 10, suffix = None, **kwds):
    df = pd.DataFrame()
    rep_wheats = get_iter_rep_wheats(year, variety, nplants, nsect, nreps)
    for rep in range(nreps):
        g, recorder = annual_loop_septo_rust(year=year, variety=variety,
                                            sowing_date=sowing_date,
                                            nplants=nplants, nsect=nsect,
                                            sporulating_fraction=sporulating_fraction,
                                            density_dispersal_units=density_dispersal_units,
                                            layer_thickness_septo=layer_thickness_septo, 
                                            layer_thickness_rust=layer_thickness_rust, 
                                            rep_wheat=next(rep_wheats), **kwds)
        df_ = recorder.data
        df_['rep'] = rep
        df = pd.concat([df, df_])
    inoc = concat_inocs(inoc_septo=sporulating_fraction, 
                        inoc_rust=density_dispersal_units)
    output_file = get_filename(fungus='septo_rust', year=year, variety=variety,
                               nplants=nplants, inoc=inoc, suffix = None)
    df.to_csv(output_file)
    
def plot_septo_rust(fungus = 'septo_rust', year = 2012,
                    variety = 'Tremie12', nplants = 15,
                    sporulating_fraction = 5e-4, density_dispersal_units=150., 
                    nreps = 10, variable = 'severity', xaxis = 'degree_days', 
                    leaves = range(1, 14), from_top = True,
                    plant_axis = ['MS'], error_bars = False, 
                    error_method = 'confidence_interval', marker = '', 
                    empty_marker = False, linestyle = '-', 
                    fixed_color = None, alpha = None, title = None, 
                    legend = True, xlabel = None, ylabel = None, 
                    xlims = None, ylims = None, ax = None,
                    return_ax = False, fig_size = (10,8)):
    inoc = concat_inocs(inoc_septo=sporulating_fraction, 
                        inoc_rust=density_dispersal_units)
    data_sim = get_data_sim(fungus=fungus, year=year, variety=variety,
                           nplants=nplants, inoc=inoc)
    if np.isnan(set(data_sim['variety'])[0]):
        data_sim['variety'] = variety
    plot_by_leaf(data_sim, variable=variable, xaxis=xaxis, leaves=leaves,
                 from_top=from_top, plant_axis=plant_axis, error_bars=error_bars,
                 error_method=error_method, marker=marker, 
                 empty_marker=empty_marker, linestyle=linestyle, 
                 fixed_color=fixed_color, alpha=alpha, title=title, 
                 legend=legend, xlabel=xlabel, ylabel=ylabel, xlims=xlims,
                 ylims=ylims, ax=ax, return_ax=return_ax, fig_size=fig_size)

def get_aggregated_data_sim(year = 2013, variety = 'Tremie13', nplants = 15, 
                            sporulating_fraction=5e-3, density_dispersal_units=150.,
                            num_leaf = 'num_leaf_top', suffix=None):
    from alinea.alep.simulation_tools.simulation_tools import get_data_sim
    inoc = concat_inocs(inoc_septo=sporulating_fraction, 
                        inoc_rust=density_dispersal_units)
    data_sim = get_data_sim(fungus='septo_rust', year=year,
                            variety=variety, nplants=nplants,
                            inoc=inoc, suffix=suffix)
    data_sim = get_data_without_death(data_sim, num_leaf=num_leaf)
    data_sim['severity'] *= 100
    data_sim['severity_on_green'] *= 100
    return data_sim

def example_climate(years = [2003,2012,2013], variety = 'Tremie13',
                    nplants = 15,  sowing_date = '10-29', 
                    inoc_septo = 5e-3, inoc_rust = 150.,
                    suffix = None, nreps=3, **kwds):
    scenarios_inoc = [(inoc_septo, 0), (0, inoc_rust), (inoc_septo, inoc_rust)]
    for yr in years:
        for inoc in scenarios_inoc:
            if yr!=2003 and not inoc in [(inoc_septo, 0), (0, inoc_rust)]:
                print 'pass'
                pass
            else:
                run_reps_septo_rust(year=yr, variety=variety, nplants=nplants,
                                    sowing_date=sowing_date,
                                    sporulating_fraction=inoc[0],
                                    density_dispersal_units=inoc[1],
                                    nreps=nreps, suffix=suffix, **kwds)
                                
def plot_example_climate(years = [2003,2012,2013], variety = 'Tremie13',
                        nplants = 15,  sowing_date = '10-15', 
                        inoc_rust = 150., inoc_septo = 5e-3, 
                        suffix = None, nreps=3):

    def plot_variable(df, variable='severity_septo', ax=None):
        plot_by_leaf(df, variable, xaxis = 'age_leaf_vs_flag_emg', 
                     ax=ax, ylims=[0, 1], xlims=[0, 1500], legend=False)

    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(4,len(years))
    scenarios_inoc = [(inoc_septo, 0), (0, inoc_rust), (inoc_septo, inoc_rust)]
    for i_yr, yr in enumerate(years):
        for inoc in scenarios_inoc:
            data_sim = get_aggregated_data_sim(year=yr, variety=variety, 
                                      nplants=nplants, 
                                      sporulating_fraction=inoc[0],
                                      density_dispersal_units=inoc[1],
                                      suffix=suffix)
            data_sim = add_leaf_dates_to_data(data_sim)
#            df_count = data_sim.groupby(['date', 'num_leaf_top']).count()
#            df_count = df_count.reset_index()
#            data_sim = data_sim[df_count['severity']==nplants*nreps]
            if inoc[0]==0:
                plot_variable(data_sim, variable='severity_rust', ax=axs[1][i_yr])
            elif inoc[1]==0:
                plot_variable(data_sim, variable='severity_septo_spo', ax=axs[0][i_yr])
            else:
                plot_variable(data_sim, variable='severity_septo_spo', ax=axs[2][i_yr])
                plot_variable(data_sim, variable='severity_rust', ax=axs[3][i_yr])
            