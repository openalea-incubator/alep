"""
Main steps to run a simulation of brown rust epidemics
"""
# General imports
import pandas as pd
import random as rd
import sys
from openalea.deploy.shared_data import shared_data

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
from alinea.alep.simulation_tools.simulation_tools import get_weather

# Imports for scheduling of simulation
from alinea.alep.alep_time_control import CustomIterWithDelays
from alinea.astk.TimeControl import (time_filter, IterWithDelays,
                                     thermal_time_filter, DegreeDayModel,
                                     time_control)
                                     
# Imports for disease
import alinea.alep
from alinea.alep.brown_rust import BrownRustFungus
from alinea.alep.simulation_tools.simulation_tools import group_duplicates_in_cohort
from alinea.alep.disease_outputs import BrownRustRecorder
from alinea.alep.growth_control import GeometricPoissonCompetition, SeptoRustCompetition
from alinea.alep.inoculation import AirborneContamination
from alinea.alep.protocol import infect, update, disperse, external_contamination
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.dispersal_transport import BrownRustDispersal

# Temp
from alinea.alep.disease_outputs import plot_by_leaf, get_data_without_death, get_synthetic_outputs_by_leaf

def setup_simu(sowing_date="2000-10-15 12:00:00", start_date = None,
               end_date="2001-05-25 01:00:00", 
               variety = 'Mercia', nplants = 30, nsect = 7,
               TT_delay = 20, dispersal_delay = 24,
               record=True, layer_thickness=1., rep_wheat = None, 
               leaf_duration=2., **kwds):
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
    
    # Manage temporal sequence  
    if start_date is None:
        start_date = sowing_date
    seq = pd.date_range(start=start_date, end=end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase=0.)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay=TT_delay)
    every_dispersal = time_filter(seq, delay=dispersal_delay)
    rust_filter = thermal_time_filter(seq, weather, TTmodel, delay=TT_delay)
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')
    dispersal_timing = IterWithDelays(*time_control(seq, every_dispersal, weather.data))
    rust_timing = CustomIterWithDelays(*time_control(seq, rust_filter, weather.data), eval_time='end')
    
    # Set up models
    if 'alinea.alep.brown_rust' in sys.modules:
        del (sys.modules['alinea.alep.brown_rust'])
    fungus = BrownRustFungus()
    fungus.parameters(**kwds)
    if record==True:
        recorder = BrownRustRecorder()
    else:
        recorder = None
#    growth_controler = GeometricPoissonCompetition()
    growth_controler = SeptoRustCompetition()
    infection_controler = BiotrophDUProbaModel(fungus=['brown_rust'])
    contaminator = AirborneContamination(fungus = fungus,
                                         group_dus = True,
                                         domain_area = adel.domain_area,
                                         layer_thickness=layer_thickness)
    dispersor = BrownRustDispersal(fungus = fungus,
                                   group_dus = True,
                                   domain_area = adel.domain_area,
                                   layer_thickness=layer_thickness)
    return (g, adel, fungus,  canopy_timing, dispersal_timing, rust_timing, 
            recorder, growth_controler, infection_controler, 
            contaminator, dispersor, it_wheat, wheat_dir, wheat_is_loaded)

def annual_loop_rust(year = 2013, variety = 'Tremie13', 
                     nplants = 15, nsect = 7, sowing_date = '10-29',
                     density_dispersal_units = 150, TT_delay=20,
                     record = True, output_file = None, layer_thickness=1.,
                     rep_wheat = True, leaf_duration=2., **kwds):
    """ Simulate an epidemics over the campaign. """
    # Setup simu
    (g, adel, fungus, canopy_timing, dispersal_timing, rust_timing, 
     recorder, growth_controler, infection_controler, 
     contaminator, dispersor, it_wheat, wheat_dir,
     wheat_is_loaded) = setup_simu(sowing_date=str(year-1)+"-"+sowing_date+" 12:00:00", 
                   end_date=str(year)+"-07-30 00:00:00",
                   variety = variety, nplants = nplants, nsect = nsect, 
                   TT_delay = TT_delay, dispersal_delay = 24, record=record, 
                   layer_thickness=layer_thickness, 
                   leaf_duration=leaf_duration, **kwds)
        
    # Simulation loop
    for i, controls in enumerate(zip(canopy_timing, 
                                     dispersal_timing, 
                                     rust_timing)):
        canopy_iter, dispersal_iter, rust_iter = controls
        # Grow wheat canopy
        if canopy_iter:
            it_wheat += 1
            g = grow_canopy(g, adel, canopy_iter, it_wheat,
                        wheat_dir, wheat_is_loaded)
        # Get weather for date and add it as properties on leaves
        if rust_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = rust_iter.value.temperature_air.tolist(),
                           wetness_sequence = rust_iter.value.wetness.tolist(),
                           dd_sequence = rust_iter.value.degree_days.tolist())
        # Simulate airborne contamination
        geom = g.property('geometry')
        if dispersal_iter and len(geom)>0:
            external_contamination(g, contaminator, contaminator, 
                                   density_dispersal_units=density_dispersal_units,
                                   domain_area=adel.domain_area)
        # Develop disease (infect for dispersal units and update for lesions)
        if rust_iter:
            infect(g, rust_iter.dt, infection_controler, label='LeafElement')
            group_duplicates_in_cohort(g) # Additional optimisation (group identical cohorts)
            update(g, rust_iter.dt, growth_controler, label='LeafElement')
        # Disperse disease
        if dispersal_iter and len(geom)>0:
            g = disperse(g, dispersor, dispersor,
                         fungus_name = "brown_rust",
                         label='LeafElement', 
                         weather_data=dispersal_iter.value,
                         domain_area=adel.domain_area)
        # Save outputs
        if rust_iter and record == True:
            date = rust_iter.value.index[-1]
            print date
            recorder.record(g, date, 
                            degree_days = rust_iter.value.degree_days[-1])
   
    if record == True:
        recorder.post_treatment(variety=variety)
        if output_file is not None:
            recorder.save(output_file)
        else:
            return g, recorder
    else:
        return g
        
def run_reps_rust(year = 2013, variety = 'Tremie13', 
                  nplants = 15, nsect = 7, sowing_date = '10-15',
                  density_dispersal_units = 150, 
                  layer_thickness = 1., nreps = 5, suffix = None, **kwds):
    df = pd.DataFrame()
    rep_wheats = get_iter_rep_wheats(year, variety, nplants, nsect, nreps)
    for rep in range(nreps):
        g, recorder = annual_loop_rust(year=year, variety=variety,
                                       sowing_date=sowing_date,
                                       nplants=nplants, nsect=nsect,
                                       density_dispersal_units=density_dispersal_units,
                                       layer_thickness=layer_thickness, 
                                       rep_wheat=next(rep_wheats), **kwds)
        df_ = recorder.data
        df_['rep'] = rep
        df = pd.concat([df, df_])
    output_file = get_filename(fungus='brown_rust', year=year, variety=variety,
                               nplants=nplants, inoc=density_dispersal_units,
                               suffix=suffix)
    df.to_csv(output_file)

def get_aggregated_data_sim(year = 2013, variety = 'Tremie13', nplants = 15, 
                            density_dispersal_units=150.,
                            num_leaf = 'num_leaf_top', suffix=None):
    from alinea.alep.simulation_tools.simulation_tools import get_data_sim
    data_sim = get_data_sim(fungus='brown_rust', year=year,
                            variety=variety, nplants=nplants,
                            inoc=density_dispersal_units,
                            suffix=suffix)
    data_sim = get_data_without_death(data_sim, num_leaf=num_leaf)
    data_sim['severity'] *= 100
    return data_sim    
    
def explore_scenarios(years = range(2000,2007), nplants=15, nreps=3,
                      parameters = {'scale_HS':0.9, 'scale_leafSenescence':0.9,
                                    'scale_stemDim':1.3, 'scale_stemRate':1.1,
                                    'tiller_probability':0.8, 'scale_leafDim_length':1.2,
                                    'scale_leafDim_width':1.2, 'scale_leafRate':1.1,
                                    'scale_fallingRate':0.8}):
    parameters['reference']=1.
    for param in parameters:
        kwds = {k:1. if k!=param else v for k,v in parameters.iteritems()}
        scale_leafRate = 1.5*kwds.pop('scale_leafRate')
        for yr in years:
            run_reps_rust(year=yr, variety='Custom', sowing_date='10-29',
                   nplants=nplants, density_dispersal_units = 50.,
                   scale_leafRate = scale_leafRate,
                   sowing_density=0.85*250., plant_density=0.85*250.,
                   suffix='scenario_20151611_'+param+'_'+str(yr), nreps=nreps, **kwds)

def plot_explore_scenarios(years = range(1999,2007), nplants=15, variable='max_severity', 
                      parameters = {'scale_HS':0.9, 'scale_leafSenescence':0.9,
                                    'scale_stemDim':1.3, 'scale_stemRate':1.1,
                                    'tiller_probability':0.8, 'scale_leafDim_length':1.2,
                                    'scale_leafDim_width':1.2, 'scale_leafRate':1.1,
                                    'scale_fallingRate':0.8},
                                    force_rename={}, title='quanti_rust'):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1,2, figsize=(16,8))
    axSev = axs[0]
    axRef = axs[1]
    parameters['reference']=1.
    markers = {'reference':'o', 
                  'scale_HS':'x',
                  'scale_leafSenescence':'^',
                  'scale_stemDim':'s',
                  'scale_stemRate':'p',
                  'tiller_probability':'*',
                  'scale_leafDim_length':'d',
                  'scale_leafDim_width':'D', 
                  'scale_leafRate':'h',
                  'scale_fallingRate':'H'}
    colors = {'reference':'b', 
              'scale_HS':'r',
              'scale_leafSenescence':(182/255., 91./255, 22/255.),
              'scale_stemDim':'m',
              'scale_stemRate':'m',
              'tiller_probability':'k',
              'scale_leafDim_length':'g',
              'scale_leafDim_width':'g', 
              'scale_leafRate':'c',
              'scale_fallingRate':'y'}
    labels = []
    proxys = []
    refs = {}
    for yr in years:
        suffix='scenario_sd_reference_'+str(yr)
        df_ref = get_aggregated_data_sim(variety='Custom', nplants=nplants,
                                         density_dispersal_units=150,
                                         suffix=suffix, year=yr)
        df = get_synthetic_outputs_by_leaf(df_ref)
        refs[yr] = df[df['num_leaf_top']==1][variable].values[0]
        del df_ref
        del df
    
    sort_refs = sorted([(value,key) for (key,value) in refs.items()], reverse=True)
    rank_years = {ref[1]:x for x, ref in enumerate(sort_refs)}
    for use_ref in [False, True]:
        for param in sorted(parameters.keys()):
            color = colors[param]
            marker = markers[param]
            for yr in years:
                x = rank_years[yr]
                suffix='scenario_sd_'+param+'_'+str(yr)
                df_sim = get_aggregated_data_sim(variety='Custom', nplants=nplants,
                                                 density_dispersal_units=150,
                                                 suffix=suffix, year=yr)
                df = get_synthetic_outputs_by_leaf(df_sim)
                y = df[df['num_leaf_top']==1][variable].values[0]
                if use_ref==False:
                    axSev.plot([x], [y], color=color, marker=marker, markersize=8)
                else:
                    y = y/refs[yr] if refs[yr]>0 else 1.
                    axRef.plot([x], [y], color=color, marker=marker, markersize=8)
                del df_sim
                del df
            if param in force_rename:
                labels += [force_rename[param]]
            elif param=='reference':
                labels += [r"$\mathit{Reference}$"]
            else:
                labels += [param]
            if use_ref==False:
                proxys += [plt.Line2D((0,1),(0,0), color=color, 
                                      marker=marker, linestyle='None')]

    str_ranked_years = [str(t[1]) for t in sorted({v:k for k,v in rank_years.iteritems()}.items())]    
    axRef.set_xlim([-1, len(years)])
    axRef.set_xticks([-1] + range(len(years))+ [len(years)+1])
    axRef.set_xticklabels(['']+str_ranked_years+[''])
    axSev.set_xlim([-1, len(years)])
    axSev.set_xticks([-1] + range(len(years))+ [len(years)+1])
    axSev.set_xticklabels(['']+str_ranked_years+[''])
    axRef.grid(alpha=0.5)
    axSev.grid(alpha=0.5)
    axSev.set_ylabel('Maximum severity (%)', fontsize=16)
    axRef.set_ylabel('Variation of maximum severity', fontsize=16)
    lgd = axRef.legend(proxys, labels, numpoints=1, 
                        bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
#    plt.subplots_adjust(hspace=0.05)
    fig.savefig(title, bbox_extra_artists=(lgd,), bbox_inches='tight')
    return fig, axs, lgd