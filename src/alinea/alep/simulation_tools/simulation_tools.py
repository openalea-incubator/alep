# -*- coding: iso-8859-15  -*- 
"""
Main tools to run a simulation of disease epidemics on wheat
"""
import os
import pandas as pd
import numpy as np
import random as rd

# Imports for weather
from alinea.alep.alep_weather import (wetness_rapilly,
                                      linear_degree_days,
                                      plot_rain_and_temp,
                                      plot_wetness_and_temp)
from alinea.alep.alep_time_control import *
from alinea.astk.TimeControl import *
from alinea.echap.weather_data import *
import alinea.septo3d
import alinea.alep
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather

# Imports for disease
from alinea.alep.disease_outputs import plot_by_leaf

# Imports for wheat
from alinea.adel.astk_interface import AdelWheat
from alinea.adel.Stand import AgronomicStand
from alinea.adel.AdelR import devCsv
import alinea.adel.plantgen_extensions as pgen_ext
from alinea.echap.architectural_reconstructions import (EchapReconstructions, pdict,
                                                        reconstruction_parameters,
                                                        HS_fit, GL_fits,
                                                        dimension_fits,
                                                        leafshape_fits)
from alinea.adel.newmtg import move_properties
from alinea.caribu.caribu_star import rain_and_light_star

# Tools for weather ###########################################################
def get_weather(start_date="2010-10-15 12:00:00", end_date="2011-08-01 01:00:00"):
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
    weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':linear_degree_days}, start_date=start_date, base_temp=0., max_temp=30.)
    return weather
    
def add_leaf_dates_to_weather(weather, variety='Tremie12'):
    weather.data['variety'] = variety
    weather.data['fnl'] = 12
    weather.data['num_leaf_bottom'] = 12
    weather.data['num_leaf_top'] = 1
    weather.data = add_leaf_dates_to_data(weather.data)

def plot_weather_annual_loop(year = 2012, variety='Tremie12',
                             sowing_date = '10-15', 
                             xlims = [0, 2500], title = None,
                             xaxis = 'degree_days', yaxis='rain_and_temp',
                             ylims = None):
    sowing_date = str(year-1)+"-"+sowing_date+" 12:00:00"
    end_date = str(year)+"-07-01 00:00:00"
    weather = get_weather(start_date=sowing_date,
                          end_date=end_date)
    if xaxis != 'degree_days':
        add_leaf_dates_to_weather(weather, variety=variety)
    if yaxis=='rain_and_temp':
        plot_rain_and_temp(weather, xaxis=xaxis, xlims=xlims, ylims_rain=ylims, title=title)
    elif yaxis=='wetness_and_temp':
        plot_wetness_and_temp(weather, xaxis=xaxis, xlims=xlims, title=title)
    elif yaxis=='relative_humidity':
        plot_relative_humidity(weather, xaxis=xaxis, xlims=xlims, title=title)

# Tools for wheat #############################################################
def count_available_canopies(year, variety, nplants, nsect):
    filepath = str(shared_data(alinea.alep)/'wheat_reconstructions')
    list_files = os.listdir(filepath)
    start = variety.lower()+'_'+str(int(year))+'_'+\
            str(nplants)+'pl_'+str(nsect)+'sect'
    return len([f for f in list_files if f.startswith(start)])

def wheat_path(year, variety, nplants, nsect, rep):
    if rep is None:
        rep = ''
    filepath = str(shared_data(alinea.alep)/'wheat_reconstructions')
    return filepath+'/'+variety.lower()+'_'+str(int(year))+'_'+\
            str(nplants)+'pl_'+str(nsect)+'sect_rep'+str(rep)

def init_canopy(adel, wheat_dir, rain_and_light=True):
    if os.path.exists(wheat_dir):
        wheat_is_loaded = True
        it_wheat = 0
        g, TT = adel.load(it_wheat, dir=wheat_dir)
    else:
        wheat_is_loaded = False
        g = adel.setup_canopy(age=0.)
        if rain_and_light==True:
            rain_and_light_star(g, light_sectors = '1', 
                                domain=adel.domain, convUnit=adel.convUnit)
    return g, wheat_is_loaded

def grow_canopy(g, adel, canopy_iter, it_wheat,
                wheat_dir, wheat_is_loaded=True, rain_and_light=True):
    if wheat_is_loaded:
        newg, TT = adel.load(it_wheat, dir=wheat_dir)
        move_properties(g, newg)
        return newg
    else:
        g = adel.grow(g, canopy_iter.value)
        if rain_and_light==True:
            rain_and_light_star(g, light_sectors = '1', 
                                domain=adel.domain, convUnit=adel.convUnit)
        return g
    
def alep_echap_reconstructions(keep_leaves=False, leaf_duration=2.,
                               single_nff=False, variability=True):
    pars = reconstruction_parameters()
    pars['density_tuning'] = pdict(None)
    pars['density_tuning']['Tremie12'] = 0.85
    pars['density_tuning']['Tremie13'] = 0.85
    if keep_leaves==True:
        pars['pgen_base'] = {'TT_hs_break':None, 'inner_params':{'DELAIS_PHYLL_SEN_DISP':20}}
    #CF : ici si tu veux ajuster la duree d'allongement des feuilles, en particulier pour avoir des emergences plus precoces:
    #    Une feuille emerge (0.8 * leafDuration * phyllochron) dd  avant HS
    pars['adel_pars']['leafDuration'] = leaf_duration #2 = valeur normale ; 3 = valeur calibree pour septo
    # Temp
    pars['adel_pars']['dHS_en'] -= 1.
#    pars['adel_pars']['stemDuration'] += 1.
    #
    
    # Manage NFF
    if single_nff==False:
        nff_in_disease_data = {'Mercia':None, 'Rht3':None,
                               'Tremie12': {'12':0.21, '13':0.79},
                               'Tremie13': {'11':23./43, '12':20./43.}}
    else:       
        nff_in_disease_data = {'Mercia':None, 'Rht3':None,
                               'Tremie12': {'12':0., '13':1.},
                               'Tremie13': {'11':0., '12':1.}}
    pars['nff_probabilities'] = nff_in_disease_data

    # Manage variability
    if variability == True:
        my_variability = {'Mercia': None,
                          'Rht3': None,
                          'Tremie12': 0,
                          'Tremie13': 10}
    else:
        my_variability = {'Mercia': 0,
                          'Rht3': 0,
                          'Tremie12': 0,
                          'Tremie13': 0} 
    pars['std_emergence'] = my_variability

    
    reconst = EchapReconstructions(reset_data=True, pars=pars)
    
    # G.Garin : temporary commented : new way to set nff
#    reconst.axepop_fits['Tremie12'].MS_probabilities = {12:0.21, 13:0.79}
#    reconst.axepop_fits['Tremie13'].MS_probabilities = {11:23./43, 12:20./43.}

    # Temp    
#    reconst.HS_fit['Tremie12'].TT_hs_0 -= 35
#    reconst.HS_fit['Tremie13'].TT_hs_0 -= 35

    # Temp
    reconst.GL_fits['Tremie12'].GL_bolting = reconst.GL_fits['Tremie13'].GL_bolting   
    return reconst
    
def alep_custom_reconstructions(variety='Tremie13', nplants=30, 
                                sowing_density=250.,
                                plant_density=250., inter_row=0.15,
                                nsect=7, seed=1, 
                                scale_HS = 1,
                                scale_leafDim_length = 1,
                                scale_leafDim_width = 1,
                                scale_leafRate=1,
                                scale_stemDim = 1,
                                scale_stemRate=1, 
                                scale_fallingRate=1,
                                scale_leafSenescence=1,
                                scale_tillering=1,
                                **kwds):
    parameters = reconstruction_parameters()
    stand = AgronomicStand(sowing_density=sowing_density,
                            plant_density=plant_density, 
                            inter_row=inter_row, noise=0.04, 
                            density_curve_data = None)       
    n_emerged = nplants
    
    MS_leaves_number_probabilities={'11': 0.3, '12': 0.7}
    tiller_proba = {'T1':1., 'T2':0.5, 'T3':0.5, 'T4':0.3}
    ears_per_plant = 2.5
    
    if 'tiller_probability' in kwds:
        tiller_proba = {k:min(1, p*kwds['tiller_probability']) for k,p in tiller_proba.items()}

    if 'proba_main_nff' in kwds:
        MS_leaves_number_probabilities = {'11':1-kwds['proba_main_nff'], '12':kwds['proba_main_nff']}
    
    # scale _tillering
    tiller_proba = {k:min(1, p * scale_tillering) for k,p in tiller_proba.items()}
    em = pgen_ext.TillerEmission(primary_tiller_probabilities=tiller_proba)
    reg = pgen_ext.TillerRegression(ears_per_plant = ears_per_plant * scale_tillering)
    
    #generate plants
    HSfit = HS_fit()[variety]
    axp = pgen_ext.AxePop(MS_leaves_number_probabilities = MS_leaves_number_probabilities,
                          Emission=em, Regression=reg, std_em=HSfit.std_TT_hs_0)
    plants = axp.plant_list(n_emerged)
    
    # Measure and save Green leaf durations on the reference
    adel_pars = parameters['adel_pars']

    HSfit.mean_nff = axp.mean_nff()
    Dimfit = dimension_fits(HS_fit(), **parameters)[variety]
    GLfit = GL_fits(HS_fit(), **parameters)[variety]
    pgen = pgen_ext.PlantGen(HSfit=HSfit, GLfit=GLfit, Dimfit=Dimfit, adel_pars = adel_pars)
    axeT, dimT, phenT = pgen.adelT(plants)
    axeT = axeT.sort(['id_plt', 'id_cohort', 'N_phytomer'])
    GreenDuration1_ref = axeT['TT_sen_phytomer1'] - axeT['TT_em_phytomer1']
    phenT = phenT.sort(['id_phen', 'index_phytomer'])
    dTTsen_ref = phenT['dTT_sen_phytomer']
    dTTem_ref = phenT['dTT_em_phytomer']

    
    # Modify fits
    
    # HS_fit
    
    # petit commentaire : scale_HS determine ligulation rate + emergence rate
    # la duree emergence-ligulation est impactee du fait de ce changement (inevitable si on impose leaf_rate)
    # mais aussi en raison de scale_leafRate(non souhaite => compensation ci dessous)
    HS_ref = HSfit.a_cohort
    leafDuration_ref = adel_pars['leafDuration']
    HSfit.a_cohort /= scale_HS
    # on compense en decalant l'emergence des plante
    HSfit.TT_hs_0 += (1 - adel_pars['fracLeaf']) / HS_ref * leafDuration_ref * (float(scale_HS) / scale_leafRate - 1) 
    
    
    # Modify leaf dimension and growth rate  
    
    adel_pars['leafDuration'] = ((scale_leafDim_length*scale_leafDim_width) / scale_HS)* \
                                (leafDuration_ref / scale_leafRate)
    Dimfit.scale['L_blade'] *= scale_leafDim_length
    Dimfit.scale['W_blade'] *= scale_leafDim_width
                                    
    # Modify stem dimension and growth rate
    stemDuration_ref = adel_pars['stemDuration']
    adel_pars['stemDuration'] = (scale_stemDim / scale_HS)* \
                                (stemDuration_ref / scale_stemRate)
    Dimfit.scale['L_internode'] *= scale_stemDim

    
    
    # Scale GL : GL flag, bolting and start_senescence : apply factor
    # TEST EXTREME WITH LEAF RATE
    #GLfit.GL_bolting *= scale_leafSenescence
    #GLfit.GL_flag *= scale_leafSenescence
    #GLfit.n0 *= scale_leafSenescence

    # Modify  green leaf durations and falling rate
    pgen = pgen_ext.PlantGen(HSfit=HSfit, GLfit=GLfit, Dimfit=Dimfit, adel_pars = adel_pars)
    axeT, dimT, phenT = pgen.adelT(plants)
    axeT = axeT.sort(['id_plt', 'id_cohort', 'N_phytomer'])    
    phenT = phenT.sort(['id_phen', 'index_phytomer'])
    
    # Modify senescence
    axeT['TT_sen_phytomer1'] = axeT['TT_em_phytomer1'] + scale_leafSenescence * GreenDuration1_ref
    phenT['dTT_sen_phytomer'] = dTTsen_ref  + (phenT['dTT_em_phytomer'] - dTTem_ref)
    devT = devCsv(axeT, dimT, phenT)
    
    leaves = leafshape_fits(**parameters)[variety]
    bins_ref = leaves.bins
    bins_ref[-1] = 21.
    bins = [x * scale_HS / scale_fallingRate if i_x!=0 else x*scale_HS for i_x,x in enumerate(bins_ref)]
    leaves.bins = bins
    
    return AdelWheat(nplants = nplants, nsect=nsect, devT=devT, stand = stand , 
                    seed=seed, sample='sequence', leaves = leaves, run_adel_pars = adel_pars)
    
def make_canopy(year = 2013, variety = 'Tremie13', sowing_date = '10-29',
                nplants = 15, nsect = 7, nreps=10, fixed_rep=None, delay = 20.):
    """ Simulate and save canopy (prior to simulation). """    
    def make_canopy_one_rep(rep=0):
        reconst = alep_echap_reconstructions()
        adel = reconst.get_reconstruction(name=variety, nplants=nplants, nsect=nsect)    
        domain = adel.domain
        convUnit = adel.convUnit
        g = adel.setup_canopy(age=0.)
        rain_and_light_star(g, light_sectors = '1', domain = domain, convUnit = convUnit)
        it_wheat = 0
        wheat_dir = wheat_path(year, variety, nplants, nsect, rep)
        adel.save(g, it_wheat, dir=wheat_dir)
        for i, canopy_iter in enumerate(canopy_timing):
            if canopy_iter:
                it_wheat += 1
                g = adel.grow(g, canopy_iter.value)
                rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
                adel.save(g, it_wheat, dir=wheat_dir)
                
    # Manage weather and scheduling
    start_date=str(year-1)+"-"+sowing_date+" 12:00:00"
    end_date=str(year)+"-07-30 00:00:00"
    weather = get_weather(start_date=start_date, end_date=end_date)
    seq = pandas.date_range(start=start_date, end=end_date, freq='H')
    TTmodel = DegreeDayModel(Tbase = 0)
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay=delay)
    canopy_timing = CustomIterWithDelays(*time_control(seq, every_dd, weather.data), eval_time='end')

    # Simulate and save wheat development
    if fixed_rep is None:
        for rep in range(nreps):
            make_canopy_one_rep(rep=rep)
    else:
        make_canopy_one_rep(rep=fixed_rep)
        
def get_iter_rep_wheats(year = 2013, variety = 'Tremie13',
                        nplants = 15, nsect = 7, nreps=5):
    nb_can = count_available_canopies(year, variety, nplants, nsect)
    if nb_can >= nreps:
        return iter(rd.sample(list(range(nb_can)), nreps))
    else:
        return iter([None for rep in range(nreps)])
        
def get_missing(numbers, start, end):
    """ Find missing numbers in list compared to range(start,end) """
    all_numbers = set(range(start, end))
    missing = all_numbers - set(numbers)
    return np.array(list(missing))
        
def get_fnl_by_plant(data):
    """ Group plants from data by final leaf number (fnl) 
    
        :Returns:
        - dict([(fnl, [plant_ids])])
    """
    # Get data for last date and top leaves
    df = data.copy()
    df = df.xs(max(set(df.index)))
    df = df[df['num_leaf_top']==1]
    # Group plant ids by final number of leaves (FNL)
    grps = {fnl:[plant for plant in df[df['num_leaf_bottom']==fnl]['num_plant']] 
            for fnl in set(df['num_leaf_bottom'])}
    # Find fnl if missing value for a leaf on last notation
    if len(sum(iter(grps.values()),[]))<max(data['num_plant']):
        miss = get_missing(np.sort(sum(iter(grps.values()),[])), start=1,
                           end=max(data['num_plant']))
        for m in miss:
            if len(data[data['num_plant']==m])>0:
                fnl = max(data['num_leaf_bottom'][data['num_plant']==m]) + \
                      min(data['num_leaf_top'][data['num_plant']==m]) - 1
                grps[fnl].append(m)
    return grps
    
def add_leaf_dates_to_data(df, correct_leaf_number=True, force_mean_fnl=False):
    def get_date_em(num_leaf_bottom=1, fnl=None):
        return hs_fit.TTemleaf(num_leaf_bottom, nff=fnl)
    def get_date_lig(num_leaf_bottom=1, fnl=None):
        return hs_fit.TTligleaf(num_leaf_bottom, nff=fnl)
    def get_date_lig_flag(fnl=None):
        return hs_fit.TTligleaf(fnl, nff=fnl)
    def add_dates(df):
        fun_em = np.frompyfunc(get_date_em, 2, 1)
        fun_lig = np.frompyfunc(get_date_lig, 2, 1)
        df['date_emergence_leaf'] = fun_em(df['num_leaf_bottom'], df['fnl'])
        df['date_ligulation_leaf'] = fun_lig(df['num_leaf_bottom'], df['fnl'])
        if force_mean_fnl==False:
            df['date_emergence_flag_leaf'] = [hs_fit.TTemleaf(fnl, nff=fnl) for fnl in df['fnl']]
            df['date_ligulation_flag_leaf'] = [hs_fit.TTligleaf(fnl, nff=fnl) for fnl in df['fnl']]
        else:
            df['date_emergence_flag_leaf'] = hs_fit.TTemleaf(hs_fit.mean_nff, nff=None)[0]
            df['date_ligulation_flag_leaf'] = hs_fit.TTligleaf(hs_fit.mean_nff, nff=None)[0]
        return df
    
    if not 'fnl' in df.columns:
        fnls = get_fnl_by_plant(df)
        df['fnl'] = [{vv:k for k,v in fnls.items() for vv in v}[x] for x in df['num_plant']]
    # Correction only for specific data of disease measurements (Grignon Tremie 2012)
    if correct_leaf_number == True:
        df['num_leaf_bottom'][df['fnl']==14] -= 1
        df['fnl'][df['fnl']==14] -= 1
    variety = np.unique(df['variety'])[0].title()
    if variety=='Custom':
        variety='Tremie13'
    hs_fit = HS_fit()[variety]
    df = add_dates(df)
    df['age_leaf'] = df['degree_days'] - df['date_emergence_leaf']
    df['age_leaf_lig'] = df['degree_days'] - df['date_ligulation_leaf']
    df['age_leaf_vs_flag_lig'] = df['degree_days'] - df['date_ligulation_flag_leaf']
    df['age_leaf_vs_flag_emg'] = df['degree_days'] - df['date_emergence_flag_leaf']
    return df

# Tools for disease ###########################################################
from collections import defaultdict
def group_duplicates_in_cohort(g):
    def _get_index_duplicates(seq):
        dd = defaultdict(list)
        for i,item in enumerate(seq):
            dd[item].append(i)
        dups = [idx for key,idx in dd.items() if len(idx)>1 and key==0.]
        if len(dups)>0:
            return dups[-1]
        else:
            return []
    
    def group_lesions(les):
        fungi = set([l.fungus.name for l in les])
        for f in fungi:
            les_f = [l for l in les if l.fungus.name.startswith(f)]
            new_l = les_f[0].fungus.lesion()
            new_l.position = sum([l.position for l in les_f], [])
        return new_l
    
    lesions = g.property('lesions')
    for vid, les in lesions.items():
        ages = [l.age_tt for l in les]
        if len(les)!=len(set(ages)):
            idxs = _get_index_duplicates(ages)
            if len(idxs)>0:
                new_les = group_lesions([les[i] for i in idxs])
                les = les[:idxs[0]]
                les.append(new_les)
                lesions[vid] = les
        
# Tools for plotting results ##################################################
def get_filename(fungus = 'brown_rust', year=2012, variety = 'Tremie12', 
                 nplants = 15, inoc=300, suffix=None):
    inoc = str(inoc)
    inoc = inoc.replace('.', '_')
    if suffix is None:
        suffix=''
    else:
        suffix='_'+suffix
    filename= variety+'_'+str(year)+'_'+str(nplants)+'pl_inoc'+inoc+suffix+'.csv'
    fungus_path = fungus + '_simulations'
    return str(shared_data(alinea.alep)/fungus_path/filename)   
        
def get_data_sim(fungus = 'brown_rust', year = 2012, variety = 'Tremie12',
                 nplants = 15, inoc=150, suffix=None):
    filename = get_filename(fungus=fungus, year=year, variety=variety, 
                            nplants=nplants, inoc=inoc, suffix=suffix)
    df = pd.read_csv(filename, parse_dates={'datetime':[1]}, index_col=[0])
    df.index.names = ['date']
    df = df.drop('Unnamed: 0', axis=1)
    if 'age_leaf' not in df.columns:
        df = add_leaf_dates_to_data(df, correct_leaf_number=False)
    return df
    
def plot_simu_results(fungus = 'brown_rust', year = 2012,
                      variety = 'Tremie12', nplants = 15, inoc=300, 
                      nreps = 5, variable = 'severity', xaxis = 'degree_days', 
                      leaves = list(range(1, 14)), from_top = True,
                      plant_axis = ['MS'], error_bars = False, 
                      error_method = 'confidence_interval', marker = '', 
                      empty_marker = False, linestyle = '-', 
                      fixed_color = None, alpha = None, title = None, 
                      legend = True, xlabel = None, ylabel = None, 
                      xlims = None, ylims = None, ax = None,
                      return_ax = False, fig_size = (10,8), suffix=None):
   data_sim = get_data_sim(fungus=fungus, year=year, variety=variety,
                           nplants=nplants, inoc=inoc, suffix=suffix)          
   plot_by_leaf(data_sim, variable=variable, xaxis=xaxis, leaves=leaves,
                from_top=from_top, plant_axis=plant_axis, error_bars=error_bars,
                error_method=error_method, marker=marker, 
                empty_marker=empty_marker, linestyle=linestyle, 
                fixed_color=fixed_color, alpha=alpha, title=title, 
                legend=legend, xlabel=xlabel, ylabel=ylabel, xlims=xlims,
                ylims=ylims, ax=ax, return_ax=return_ax, fig_size=fig_size)