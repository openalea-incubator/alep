# params = {'sporulating_fraction':0.7777800000005, 'degree_days_to_chlorosis': 164.44444444444446, 'threshold_spores': 1000,
 # 'production_rate': 100000, 'DEAD': 5, 'reduction_by_rain': 0.77777777777777779, 
 # 'Smin': 0.0055444444444444456, 'temp_min': 10.0, 'rain_events_to_empty': 3, 'EMPTY':4, 'temp_max': 30.0,
 # 'wd_min': 10.0, 'degree_days_to_sporulation': 20.0, 'INCUBATING': 0, 'NECROTIC': 2, 'growth_rate': 0.0011499999999999995,
 # 'group_dus': True, 'density_dus_emitted_ref': 2555.5555555555557, 'SPORULATING': 3, 'epsilon': 0.001,
 # 'loss_rate': 0.008333333333333333, 'threshold_spo': 0.0001, 'CHLOROTIC': 1, 'basis_for_dday': -2.0,
 # 'rh_min': 85.0, 'degree_days_to_necrosis': 110.0, 'Smax': 1.0, 'loss_delay': 120.0, 'delta_age_ring': 20.0,
 # 'nb_rings_by_state': 1, 'age_physio_switch_senescence': 0.56000000000000005}
 
 # General imports
import numpy as np
import random as rd
import pandas
import pickle
import sys

# Imports for sensitivity analysis
from SALib.sample import morris_oat
from SALib.analyze import morris
from SALib.util import scale_samples, read_param_file

# Imports for wheat
from alinea.adel.newmtg import move_properties
from alinea.echap.architectural_reconstructions import reconst_db
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves

# Imports for weather
from alinea.alep.alep_weather import wetness_rapilly, basic_degree_days
from alinea.alep.alep_time_control import *
from alinea.astk.TimeControl import *
from alinea.echap.weather_data import *

# Imports for alep septoria
from alinea.alep.protocol import *
from septo_decomposed import septo_disease, load_ids
from alinea.alep.disease_outputs import SeptoRecorder

def get_flag_leaves(g):
    labels = g.property('label')
    stems = [id for id,lb in labels.iteritems() if lb.startswith('MS')]
    blades = [id for id,lb in labels.iteritems() if lb.startswith('blade')]
    flag_leaves_ids = {}
    ind_plant = 0
    for st in stems:
        ind_plant += 1
        flag_leaves_ids['P%d' % ind_plant] = {}
        nff = int(g.node(st).properties()['nff'])
        lf = [bl for bl in blades if bl>st][nff-1]
        stem_elt = g.node(lf).components()[0].index()
        flag_leaves_ids['P%d' % ind_plant]['F1'] = range(stem_elt+1, stem_elt+nsect+1)
    return flag_leaves_ids
        
def run_disease(sporulating_fraction=1e-2,
                degree_days_to_chlorosis=220.,
                Smin = 0.03,
                Smax = 0.3,
                growth_rate = 0.0006,
                age_physio_switch_senescence=1,
                density_dus_emitted_ref = 1.79e3,
                reduction_by_rain=0.5):
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
        
    inoculum, contaminator, infection_controler, growth_controler, emitter, transporter = septo_disease(domain, 
    domain_area, sporulating_fraction=sporulating_fraction,degree_days_to_chlorosis=degree_days_to_chlorosis,
    Smin = Smin, Smax = Smax, growth_rate = growth_rate, age_physio_switch_senescence=age_physio_switch_senescence,
    density_dus_emitted_ref = density_dus_emitted_ref, reduction_by_rain=reduction_by_rain)
    
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    septo_timing = CustomIterWithDelays(*time_control(seq, septo_filter, weather.data), eval_time='end')
    it = 0
    g,TT = adel.load(it, dir = './adel/adel_saved2')
    
    recorders = {}
    flag_ids = get_flag_leaves(g)
    for plant, flag in flag_ids.iteritems():
        recorders[plant] = SeptoRecorder(vids=flag.values()[0], group_dus=True)
    
    # recorders = {}
    # it_septo = 0
    # leaf_sectors = load_ids(it_septo, dir = './adel/adel_saved2')
    # print 'leaf_sectors: '+ str(leaf_sectors['P1']['F1'])
    # for plant in leaf_sectors:
        # recorders[plant] = {}
        # for leaf, lf_sectors in leaf_sectors[plant].iteritems():
            # recorders[plant][leaf] = SeptoRecorder(vids=lf_sectors, group_dus=True)
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, septo_timing)):
        canopy_iter, rain_iter, septo_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            print canopy_iter.value.index[0]
            it += 1
            newg,TT = adel.load(it, dir = './adel/adel_saved2')
            move_properties(g, newg)
            g = newg
        
        # Get weather for date and add it as properties on leaves
        if septo_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = septo_iter.value.temperature_air,
                           wetness = septo_iter.value.wetness.mean(),
                           relative_humidity = septo_iter.value.relative_humidity.mean())
        if rain_iter:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = rain_iter.value.rain.mean(),
                           rain_duration = len(rain_iter.value.rain) if rain_iter.value.rain.sum() > 0 else 0.)
        
        # External contamination
        geom = g.property('geometry')
        if rain_iter and len(geom)>0:
            if rain_iter.value.rain.mean()>0. and rain_iter.value.degree_days[-1]<1000:
                g = external_contamination(g, inoculum, contaminator, rain_iter.value)

        # Develop disease (infect for dispersal units and update for lesions)
        if septo_iter:
            infect(g, septo_iter.dt, infection_controler, label='LeafElement')
            update(g, septo_iter.dt, growth_controler, senescence_model=None, label='LeafElement')
            
        # Disperse and wash
        if rain_iter and len(geom)>0:
            if rain_iter.value.rain.mean()>0.:
                g = disperse(g, emitter, transporter, "septoria", label='LeafElement', weather_data=rain_iter.value)
    
        # if septo_iter:
            # it_septo += 1
            # date = septo_iter.value.index[0]
            # print date
            # leaf_sectors = load_ids(it_septo, dir = './adel/adel_saved2')
            # print 'leaf_sectors: ' + str(leaf_sectors['P1']['F1'])
            # for plant in recorders:
                # for lf, recorder in recorders[plant].iteritems():
                    # recorder.update_vids(vids=leaf_sectors[plant][lf])
                    # recorder.record(g, date, degree_days = septo_iter.value.degree_days[-1])
                    
    # for plant in recorders:
        # for recorder in recorders[plant].itervalues():
            # recorder.get_complete_dataframe()
            # recorder.get_audpc(variable = 'necrosis_percentage')

        # Save outputs
        if septo_iter:
            date = septo_iter.value.index[0]
            flag_ids = get_flag_leaves(g)
            for plant, flag in flag_ids.iteritems():
                recorders[plant].update_vids(vids=flag.values()[0])
                recorders[plant].record(g, date, degree_days = septo_iter.value.degree_days[-1])
    
    for recorder in recorders.itervalues():
        recorder.get_complete_dataframe()
        recorder.get_audpc()

    return np.mean([recorder.audpc for recorder in recorders.itervalues()])
    # return recorders
 
# Set random seed (does not affect quasi-random Sobol sampling)
# seed = 1
seed = 0
np.random.seed(seed)
rd.seed(seed)

# Initialize wheat plant
Mercia = reconst_db['Mercia']
nsect = 5
pgen, adel, domain, domain_area, convUnit, nplants = Mercia(nplants = 3, nsect=nsect)
    
# Manage weather
weather = Boigneville_2010_2011()
start_date = "2010-10-15 12:00:00"
end_date = "2011-08-31 01:00:00"
weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
weather.check(varnames=['degree_days'], models={'degree_days':basic_degree_days}, start_date=start_date)
weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':basic_degree_days}, start_date=start_date, base_temp=-2.)

# Define the schedule of calls for each model
seq = pandas.date_range(start = start_date, end = end_date, freq='H')
TTmodel = DegreeDayModel(Tbase = 0)
every_dd = thermal_time_filter(seq, weather, TTmodel, delay = 10)
every_rain = rain_filter(seq, weather)
every_dd_or_rain = filter_or([every_dd, every_rain])
canopy_timing = IterWithDelays(*time_control(seq, every_dd_or_rain, weather.data))
septo_filter = septo_infection_filter(seq, weather, every_rain)
 
# run_disease(sporulating_fraction=0.7777800000005,
                # degree_days_to_chlorosis=164.44444444444446,
                # Smin = 0.0055444444444444456,
                # Smax = 1.0,
                # growth_rate = 0.0011499999999999995,
                # age_physio_switch_senescence=0.56000000000000005,
                # density_dus_emitted_ref = 2555.5555555555557,
                # reduction_by_rain=0.77777777777777779)

recorders = run_disease(sporulating_fraction=7.780000000000000457e-03,
                        degree_days_to_chlorosis=2.311111111111111143e+02,
                        Smin = 1.000000000000000048e-04,
                        Smax = 3.399999999999999689e-01,
                        growth_rate = 4.449999999999999963e-03,
                        age_physio_switch_senescence = 7.800000000000000266e-01,
                        density_dus_emitted_ref = 1.222222222222222172e+03,
                        reduction_by_rain = 1.000000000000000000e+00)                
                        
# recorders = run_disease()