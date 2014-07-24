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
from septo_decomposed import septo_disease
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
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, septo_timing)):
        canopy_iter, rain_iter, septo_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
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
    
        # Save outputs
        if septo_iter:
            date = septo_iter.value.index[0]
            flag_ids = get_flag_leaves(g)
            for plant, flag in flag_ids.iteritems():
                recorders[plant].update_vids(vids=flag.values()[0])
                recorders[plant].record(g, date, degree_days = septo_iter.value.degree_days[-1])
    
    for recorder in recorders.itervalues():
        recorder.get_complete_dataframe()
        recorder.get_audpc('necrosis_percentage')

    return np.mean([recorder.audpc for recorder in recorders.itervalues()])

def evaluate(values): 
    Y = np.empty([values.shape[0]])
    for i, X in enumerate(values):
        print i
        Y[i] = run_disease(sporulating_fraction=X[0],
                           degree_days_to_chlorosis=X[1],
                           Smin = X[2],
                           Smax = X[3],
                           growth_rate = X[4],
                           age_physio_switch_senescence=X[5],
                           density_dus_emitted_ref = X[6],
                           reduction_by_rain=X[7])
    return Y

# Set random seed (does not affect quasi-random Sobol sampling)
seed = 1
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

# Read the parameter range file and generate samples
param_file = 'poc_sensi.txt'
pf = read_param_file(param_file)

# Generate samples 
# param_values = morris_oat.sample(20, pf['num_vars'], num_levels = 10, grid_jump = 5)
param_values = morris_oat.sample(10, pf['num_vars'], num_levels = 10, grid_jump = 5)

# Samples are given in range [0, 1] by default. Rescale them to your parameter bounds. (If using normal distributions, use "scale_samples_normal" instead)
scale_samples(param_values, pf['bounds'])

# For Method of Morris, save the parameter values in a file (they are needed in the analysis)
np.savetxt('SGInput.txt', param_values, delimiter=' ')

# Run the model and save the output in a text file
# This will happen offline for external models
Y = evaluate(param_values)
np.savetxt("SGOutput.txt", Y, delimiter=' ')

# Perform the sensitivity analysis using the model output
Si = morris.analyze(param_file, 'SGInput.txt', 'SGOutput.txt', column = 0, conf_level = 0.95)
