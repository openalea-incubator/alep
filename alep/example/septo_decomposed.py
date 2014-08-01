"""

Run annual loop for the model of septoria in alep on the basis of 'annual_loop_decomposed' in echap

"""
import pandas
try:
    import cPickle as pickle
except:
    import pickle
import sys

# Imports for wheat and rain/light interception
from alinea.adel.newmtg import move_properties, adel_ids
from alinea.echap.architectural_reconstructions import reconst_db
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves
from alinea.caribu.caribu_star import rain_and_light_star

# Imports for weather
from alinea.alep.alep_weather import wetness_rapilly, basic_degree_days
from alinea.alep.alep_time_control import *
from alinea.astk.TimeControl import *
from alinea.echap.weather_data import *
import alinea.septo3d
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather

# Imports for alep septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum, Septo3DSoilContamination
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport
from alinea.alep.growth_control import PriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.disease_outputs import initiate_all_adel_septo_recorders, plot_severity_by_leaf

def get_weather(start_date="2010-10-15 12:00:00", end_date="2011-06-20 01:00:00"):
    start_yr = start_date[2:4]
    end_yr = end_date[2:4]
    weather_file = 'meteo'+ start_yr + '-' + end_yr + '.txt'
    meteo_path = shared_data(alinea.septo3d, weather_file)
    return Weather(data_file=meteo_path)

def setup(start_date="2010-10-15 12:00:00", end_date="2011-06-20 01:00:00", nplants = 3, nsect = 5, disc_level = 20):
    # Initialize wheat plant
    Mercia = reconst_db['Mercia']
    pgen, adel, domain, domain_area, convUnit, nplants = Mercia(nplants = nplants, nsect = nsect, disc_level = disc_level)

    # Manage weather
    if start_date[:4]=='2010':
        weather = Boigneville_2010_2011()
    else:
        weather = get_weather(start_date, end_date)
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
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    septo_timing = CustomIterWithDelays(*time_control(seq, septo_filter, weather.data), eval_time='end')
    return adel, domain, domain_area, convUnit, weather, seq, rain_timing, canopy_timing, septo_timing

def septo_disease(domain, domain_area, sporulating_fraction, **kwds):
    fungus = plugin_septoria()
    fungus.parameters(group_dus=True, nb_rings_by_state=1, **kwds)
    inoculum = SoilInoculum(fungus, sporulating_fraction=sporulating_fraction, domain_area=domain_area)
    contaminator = Septo3DSoilContamination(domain=domain, domain_area=domain_area)
    growth_controler = PriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    emitter = PopDropsEmission(domain=domain)
    transporter = PopDropsTransport(domain=domain, domain_area=domain_area)
    return inoculum, contaminator, infection_controler, growth_controler, emitter, transporter

def save_ids(ids, it, dir = './adel/adel_saved'):
    stored_ids = dir+'/ids_'+str(it)+'.pckl'
    f_ids = open(stored_ids, 'w')
    pickle.dump(ids, f_ids)
    f_ids.close()

def load_ids(it, dir = './adel/adel_saved'):
    stored_ids = dir+'/ids_'+str(it)+'.pckl'
    f_ids = open(stored_ids)
    ids = pickle.load(f_ids)
    f_ids.close()
    return ids
    
def get_leaf_ids(g, nsect = 5):
    labels = g.property('label')
    stems = [id for id,lb in labels.iteritems() if lb.startswith('MS')]
    blades = [id for id,lb in labels.iteritems() if lb.startswith('blade')]
    leaf_sectors = {}
    ind_plant = 0
    for st in stems:
        ind_plant += 1
        leaf_sectors['P%d' % ind_plant] = {}
        nff = int(g.node(st).properties()['nff'])
        leaves = [bl for bl in blades if bl>st][:nff]
        ind_lf = nff+1
        for lf in leaves:
            ind_lf -= 1
            stem_elt = g.node(lf).components()[0].index()
            leaf_sectors['P%d' % ind_plant]['F%d' % ind_lf] = range(stem_elt+1, stem_elt+nsect+1)
    return leaf_sectors
    
def make_canopy(start_date = "2010-10-15 12:00:00", end_date = "2011-06-20 01:00:00",
                nplants = 3, nsect = 5, disc_level = 20, dir = './adel/adel_3'):
                
    adel, domain, domain_area, convUnit, weather, seq, rain_timing, canopy_timing, septo_timing = setup(
            start_date = start_date, end_date = end_date, nplants = nplants, nsect = nsect, disc_level = disc_level)
            
    g = adel.setup_canopy(age=0.)
    rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
    it_wheat = 0
    # it_septo = 0
    adel.save(g, it_wheat, dir=dir)
    # ids = get_leaf_ids(g, nsect=nsect)
    # save_ids(ids, it_septo, dir=dir)
    for i, controls in enumerate(zip(canopy_timing, septo_timing)):
        canopy_iter, septo_iter = controls
        if canopy_iter:
            it_wheat += 1
            g = adel.grow(g, canopy_iter.value)
            rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            adel.save(g, it_wheat, dir=dir)
        
        # if septo_iter:
            # it_septo += 1
            # ids = get_leaf_ids(g, nsect=nsect)
            # save_ids(ids, it_septo, dir=dir)          
            
######### TEMP #####################

# def test_pickle_dump_canopy(adel, g):
    # for ind in range(10):
        # filename = './test_pickle'
        # adel.save(g, ind, dir=filename)

# @profile
# def test_pickle_load_canopy(adel):
    # for ind in range(10):
        # filename = './test_pickle'
        # adel.load(ind, dir=filename)
        
# adel, domain, domain_area, convUnit, weather, seq, rain_timing, canopy_timing, septo_timing = setup(
    # start_date = "2010-10-15 12:00:00", end_date = "2011-06-20 01:00:00", nplants = 3, nsect = 5, disc_level = 30)
# g = adel.setup_canopy(age=1500.)
# test_pickle_dump_canopy(adel, g)
# test_pickle_load_canopy(adel)
    
def get_leaf_ids(g):
    nsect = 5
    labels = g.property('label')
    stems = [id for id,lb in labels.iteritems() if lb.startswith('MS')]
    blades = [id for id,lb in labels.iteritems() if lb.startswith('blade')]
    leaf_sectors = {}
    ind_plant = 0
    for st in stems:
        ind_plant += 1
        leaf_sectors['P%d' % ind_plant] = {}
        nff = int(g.node(st).properties()['nff'])
        leaves = [bl for bl in blades if bl>st][:nff]
        ind_lf = nff+1
        for lf in leaves:
            ind_lf -= 1
            leaf_sectors['P%d' % ind_plant]['F%d' % ind_lf] = range(lf+2, lf+nsect+2)
    return leaf_sectors

# @profile    
def run_canopy(start_date = "2011-03-01 12:00:00", end_date = "2011-04-01 01:00:00",
                nplants = 3, nsect = 5, disc_level = 20, get_ids = 0):
    Mercia = reconst_db['Mercia']
    pgen, adel, domain, domain_area, convUnit, nplants = Mercia(nplants = nplants, nsect = nsect, disc_level = disc_level)

    # Manage weather
    weather = Boigneville_2010_2011()
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
    
    g = adel.setup_canopy(age=1000.)
    it_wheat = 0
    for i, canopy_iter in enumerate(canopy_timing):
        if canopy_iter:
            it_wheat += 1
            print canopy_iter.value.index[-1]
            g = adel.grow(g, canopy_iter.value)
            vids = adel_ids(g)
            vids2 = get_leaf_ids(g)

# run_canopy()
            
######### TEMP #####################
            
def save_leaf_ids(start_date = "2010-10-15 12:00:00", end_date = "2011-06-20 01:00:00",
                nplants = 3, nsect = 5, dir = './adel/adel_3'):
       
    adel, domain, domain_area, convUnit, weather, seq, rain_timing, canopy_timing, septo_timing = setup(
            start_date = start_date, end_date = end_date, nplants = nplants, nsect = nsect)
    g = adel.setup_canopy(age=0.)
    it_septo = 0
    ids = get_leaf_ids(g, nsect=nsect)
    save_ids(ids, it_septo, dir=dir)
    for i, controls in enumerate(zip(canopy_timing, septo_timing)):
        canopy_iter, septo_iter = controls
        if canopy_iter:
            g = adel.grow(g, canopy_iter.value)
            
        if septo_iter:
            it_septo += 1
            print it_septo
            ids = get_leaf_ids(g, nsect=nsect)
            save_ids(ids, it_septo, dir=dir)

def run_disease(start_date = "2010-10-15 12:00:00", end_date = "2011-06-20 01:00:00", nplants = 3, nsect = 5,
                disc_level = 20, dir = './adel/adel_saved', sporulating_fraction = 1e-4, adel = None, 
                domain = None, domain_area = None, convUnit = None, weather = None, seq = None, 
                rain_timing = None, canopy_timing = None, septo_timing = None, **kwds):

    if any(x==None for x in [adel, domain, domain_area, convUnit, weather, seq, rain_timing, canopy_timing, septo_timing]):
        adel, domain, domain_area, convUnit, weather, seq, rain_timing, canopy_timing, septo_timing = setup(
                start_date = start_date, end_date = end_date, nplants = nplants, nsect = nsect, disc_level = disc_level)
            
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
    inoculum, contaminator, infection_controler, growth_controler, emitter, transporter = septo_disease(domain, domain_area, sporulating_fraction, **kwds)
    it_wheat = 0
    it_septo = 0
    g,TT = adel.load(it_wheat, dir=dir)
    
    # Prepare saving of outputs
    recorders = initiate_all_adel_septo_recorders(g, nsect)
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, septo_timing)):
        canopy_iter, rain_iter, septo_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            it_wheat += 1
            print it_wheat
            newg,TT = adel.load(it_wheat, dir=dir)
            move_properties(g, newg)
            g = newg
            leaf_ids = adel_ids(g)
        
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
        
        # if canopy_iter:
            # scene = plot_severity_by_leaf(g)

        # Save outputs
        if septo_iter:
            # it_septo += 1
            date = septo_iter.value.index[-1]
            # leaf_sectors = load_ids(it_septo, dir=dir)
            for plant in recorders:
                for lf, recorder in recorders[plant].iteritems():
                    recorder.update_vids_with_labels(adel_ids = leaf_ids)
                    recorder.record(g, date, degree_days = septo_iter.value.degree_days[-1])
                    
    for plant in recorders:
        for recorder in recorders[plant].itervalues():
            recorder.get_complete_dataframe()
            recorder.get_audpc()
            
    return g, recorders

######### TEMP #####################
# run_disease()    
######### TEMP #####################

def run_and_save():
    frac = ['1e-4', '1e-3', '1e-2', '1e-1']
    for fr in frac:
        for i_sim in range(15):
            print fr, i_sim
            g, recorder = run_disease(sporulating_fraction = float(fr))
            stored_rec = '.\mercia\\recorder_fspo_'+fr+'_'+str(i_sim)+'.pckl'
            f_rec = open(stored_rec, 'w')
            pickle.dump(recorder, f_rec)
            f_rec.close()
            del recorder
            
def suite_run_and_save():
    frac = ['1e-5', '1e-4', '1e-3', '1e-2', '1e-1']
    for fr in frac:
        for i_sim in range(15):
            if not (fr == '1e-3' and i_sim < 13):
                print fr, i_sim
                g, recorder = run_disease(sporulating_fraction = float(fr))
                stored_rec = '.\mercia\\recorder_fspo_'+fr+'_'+str(i_sim)+'.pckl'
                f_rec = open(stored_rec, 'w')
                pickle.dump(recorder, f_rec)
                f_rec.close()
                del recorder
                
    incub = ['120', '220', '320']
    for inc in incub:
        for i_sim in range(10):
            print inc, i_sim
            g, recorder = run_disease(sporulating_fraction = 1e-2, degree_days_to_chlorosis = float(inc))
            stored_rec = '.\mercia\\recorder_incub_'+inc+'_'+str(i_sim)+'.pckl'
            f_rec = open(stored_rec, 'w')
            pickle.dump(recorder, f_rec)
            f_rec.close()
            del recorder
            
    growth_rates = ['6e-10', '220', '320']
    for inc in incub:
        for i_sim in range(10):
            print inc, i_sim
            g, recorder = run_disease(sporulating_fraction = 1e-2, degree_days_to_chlorosis = float(inc))
            stored_rec = '.\mercia\\recorder_incub_'+inc+'_'+str(i_sim)+'.pckl'
            f_rec = open(stored_rec, 'w')
            pickle.dump(recorder, f_rec)
            f_rec.close()
            del recorder

def run_and_save_years():
    years = [1998, 2003, 2004]
    for yr in years:
        start_date = str(yr)+"-10-15 12:00:00"
        end_date = str(yr+1)+"-06-20 01:00:00"
        path = './adel/adel_'+str(yr)
        for i_sim in range(15):
            print yr, i_sim
            g, recorder = run_disease(start_date = start_date,
                                      end_date = end_date,
                                      dir = path,
                                      sporulating_fraction = 1e-3)
            stored_rec = '.\mercia\\recorder_year_'+yr+'_'+str(i_sim)+'.pckl'
            f_rec = open(stored_rec, 'w')
            pickle.dump(recorder, f_rec)
            f_rec.close()
            del recorder
            
def load_out(params=None, file_path='.\mercia\\recorder_', nb_rep=15):
    out = {}
    if params is None:
        for i_sim in range(1,nb_rep):
            print i_sim
            stored_rec = file_path+str(i_sim)+'.pckl'
            f_rec = open(stored_rec)
            out[i_sim] = pickle.load(f_rec)
            f_rec.close()
    else:
        for par in params:
            out[par] = {}
            for i_sim in range(15):
                print par, i_sim
                stored_rec = file_path+par+'_'+str(i_sim)+'.pckl'
                f_rec = open(stored_rec)
                out[par][i_sim] = pickle.load(f_rec)
                f_rec.close()
    return out
    
def load_out_lat():
    latency = ['120', '220', '320']
    out = {}
    for lat in latency:
        out[lat] = {}
        for i_sim in range(15):
            print lat, i_sim
            stored_rec = '.\mercia\\recorder_lat_'+lat+'_'+str(i_sim)+'.pckl'
            f_rec = open(stored_rec)
            out[lat][i_sim] = pickle.load(f_rec)
            f_rec.close()
    return out
    
def save_canopies_size(nplants=[1, 3, 5, 10, 30]):
    for nb in nplants:
        print '---------------------------------------------------'
        print 'number of plants: %d' % nb
        print '---------------------------------------------------'
        make_canopy(nplants = nb, dir = './adel/adel_%d' % nb)

def save_canopies_date(years=[1998, 2003, 2004]):
    for yr in years:
        print '---------------------------------------------------'
        print 'Year: %d' % yr
        print '---------------------------------------------------'
        make_canopy(start_date = str(yr)+"-10-15 12:00:00", end_date = str(yr+1)+"-06-20 01:00:00",
                    nplants = 3, dir = './adel/adel_%d' % yr)
                    
def save_ids_date(years=[1998, 2003, 2004, 2010]):
    for yr in years:
        print '---------------------------------------------------'
        print 'Year: %d' % yr
        print '---------------------------------------------------'
        if yr!=2010:
            save_leaf_ids(start_date = str(yr)+"-10-15 12:00:00", end_date = str(yr+1)+"-06-20 01:00:00",
                        nplants = 3, dir = './adel/adel_%d' % yr)
        else:
            save_leaf_ids(start_date = str(yr)+"-10-15 12:00:00", end_date = str(yr+1)+"-06-20 01:00:00",
                        nplants = 3, dir = './adel/adel_saved')
                        
def stat_profiler(call='run_disease()'):
    import cProfile
    import pstats
    cProfile.run(call, 'restats')
    return pstats.Stats('restats')   