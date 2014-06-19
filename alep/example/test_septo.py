from septo_imports import *
from alinea.alep.disease_outputs import plot_severity_by_leaf
import collections
import pickle
import sys

# Temp
from alinea.alep.disease_outputs import count_lesions

def run_simulation(start_year=1998, **kwds):

    # Initiation of the simulation ##########################################
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Read weather and adapt it to septoria (add wetness)
    weather_file = 'meteo'+ str(start_year)[-2:] + '-' + str(start_year+1)[-2:] + '.txt'
    meteo_path = shared_data(alinea.septo3d, weather_file)
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    weather.check(varnames=['degree_days'], models={'degree_days':basic_degree_days}, start_date=str(start_year)+"-10-01 01:00:00")
    seq = pandas.date_range(start = str(start_year)+"-10-01 01:00:00",
                            end = str(start_year+1)+"-07-01 01:00:00", 
                            freq='H')
    print start_year
    
    # seq = pandas.date_range(start = str(start_year)+"-10-01 01:00:00",
                            # end = str(start_year)+"-11-10 01:00:00", 
                            # freq='H')
    # g, wheat, domain_area, domain = initialize_stand(age=1250., length=0.1, width=0.1,
        # sowing_density=150, plant_density=150, inter_row=0.12, nsect=5, seed=3)
                            
    # Initialize a wheat canopy
    g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, width=0.1,
        sowing_density=150, plant_density=150, inter_row=0.12, nsect=5, seed=3)
        
    # Initialize the models for septoria
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
    septoria = plugin_septoria()
    generator = DU_Generator(septoria, group_dus=True, **kwds)
    # generator = DU_Generator(septoria, group_dus=True, nb_rings_by_state=1, **kwds)
    # generator = DU_Generator(septoria, group_dus=False)
    if 'sporulating_fraction' in kwds:
        frac = kwds['sporulating_fraction']
        del kwds['sporulating_fraction']
    else:
        # See Baccar et al. for parameters
        # frac = 0.65e-4
        frac = 0.01
    inoc = SoilInoculum(DU_generator=generator, sporulating_fraction=frac, domain_area=domain_area)
    contaminator = Septo3DSoilContamination(domain=domain, domain_area=domain_area)

    growth_controler = PriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    # sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
    # emitter = SeptoriaRainEmission(domain_area=domain_area)
    # transporter = Septo3DSplash()
    emitter = PopDropsEmission(domain=domain)
    transporter = PopDropsTransport(domain=domain, domain_area=domain_area)
    # transporter = Septo3DTransport(wash=True)
    # washor = RapillyWashing()
    
    # Temp
    inoculator = InoculationLowerLeaves()

    # Define the schedule of calls for each model
    every_h = time_filter(seq, delay=1)
    every_24h = time_filter(seq, delay=24)
    every_rain = rain_filter(seq, weather)
    weather_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    wheat_timing = IterWithDelays(*time_control(seq, every_24h, weather.data))
    septo_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))

    # Call leaf inspectors for target blades (top 3)
    inspectors = {}
    #first_blade = 80
    ind = 4.
    # for blade in range(first_blade,104,8):
    for blade in [116, 128, 140]:
            ind -= 1
            inspectors['F %d' % ind] = LeafInspector(g, blade_id=blade)    
    inspectors['rosette'] = LeafInspector(g, blade_id=8)
    dates = []
    
    # Simulation ############################################################
    for i, controls in enumerate(zip(weather_timing, wheat_timing, septo_timing, rain_timing)):
        weather_eval, wheat_eval, septo_eval, rain_eval = controls
        date = weather_eval.value.index[0]
        if wheat_eval:
            print date
        # dates.append(date)

        # Get weather for date and add it as properties on leaves
        if weather_eval:
            set_properties(g,label = 'LeafElement',
                           temp = weather_eval.value.temperature_air[0],
                           wetness = weather_eval.value.wetness[0],
                           relative_humidity = weather_eval.value.relative_humidity[0],
                           wind_speed = weather_eval.value.wind_speed[0])
        if rain_eval:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = rain_eval.value.rain.mean(),
                           rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)
        
        # Grow wheat canopy
        if wheat_eval:
            g,_ = grow_canopy(g, wheat, wheat_eval.value)
            # # Note : The position of senescence goes back to its initial value after
            # # a while for undetermined reason
            # # --> temporary hack for keeping senescence position low when it is over
            # positions = g.property('position_senescence')
            # greens = g.property('is_green')
            # areas = g.property('area')
            # senesced_areas = g.property('senesced_area')
            # leaves = get_leaves(g, label = 'LeafElement')
            # vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
            # # for vid in vids:
                # # if ('dispersal_units' in g.node(vid).properties() and len(g.node(vid).dispersal_units)>0. and
                    # # ((positions[vid]==1 and not greens[vid]) or
                   # # (positions[vid]>0 and round(areas[vid],5)==round(senesced_areas[vid],5)))):
                    # # import pdb
                    # # pdb.set_trace()
            # positions.update({vid:(0 if (positions[vid]==1 and not greens[vid]) or
                                        # (positions[vid]>0 and round(areas[vid],5)==round(senesced_areas[vid],5))
                                        # else positions[vid]) for vid in vids})


        # Update g for the disease:
        if septo_eval:
            #sen_model.find_senescent_lesions(g, label = 'LeafElement')
            update_healthy_area(g, label = 'LeafElement')
        
        # External contamination
        if rain_eval:        
            if rain_eval.value.rain.mean()>0. and rain_eval.value.degree_days[-1]<1000:
                g = external_contamination(g, inoc, contaminator, weather_eval.value)
                # dus = generate_stock_du(10, septoria, **kwds)
                # initiate(g, dus, inoculator)

        # Develop disease
        if septo_eval:
            # Update dispersal units and lesions
            infect(g, septo_eval.dt, infection_controler, label='LeafElement')
            update(g, septo_eval.dt, growth_controler, senescence_model=None, label='LeafElement')
            
        # Disperse and wash
        if rain_eval:
            if rain_eval.value.rain.mean()>0.:
                g = disperse(g, emitter, transporter, "septoria", label='LeafElement', weather_data=weather_eval.value)
                # wash(g, washor, rain_eval.value.rain.mean(), label='LeafElement')
        
        # if wheat_eval:
            # scene = plot_severity_by_leaf(g)
            # print 'nb lesions %d' % count_lesions(g)
        
        # # Save outputs
        for inspector in inspectors.itervalues():
            inspector.update_variables(g)
            # inspector.update_du_variables(g)

    for inspector in inspectors.itervalues():
        inspector.update_audpc()
        inspector.dates = dates
        
    return g, inspectors

def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)
    
def run_and_save(year, **kwds):
    g, inspector = run_simulation(start_year=year, **kwds)
    ind = 0
    if len(kwds)>0:
        stored_insp = '.\sensitivity\inspector_'+str(year)+'_'+kwds.keys()[0][:3]+'-'+str(kwds.values()[0])+'.pckl'
    else:
        stored_insp = '.\sensitivity\inspector_'+str(year)+'.pckl'
    f_insp = open(stored_insp, 'w')
    pickle.dump(inspector, f_insp)
    f_insp.close()
    del inspector

def load_out(year, **kwds):
    if len(kwds)>0:
        stored_insp =  '.\sensitivity\inspector_'+str(year)+'_'+kwds.keys()[0][:3]+'-'+str(kwds.values()[0])+'.pckl'
    else:
        stored_insp = '.\sensitivity\inspector_'+str(year)+'.pckl'
    f_insp = open(stored_insp)
    out = pickle.load(f_insp)
    f_insp.close()
    return out
    
def run_all_simulations(start_years = [1998, 2001, 2002], **kwds):
    # Run simulation for each year with variable and fixed chlorosis threshold
    for year in start_years:
        if len(kwds)>0:
            for k in kwds.keys():
                if not is_iterable(kwds[k]):
                    kwds[k] = [kwds[k]]
                for v in kwds[k]:
                    # Compute and store results for simulation with variability
                    run_and_save(year, **{k:v})
        else:
            run_and_save(year)
            
def read_outputs(start_years = [1998, 2001, 2002], **kwds):
    out = {}
    for year in start_years:
        if len(kwds)>0:
            for k in kwds.keys():
                if not is_iterable(kwds[k]):
                    kwds[k] = [kwds[k]]
                for v in kwds[k]:
                    out[str(year)+'_'+k[:3]+'-'+str(v)] = load_out(year, **{k:v})

        else:
            out[str(year)] = load_out(year)
    return out
    
def run_sensi():
    run_all_simulations(start_years = [1998, 2001, 2002], degree_days_to_chlorosis=[120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320])
    run_all_simulations(start_years = [1998, 2001, 2002], Smin=[0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3])
    run_all_simulations(start_years = [1998, 2001, 2002], Smax=[0.04, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 3, 5, 10])
    run_all_simulations(start_years = [1998, 2001, 2002], frac=[0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1])

def run_sensi_smax():
    run_all_simulations(start_years = [1998, 2001, 2002], Smax=[0.04, 0.3, 3])