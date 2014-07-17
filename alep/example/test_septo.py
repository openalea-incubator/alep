from septo_imports import *
from alinea.alep.disease_outputs import plot_severity_by_leaf
import collections
import pickle
import sys

# Temp
from alinea.alep.disease_outputs import count_lesions
from alinea.alep.architecture import get_leaves
import matplotlib.pyplot as plt
plt.ion()

class DummyEmission():
    def __init__(self, **kwds):
        self.kwds = kwds
    def get_dispersal_units(self, g, fungus_name="septoria", label='LeafElement', weather_data=None):
        leaves = get_leaves(g)
        du = plugin_septoria().dispersal_unit(**self.kwds)
        du.position=[[10.,0.] for i in range(100)]
        return {lf:[du] for lf in leaves}

class DummyTransport():
    def __init__(self, **kwds):
        self.kwds = kwds
    def disperse(self, g, DU, weather_data=None):
        return DU

def run_simulation(start_year=1998, **kwds):

    # Initiation of the simulation ##########################################
    # Set the seed of the simulation
    # rd.seed(0)
    # np.random.seed(0)

    # Read weather and adapt it to septoria (add wetness)
    weather_file = 'meteo'+ str(start_year)[-2:] + '-' + str(start_year+1)[-2:] + '.txt'
    meteo_path = shared_data(alinea.septo3d, weather_file)
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    start_date = str(start_year)+"-10-01 01:00:00"
    weather.check(varnames=['degree_days'], models={'degree_days':basic_degree_days}, start_date=start_date)
    weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':basic_degree_days}, start_date=start_date, base_temp=-2.)
    seq = pandas.date_range(start = start_date,
                            end = str(start_year+1)+"-07-01 01:00:00", 
                            freq='H')
    print start_year
    
    # seq = pandas.date_range(start = str(start_year)+"-10-01 01:00:00",
                            # end = str(start_year)+"-12-31 01:00:00", 
                            # freq='H')
    # g, wheat, domain_area, domain = initialize_stand(age=1250., length=0.1, width=0.1,
        # sowing_density=150, plant_density=150, inter_row=0.12, nsect=5, seed=3)
                            
    # Initialize a wheat canopy
    nsect=5
    g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, width=0.1,
        sowing_density=150, plant_density=150, inter_row=0.12, nsect=nsect, seed=3)
        
    # Initialize the models for septoria
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
    septoria = plugin_septoria()
    # generator = DU_Generator(septoria, group_dus=True, **kwds)
    # generator = DU_Generator(septoria, group_dus=True, nb_rings_by_state=1, **kwds)
    # generator = DU_Generator(septoria, group_dus=False)
    if 'sporulating_fraction' in kwds:
        frac = kwds['sporulating_fraction']
        del kwds['sporulating_fraction']
    else:
        # See Baccar et al. for parameters
        # frac = 0.65e-4
        frac = 0.01
    septoria.parameters(group_dus=True, nb_rings_by_state=1, **kwds)
    inoc = SoilInoculum(septoria, sporulating_fraction=frac, domain_area=domain_area)
    contaminator = Septo3DSoilContamination(domain=domain, domain_area=domain_area)

    growth_controler = PriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    # sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
    # emitter = SeptoriaRainEmission(domain_area=domain_area)
    # transporter = Septo3DSplash()
    emitter = PopDropsEmission(domain=domain)
    transporter = PopDropsTransport(domain=domain, domain_area=domain_area)
    # emitter = DummyEmission(group_dus=True, **kwds)
    # transporter = DummyTransport()
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

    # Prepare saving of outputs
    recorders = {}
    ind=4
    for blade in [116, 128, 140]:
        ind -= 1
        recorders['F%d' % ind] = SeptoRecorder(vids=[blade+2+sect for sect in range(nsect)],group_dus=True)
    recorders['rosette'] = SeptoRecorder(vids=[8+2+sect for sect in range(nsect)],group_dus=True)
    
    # Simulation ############################################################
    for i, controls in enumerate(zip(weather_timing, wheat_timing, septo_timing, rain_timing)):
        weather_eval, wheat_eval, septo_eval, rain_eval = controls
        date = weather_eval.value.index[0]
        # if wheat_eval:
            # print date
        # dates.append(date)

        # Get weather for date and add it as properties on leaves
        if weather_eval:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = weather_eval.value.temperature_air,
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
            # Note : The position of senescence goes back to its initial value after
            # a while for undetermined reason
            # --> temporary hack for keeping senescence position low when it is over
            positions = g.property('position_senescence')
            greens = g.property('is_green')
            areas = g.property('area')
            senesced_areas = g.property('senesced_area')
            leaves = get_leaves(g, label = 'LeafElement')
            vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
            positions.update({vid:(0 if (positions[vid]==1 and not greens[vid]) or
                                        (positions[vid]>0 and round(areas[vid],5)==round(senesced_areas[vid],5))
                                        else positions[vid]) for vid in vids})

        # Update g for the disease:
        if septo_eval:
            #sen_model.find_senescent_lesions(g, label = 'LeafElement')
            update_healthy_area(g, label = 'LeafElement')
        
        # External contamination
        if rain_eval:        
            if rain_eval.value.rain.mean()>0. and rain_eval.value.degree_days[-1]<1000:
                g = external_contamination(g, inoc, contaminator, weather_eval.value)

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
        
        # Save outputs
        for recorder in recorders.itervalues():
            recorder.record(g, date)

    for recorder in recorders.itervalues():
        recorder.get_complete_dataframe()
        recorder.get_audpc()
        
    return g, recorders

def run_simulation_opti(start_year=1998, **kwds):
    # Initiation of the simulation ##########################################
    # Set the seed of the simulation
    # rd.seed(0)
    # np.random.seed(0)

    # Read weather and adapt it to septoria (add wetness)
    weather_file = 'meteo'+ str(start_year)[-2:] + '-' + str(start_year+1)[-2:] + '.txt'
    meteo_path = shared_data(alinea.septo3d, weather_file)
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    start_date = str(start_year)+"-10-01 01:00:00"
    weather.check(varnames=['degree_days'], models={'degree_days':basic_degree_days}, start_date=start_date)
    weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':basic_degree_days}, start_date=start_date, base_temp=-2.)
    seq = pandas.date_range(start = start_date, end = str(start_year+1)+"-07-01 01:00:00", freq='H')
    print start_year
    
    # seq = pandas.date_range(start = str(start_year)+"-10-01 01:00:00",
                            # end = str(start_year)+"-12-31 01:00:00", 
                            # freq='H')
    # g, wheat, domain_area, domain = initialize_stand(age=1250., length=0.1, width=0.1,
        # sowing_density=150, plant_density=150, inter_row=0.12, nsect=5, seed=3)
                            
    # Initialize a wheat canopy
    nsect=5
    g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, width=0.1,
        sowing_density=150, plant_density=150, inter_row=0.12, nsect=nsect, seed=3)
        
    # Initialize the models for septoria
    if 'alinea.alep.septoria_age_physio' in sys.modules:
        del(sys.modules['alinea.alep.septoria_age_physio'])
    septoria = plugin_septoria()
    if 'sporulating_fraction' in kwds:
        frac = kwds['sporulating_fraction']
        del kwds['sporulating_fraction']
    else:
        # See Baccar et al. for parameters
        # frac = 0.65e-4
        frac = 0.01
    septoria.parameters(group_dus=True, nb_rings_by_state=1, **kwds)
    inoc = SoilInoculum(septoria, sporulating_fraction=frac, domain_area=domain_area)
    contaminator = Septo3DSoilContamination(domain=domain, domain_area=domain_area)

    growth_controler = PriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    emitter = PopDropsEmission(domain=domain)
    transporter = PopDropsTransport(domain=domain, domain_area=domain_area)

    # Define the schedule of calls for each model
    every_h = time_filter(seq, delay=1)
    every_24h = time_filter(seq, delay=24)
    every_rain = rain_filter(seq, weather)
    septo_filter = septo_infection_filter(seq, weather, every_rain)
    weather_timing = IterWithDelays(*time_control(seq, every_h, weather.data))
    wheat_timing = IterWithDelays(*time_control(seq, every_24h, weather.data))
    septo_timing = CustomIterWithDelays(*time_control(seq, septo_filter, weather.data), eval_time='end')
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    
    # Prepare saving of outputs
    recorders = {}
    ind=4
    for blade in [116, 128, 140]:
        ind -= 1
        recorders['F%d' % ind] = SeptoRecorder(vids=[blade+2+sect for sect in range(nsect)],group_dus=True)
    recorders['rosette'] = SeptoRecorder(vids=[8+2+sect for sect in range(nsect)],group_dus=True)
    
    # Simulation ############################################################
    for i, controls in enumerate(zip(weather_timing, wheat_timing, septo_timing, rain_timing)):
        weather_eval, wheat_eval, septo_eval, rain_eval = controls
        date = weather_eval.value.index[0]
        if wheat_eval:
            print date
        # dates.append(date)

        # Get weather for date and add it as properties on leaves
        if septo_eval:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = septo_eval.value.temperature_air,
                           wetness = septo_eval.value.wetness.mean(),
                           relative_humidity = septo_eval.value.relative_humidity.mean())
        if rain_eval:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = rain_eval.value.rain.mean(),
                           rain_duration = len(rain_eval.value.rain) if rain_eval.value.rain.sum() > 0 else 0.)
        
        # Grow wheat canopy
        if wheat_eval:
            g,_ = grow_canopy(g, wheat, wheat_eval.value)
            # Note : The position of senescence goes back to its initial value after
            # a while for undetermined reason
            # --> temporary hack for keeping senescence position low when it is over
            positions = g.property('position_senescence')
            senesced_length = g.property('senesced_length')
            length = g.property('length')
            greens = g.property('is_green')
            areas = g.property('area')
            senesced_areas = g.property('senesced_area')
            leaves = get_leaves(g, label = 'LeafElement')
            vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
            positions.update({vid:(0 if (positions[vid]==1 and not greens[vid]) or
                                        (positions[vid]>0 and round(areas[vid],5)==round(senesced_areas[vid],5))
                                        else positions[vid]) for vid in vids})
            # senesced_length.update({vid:(length[vid] if (senesced_length[vid]==0. and not greens[vid]) or
                                        # (senesced_length[vid]<length[vid] and round(areas[vid],5)==round(senesced_areas[vid],5))
                                        # else senesced_length[vid]) for vid in vids})

        # Update g for the disease:
        if septo_eval:
            # print 'dt : %d' % septo_eval.dt
            update_healthy_area(g, label = 'LeafElement')
        
        # External contamination
        if rain_eval:        
            if rain_eval.value.rain.mean()>0. and rain_eval.value.degree_days[-1]<1000:
                g = external_contamination(g, inoc, contaminator, weather_eval.value)

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
        
        if wheat_eval:
            scene = plot_severity_by_leaf(g)
            # print 'nb lesions %d' % count_lesions(g)
        
        # Save outputs
        if septo_eval:
            for recorder in recorders.itervalues():
                recorder.record(g, date, degree_days = septo_eval.value.degree_days[-1])

    for recorder in recorders.itervalues():
        recorder.get_complete_dataframe()
        recorder.get_audpc()
        
    return g, recorders
    
def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)
    
def run_and_save(year, **kwds):
    g, recorder = run_simulation(start_year=year, **kwds)
    ind = 0
    if len(kwds)>0:
        stored_rec = '.\sensitivity\\recorder_'+str(year)+'_'+kwds.keys()[0][:3]+'-'+str(kwds.values()[0])+'.pckl'
    else:
        stored_rec = '.\sensitivity\\recorder_'+str(year)+'.pckl'
    f_rec = open(stored_rec, 'w')
    pickle.dump(recorder, f_rec)
    f_rec.close()
    del recorder

def run_stability(start_years = [1998, 2001, 2002], nb_rep=60):
    for i_sim in range(nb_rep):
        for year in start_years:
            # Compute and store results for simulation with variability
            g, recorder = run_simulation(year)
            stored_rec = '.\stability\\recorder_'+str(year)+'_'+str(40+i_sim)+'.pckl'
            f_rec = open(stored_rec, 'w')
            pickle.dump(recorder, f_rec)
            f_rec.close()
            del recorder

def load_stability(year = 1998, nb_rep=71):
    recs = {}
    for i_sim in range(nb_rep):
        print i_sim
        stored_rec = '.\stability\\recorder_'+str(year)+'_'+str(i_sim)+'.pckl'
        f_rec = open(stored_rec)
        recs[i_sim] = pickle.load(f_rec)
        f_rec.close()
    return recs
    
def test_stability(recs, nb_tested=15, nb_combi=50, nb_rep=71):
    fig, axs = plt.subplots(1,2)
    leaves = ['F1', 'F2', 'F3']
    colors = ['b', 'r', 'k']
    for i_leaf in range(len(leaves)):
        df_sev = pandas.DataFrame()
        df_nec = pandas.DataFrame()
        for k,reco in recs.iteritems():
            df_sev[k] = reco[leaves[i_leaf]].data['severity']
            df_nec[k] = reco[leaves[i_leaf]].data['necrosis_percentage']
        for i in range(nb_combi):
            df_sev_rand = df_sev[df_sev.columns[random.sample(range(nb_rep), nb_tested)]]
            axs[0].plot(df_sev_rand.index, df_sev_rand.mean(axis=1), color=colors[i_leaf])
            axs[0].set_ylabel('severity (in %)', fontsize = 20)
            axs[0].set_ylim([0.,1.])
            
            df_nec_rand = df_nec[df_nec.columns[random.sample(range(nb_rep), nb_tested)]]
            axs[1].plot(df_nec_rand.index, df_nec_rand.mean(axis=1), color=colors[i_leaf])
            axs[1].set_ylabel('necrotic percentage (in %)', fontsize = 20)
            axs[1].set_ylim([0.,1.])
    fig.text(0.5, 0.95, 'Means of '+ str(nb_tested)+' repetitions', fontsize = 20, horizontalalignment='center')

def test_stability2(recs, nb_tested=15, nb_rep=71):
    fig, axs = plt.subplots(1,2)
    leaves = ['F1', 'F2', 'F3']
    colors = ['b', 'r', 'k']
    for i_leaf in range(len(leaves)):
        df_sev = pandas.DataFrame()
        df_nec = pandas.DataFrame()
        for k,reco in recs.iteritems():
            df_sev[k] = reco[leaves[i_leaf]].data['severity']
            df_nec[k] = reco[leaves[i_leaf]].data['necrosis_percentage']
        for i in range(int(nb_rep/nb_tested)):
            df_sev_rand = df_sev[df_sev.columns[i*int(nb_rep/nb_tested):(i+1)*int(nb_rep/nb_tested)]]
            axs[0].plot(df_sev_rand.index, df_sev_rand.mean(axis=1), color=colors[i_leaf])
            axs[0].set_ylabel('severity (in %)', fontsize = 20)
            axs[0].set_ylim([0.,1.])
            
            df_nec_rand = df_nec[df_nec.columns[i*int(nb_rep/nb_tested):(i+1)*int(nb_rep/nb_tested)]]
            axs[1].plot(df_nec_rand.index, df_nec_rand.mean(axis=1), color=colors[i_leaf])
            axs[1].set_ylabel('necrotic percentage (in %)', fontsize = 20)
            axs[1].set_ylim([0.,1.])
    fig.text(0.5, 0.95, 'Means of '+ str(nb_tested)+' repetitions', fontsize = 20, horizontalalignment='center')

    
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
    run_all_simulations(start_years = [1998, 2001, 2002], nb_rings_by_state=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    run_all_simulations(start_years = [1998, 2001, 2002], frac=[0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5])

def run_sensi_smax():
    run_all_simulations(start_years = [1998, 2001, 2002], Smax=[0.04, 0.3, 3])