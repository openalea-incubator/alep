""" Example of use of the model of continuous lesions of septoria with real weather data.

"""
# Imports #########################################################################
import random as rd
import numpy as np
import pandas
from pylab import *
import matplotlib.pyplot as plt
import string

from alinea.alep.wheat_examples import adel_mtg, adel_mtg2, adel_mtg3, adel_one_leaf
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

from alinea.alep.protocol import *
from alinea.alep.septoria import compute_dispersal_event
from alinea.alep.septoria import *
from alinea.echap.microclimate_leaf import *
from alinea.alep.architecture import *
from alinea.alep.disease_operation import *
from alinea.alep.disease_outputs import *
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.du_position_checker import BiotrophDUProbaModel
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.dispersal import RandomDispersal
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.washing import RapillyWashing

from alinea.alep.alep_weather import *
from alinea.weather.global_weather import *
from alinea.weather.mini_models import leaf_wetness_rapilly
from datetime import datetime, timedelta
import time

# Plant ###########################################################################
def set_initial_properties_g(g):
    """ Basic setting for initial properties on wheat """
    set_properties(g,label = 'LeafElement',
               surface=5., healthy_surface=5., position_senescence=None)

# Climate #########################################################################
def read_weather(filename, update_events="step_by_step"):
    """ Use weather toolkit to read weather.
    
    update_events="step_by_step" or "dispersal"
    
    """
    weather = Weather(data_file = filename)
    weather = add_wetness(weather)
    weather = add_dispersal_events(weather, dispersal_event_model = compute_dispersal_event)  
    weather = add_calls_update(weather, update_events)  
    if update_events == "dispersal":
        weather = sum_temp_calls(weather)
    return weather

def add_calls_update(weather, update_events="step_by_step"):
    """ Add a column to weather data to indicate the hours when 
        the update of the lesion must be achieved.
        
    Either it is done at each hour ("step_by_step"),
    or it is done at each dispersal event ("dispersal").
    """
    if update_events == "dispersal":
        call = zeros(len(weather.data))
    elif update_events == "step_by_step":
        call = ones(len(weather.data))
    else:
        raise Exception("Wrong str for update_events")
    
    if update_events == "dispersal":
        for i_line in range(len(weather.data)):
            if (weather.data.rain[i_line] > 0. and
                weather.data.dispersal_event[i_line]==True):
                call[i_line] = True
        
    update_events = pandas.DataFrame(dict(update_events=call))
    update_events.index = weather.data.temperature_air.index
    weather.data = weather.data.join(update_events)
    return weather

def sum_temp_calls(weather):
    """ Replace temperature by sum of temperature between calls.
    """
    # Note : This method will be problematic for hours when temp<0
    temp = zeros(len(weather.data))
    sum_temp = 0.
    for i_line in range(len(weather.data)):
        sum_temp += weather.data.temperature_air[i_line]
        if weather.data.update_events[i_line]:
            temp[i_line] = sum_temp
            sum_temp = 0.
    
    del weather.data['temperature_air']
    temperature_air = pandas.DataFrame(dict(temperature_air=temp))
    temperature_air.index = weather.data.update_events.index
    weather.data = weather.data.join(temperature_air)
    return weather

def compute_hourly_delta_ddays(temp=0., basis_for_dday=-2):
    """ Compute delta degree days in the hour.
    
    Parameters
    ----------
    dt: int
        Time step of the simulation (in hours)
    temp: float
        Temperature from weather data
    basis_for_dday: float
        Basis temperature for degree days calculation
    """
   
    return max(0,(temp - basis_for_dday)/24.)
    
# Fungus ##########################################################################
def execute_one_step(g, weather, date, dt,
                    position_checker=BiotrophDUProbaModel(),
                    controler=NoPriorityGrowthControl(),
                    dispersor=RandomDispersal(), 
                    washor=RapillyWashing()):
    """ Execute one time step of complete simulation for septoria.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    weather: object weather
        Object for meteorological data (see 'read weather')
    date: format datetime
        Date of the simulation
    dt: int
        Time step of the simulation (in hours)
    position_checker: model
        Model that disable the DU if it falls on an existing lesion or senescent tissue
        Requires a method 'check position' (see doc)
    controler: model
        Model with rules of competition between the lesions
    dispersor: model
        Model used to position each DU in stock on g
    washor: model
        Model used to wash the DUs out of the leaf
       
    Returns
    -------
        Number of dus on the MTG
    """
    # Get weather for date advance date for next simulation step
    mgc, globalclimate = weather.get_weather(dt, date)
    date = weather.next_date(dt, date)
    
    set_properties(g,label = 'LeafElement',
                    wetness=globalclimate.wetness.values[0],
                    temp=globalclimate.temperature_air.values[0],
                    rain_intensity=globalclimate.rain.values[0],
                    rain_duration=globalclimate.rain_duration.values[0],
                    relative_humidity=globalclimate.relative_humidity.values[0],
                    wind_speed=globalclimate.wind_speed.values[0])  
    
    # Run a disease cycle
    # grow(g)
    infect(g, dt, position_checker)
    if globalclimate.update_events.values[0]:
        update(g, dt, growth_control_model=controler)
           
    global_rain_intensity = globalclimate.rain.values[0]
    if global_rain_intensity != 0. :
        seed(1)
        disperse(g, dispersor, "septoria")
        seed(1)
        wash(g, washor, global_rain_intensity)
    
    return g

# Useful function ##################################################################
def int2str(integer):
    """ Convert an integer to a string.
    """
    return "%d" % integer
    
# septoria ########################################################################
def test_microclimate():
    """
    """
    # Read weather data : 
    filename = 'meteo01.csv'
    weather_data = read_weather(filename, "step_by_step")
    # Generate a MTG with required properties :
    g = adel_one_leaf()
    # g = adel_mtg2()
    # Loop of simulation
    start_date = datetime(2000, 10, 1, 01, 00)
    end_date = datetime(2001, 7, 31)
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        # Update climate on leaf elements : 
        update_on_leaves(weather_data, date, g)
    
def test_all(model="septoria_exchanging_rings"):
    # Generate a MTG with required properties :
    g = adel_mtg2()
    # g = adel_one_leaf()
    set_initial_properties_g(g)
    
    # Generate a stock of septoria dispersal units and distribute it on g
    nb_dus_in_stock = 100
    distribute_dispersal_units(g, nb_dus=nb_dus_in_stock, disease_model=model)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

    # Prepare the simulation loop
    dt = 1
    nb_steps = 1000
    severity = np.zeros(nb_steps)
    nb_DUs = np.zeros(nb_steps)
    nb_lesions = np.zeros(nb_steps)
    for i in range(0,nb_steps,dt):
        print('time step %d' % i)
        dt_pluie = 1
        # Update climate and force rain occurences       
        if i>0 and i%dt_pluie == 0:
            global_rain_intensity = 4.
        else:
            global_rain_intensity = 0.
        set_properties(g, wetness=True, temp=22., 
                          rain_intensity = global_rain_intensity*0.75, 
                          rain_duration = 2.,
                          relative_humidity = 85.)
        
        # grow(g)
        infect(g, dt)
        update(g, dt_pluie, growth_control_model=controler)

        if global_rain_intensity != 0.:
            disperse(g, dispersor, "septoria")
            wash(g, washor, global_rain_intensity)
        
        # Measure severity
        severity[i]=compute_total_severity(g)
        
        nb_DUs[i] = count_dispersal_units(g)
        nb_lesions[i] = count_lesions(g)
    
    # Display results
    fig = plt.figure()
    
    ax1 = fig.add_subplot(3,1,1)
    plot(range(0,nb_steps,dt), severity)
    ylim([0, 105])
    ylabel('Total disease severity')
    
    ax2 = fig.add_subplot(3,1,2)
    plt.bar(range(0,nb_steps,dt), nb_DUs, color='r')
    ylabel('Number of DUs')

    ax3 = fig.add_subplot(3,1,3)
    plot(range(0,nb_steps,dt), nb_lesions, color='r')
    ylabel('Number of lesions')
    xlabel('Simulation time step')
    
    fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g

def compare_models():
    """
    """
    start_date = datetime(2000, 10, 1, 01, 00)
    end_date=datetime(2000, 11, 30)
    nb_lesions = 100
    # Run the simulation with each model
    severity1, nb_DUs1, nb_lesions1 = simulate_severity(start_date,end_date, nb_lesions,
                                                        model="septoria_exchanging_rings")
    severity2, nb_DUs2, nb_lesions2 = simulate_severity(start_date,end_date, nb_lesions,
                                                        model="septoria_with_rings")
    severity3, nb_DUs3, nb_lesions3 = simulate_severity(start_date,end_date, nb_lesions,
                                                        model="septoria_continuous")
    
    # Display results
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(3,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), severity1, color='b')
    plot(pandas.date_range(start_date, end_date, freq='H'), severity2, color='r')
    plot(pandas.date_range(start_date, end_date, freq='H'), severity3, color='g')
    ylim([0, 105])
    ylabel('Total disease severity')
    
    ax2 = fig.add_subplot(3,1,2)
    plt.bar(range(nb_steps), nb_DUs1, color='b')
    plt.bar(range(nb_steps), nb_DUs2, color='r')
    plt.bar(range(nb_steps), nb_DUs3, color='g')
    ylabel('Number of DUs')
    xlim([0, nb_steps])

    ax3 = fig.add_subplot(3,1,3)
    plot(range(nb_steps), nb_lesions1, color='b')
    plot(range(nb_steps), nb_lesions2, color='r')
    plot(range(nb_steps), nb_lesions3, color='g')
    ylabel('Number of lesions')
    xlim([0, nb_steps])
    xlabel('Simulation time step')
    
    fig.subplots_adjust(hspace=1)
    
    plt.show()

def compare_time_steps(model="septoria_continuous"):
    """
    """
    start_date = datetime(2000, 10, 1, 01, 00)
    # end_date=datetime(2000, 11, 30)
    end_date=datetime(2001, 07, 31)
    nb_lesions = 1
    # Run the simulation with each model
    severity1, nb_DUs1, nb_lesions1 = simulate_severity(start_date, end_date, nb_lesions,
                                                        model=model,
                                                        update_events="step_by_step")                                                        
    severity2, nb_DUs2, nb_lesions2 = simulate_severity(start_date, end_date, nb_lesions,
                                                        model=model,
                                                        update_events = "dispersal")
    
    # Display results
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(4,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), severity1, color='b')
    plot(pandas.date_range(start_date, end_date, freq='H'), severity2, color='r')
    ylim([0, 105])
    ylabel('Total disease severity')
    
    ax2 = fig.add_subplot(4,1,2)
    plt.bar(range(nb_steps), nb_DUs1, color='b')
    ylabel('Number of DUs')
    xlim([0, nb_steps])
    ax3 = fig.add_subplot(4,1,3)
    plt.bar(range(nb_steps), nb_DUs2, color='r')
    ylabel('Number of DUs')
    xlim([0, nb_steps])

    ax4 = fig.add_subplot(4,1,4)
    plot(range(nb_steps), nb_lesions1, color='b')
    plot(range(nb_steps), nb_lesions2, color='r')
    ylabel('Number of lesions')
    xlim([0, nb_steps])
    xlabel('Simulation time step')
    
    fig.subplots_adjust(hspace=1)
    
    plt.show()
    
def simulate_severity(start_date=datetime(2000, 10, 1, 01, 00),
                      end_date=datetime(2000, 11, 30),
                      nb_lesions = 100,
                      model="septoria_exchanging_rings",
                      update_events = "step_by_step"):
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.csv'
    weather = read_weather(filename, update_events)
    
    # Generate a MTG with required properties :
    g = adel_one_leaf()
    # g = adel_mtg2()
    set_initial_properties_g(g)
    
    # Generate a stock of septoria lesions and distribute them on g
    nb_lesions_in_stock = nb_lesions
    # distribute_dispersal_units(g, nb_dus=nb_lesions_in_stock, disease_model=model)
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, disease_model=model)

    # Prepare the simulation loop
    dt = 1
    severity = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    nb_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    nb_lesions = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try: ind+=1.
        except: ind = 0.

        g = execute_one_step(g, weather, date, dt)
        
        # Measure severity
        severity[ind]=compute_total_severity(g)
        
        nb_DUs[ind] = count_dispersal_units(g)
        nb_lesions[ind] = count_lesions(g)
    
    return severity, nb_DUs, nb_lesions

def imitate_septo3D(start_date=datetime(2000, 10, 1, 01, 00),
                      end_date=datetime(2000, 11, 30),
                      nb_lesions = 100,
                      model="septoria_exchanging_rings",
                      update_events = "step_by_step"):
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.csv'
    weather = read_weather(filename, update_events)
    
    # Generate a MTG with required properties :
    g = adel_one_leaf()
    # g = adel_mtg2()
    set_initial_properties_g(g, surface_leaf_element=10.)
    
    # Generate a stock of septoria lesions and distribute them on g
    nb_lesions_in_stock = nb_lesions
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, disease_model=model)
    
    # Call the models that will be used during the simulation :
    position_checker=BiotrophDUProbaModel()
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()
    leaf_inspector1 = LeafInspector(g, leaf_number=1)
    # leaf_inspector2 = LeafInspector(g, leaf_number=2)
    # leaf_inspector3 = LeafInspector(g, leaf_number=3)
    # leaf_inspector4 = LeafInspector(g, leaf_number=4)

    # Prepare the simulation loop
    dt = 1
    degree_days = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    healthy_surface = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    total_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    unwashed_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    on_green_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    incubating_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
       
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try: 
            ind+=1.
            # Estimate delta degree days
            degree_days[ind] = degree_days[ind-1] + compute_hourly_delta_ddays(temp=weather.date.temperature_air[date])
        except: 
            ind = 0.
            # Estimate delta degree days
            degree_days[ind] = compute_hourly_delta_ddays(temp=weather.date.temperature_air[date])
        
        # Get weather for date advance date for next simulation step
        mgc, globalclimate = weather.get_weather(dt, date)
        date = weather.next_date(dt, date)
        
        set_properties(g,label = 'LeafElement',
                        wetness=globalclimate.wetness.values[0],
                        temp=globalclimate.temperature_air.values[0],
                        rain_intensity=globalclimate.rain.values[0],
                        rain_duration=globalclimate.rain_duration.values[0],
                        relative_humidity=globalclimate.relative_humidity.values[0],
                        wind_speed=globalclimate.wind_speed.values[0])
        
        # Run a disease cycle
        infect(g, dt, position_checker)
        if globalclimate.update_events.values[0]:
            update(g, dt, growth_control_model=controler)

        # Healthy surface
        healthy_surface[ind] = leaf_inspector1.compute_healthy_surface(g)
        
        # Surfaces in state
        leaf_inspector1.compute_states(g)
        # leaf_inspector2.compute_states(g)
        # leaf_inspector3.compute_states(g)
        # leaf_inspector4.compute_states(g)

        global_rain_intensity = globalclimate.rain.values[0]
        if global_rain_intensity != 0. :
            initial_DUs = leaf_inspector1.count_dispersal_units(g)
            scene = plot3d(g)
            seed(1)
            disperse(g, dispersor, "septoria")
            seed(1)
            total_DUs[ind] = leaf_inspector1.count_dispersal_units(g) - initial_DUs
            wash(g, washor, global_rain_intensity)
            unwashed_DUs[ind] = leaf_inspector1.count_dispersal_units(g) - initial_DUs
            on_green_DUs[ind] = leaf_inspector1.count_DU_on_green(g, unwashed_DUs[ind])
            # print('pluie')

    # Draw the outputs:
    # 1. Build the frame
    fig = plt.figure()
    nb_graphs = 10
    length = 5
    width = nb_graphs/length
    all_letters = string.ascii_uppercase
    for i in range(nb_graphs):
        graph_handle = 'ax'+int2str(i+1)
        # if i < 4:
            # globals()[graph_handle] = fig.add_subplot(length,width,i+1,
                                                      # xticklabels=[], yscale='log')
        # elif i < 8:
        if i < 8:
            globals()[graph_handle] = fig.add_subplot(length,width,i+1, xticklabels=[])
        else:
            globals()[graph_handle] = fig.add_subplot(length,width,i+1)
        graph_letter = all_letters[i]
        globals()[graph_letter] = plt.text(0.05, 0.9, graph_letter, 
                                        fontsize=18, ha='center',
                                        va='center', transform=globals()[graph_handle].transAxes)

    # Titles
    droplets = plt.text(-0.275, 2.1, 'Number of Infectious Droplets', 
                        fontsize=18, rotation='vertical', transform=ax3.transAxes)
    leaf_area = plt.text(-0.275, 1.4, 'Leaf Area', 
                        fontsize=18, rotation='vertical', transform=ax7.transAxes)
    coverage = plt.text(-0.275, 0.75, '% Recovered', 
                        fontsize=18, rotation='vertical', transform=ax9.transAxes)
        
    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    # ax1.bar(pandas.date_range(start_date, end_date, freq='H'), total_DUs,  color='b')
    # ax1.set_yscale('log')
    ax1.bar(degree_days, total_DUs,  color='b')
    # ax1.plot(pandas.date_range(start_date, end_date, freq='H'), washing_rate)   
    # ax1.plot(degree_days, washing_rate)   
    ax1.set_ylabel('Total intercepted')
    ax1.set_xlim([0, max(degree_days)])
    ax1.set_ylim([0, max(total_DUs)+20])
    
    ax2.bar(degree_days, unwashed_DUs,  color='b')
    ax2.set_xlim([0, max(degree_days)])
    ax2.set_ylim([0, max(total_DUs)+20])
    ax2.set_ylabel('Total unwashed')
    
    ax3.bar(degree_days, on_green_DUs,  color='b')
    ax3.set_xlim([0, max(degree_days)])
    ax3.set_ylim([0, max(total_DUs)+20])
    ax3.set_ylabel('Total on Green')
    
    ax33 = ax3.twinx()
    ax33.plot(degree_days, healthy_surface,  color='b', linestyle='--')
    ax33.set_xticklabels([])
    ax33.set_xlim([0, max(degree_days)])
    
    ax5.plot(degree_days, leaf_inspector1.ratio_inc, color='b')
    ax5.set_ylim([0, 105])
    ax5.set_xlim([0, max(degree_days)])
    
    ax6.plot(degree_days, leaf_inspector1.ratio_chlo, color='b')
    ax6.set_ylim([0, 105])
    ax6.set_xlim([0, max(degree_days)])
    
    ax7.plot(degree_days, leaf_inspector1.ratio_nec, color='b')
    ax7.set_ylim([0, 105])
    ax7.set_xlim([0, max(degree_days)])
    
    ax8.plot(degree_days, leaf_inspector1.ratio_spo, color='b')
    ax8.set_ylim([0, 105])
    ax8.set_xlim([0, max(degree_days)])
    
    # raise Exception('')
    
    # ax2 = fig.add_subplot(3,1,2)
    # plot(pandas.date_range(start_date, end_date, freq='H'), weather_data.Pluie[start_date:end_date])
    # ylabel('Rain intensity')
    
    # ax3 = fig.add_subplot(3,1,3)
    # plot(pandas.date_range(start_date, end_date, freq='H'), weather_data.rain_duration[start_date:end_date])
    # ylabel('Rain duration')
    
    # ax2 = fig.add_subplot(1,2,2)
    # plt.bar(pandas.date_range(start_date, end_date, freq='H'), unwashed_DUs, color='b')
    # ylabel('Total unwashed')
    # ax2.set_yscale('log')
    # xlim([0, nb_steps])

    # ax3 = fig.add_subplot(3,1,3)
    # plot(range(nb_steps), nb_lesions, color='r')
    # ylabel('Number of lesions')
    # xlim([0, nb_steps])
    # xlabel('Simulation time step')
    
    # fig.subplots_adjust(hspace=1)
    
    plt.show(False)

def simulate_necrosis():
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.csv'
    weather = read_weather(filename, "step_by_step")
    
    # Generate a MTG with required properties :
    # g = adel_one_leaf()
    g = adel_mtg2()
    set_initial_properties_g(g)
    
    # Generate a stock of septoria lesions and distribute them on g
    nb_lesions_in_stock = 100
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, disease_model='septoria_exchanging_rings')
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

    # Prepare the simulation loop
    dt = 1
    start_date = datetime(2000, 10, 1, 01, 00)
    end_date = datetime(2001, 7, 31)
    # end_date = datetime(2000, 11, 30)
    necrosis = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try:
            ind+=1.
        except:
            ind = 0.

        # Get weather for date advance date for next simulation step
        mgc, globalclimate = weather.get_weather(dt, date)
        date = weather.next_date(dt, date)
        
        set_properties(g,label = 'LeafElement',
                        wetness=globalclimate.wetness.values[0],
                        temp=globalclimate.temperature_air.values[0],
                        rain_intensity=globalclimate.rain.values[0],
                        rain_duration=globalclimate.rain_duration.values[0],
                        relative_humidity=globalclimate.relative_humidity.values[0],
                        wind_speed=globalclimate.wind_speed.values[0])
        
        # Run a disease cycle
        #grow(g)
        infect(g, dt)
        if globalclimate.update_events.values[0]:
            update(g, dt, growth_control_model=controler)
               
        global_rain_intensity = globalclimate.rain.values[0]
        if global_rain_intensity != 0. :
            scene = plot3d(g)
            disperse(g, scene, dispersor, "septoria") 
            wash(g, washor, global_rain_intensity)
    
        # Measure severity
        necrosis[ind]=compute_total_necrosis(g)
    
    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    # fig = plt.figure()
    
    # ax1 = fig.add_subplot(3,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), necrosis)
    ylim([0, 105])
    ylabel('Total necrosis percentage')
    
    # ax2 = fig.add_subplot(3,1,2)
    # plt.bar(range(nb_steps), nb_DUs, color='r')
    # ylabel('Number of DUs')
    # xlim([0, nb_steps])

    # ax3 = fig.add_subplot(3,1,3)
    # plot(range(nb_steps), nb_lesions, color='r')
    # ylabel('Number of lesions')
    # xlim([0, nb_steps])
    # xlabel('Simulation time step')
    
    # fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g
    
def simulate_states_one_leaf(model="septoria_exchanging_rings"):
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.csv'
    weather = read_weather(filename, "step_by_step")
    
    # Generate a MTG with required properties :
    g = adel_one_leaf()
    set_initial_properties_g(g)
    
    # Generate a stock of septoria lesions and distribute them on g
    nb_lesions_in_stock = 10
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, disease_model=model)
    
    # Call the models that will be used during the simulation :
    position_checker = BiotrophDUProbaModel()
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = Septo3DSplash()
    
    # Call the variable saver
    outputs = LeafInspector(g, leaf_number=1)
    
    # Prepare the simulation loop
    dt = 1
    start_date = datetime(2000, 10, 1, 01, 00)
    end_date = datetime(2001, 2, 1)
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try:
            ind+=1.
        except:
            ind = 0.
        # Get weather for date advance date for next simulation step
        mgc, globalclimate = weather.get_weather(dt, date)
        date = weather.next_date(dt, date)
        
        set_properties(g,label = 'LeafElement',
                        wetness=globalclimate.wetness.values[0],
                        temp=globalclimate.temperature_air.values[0],
                        rain_intensity=globalclimate.rain.values[0],
                        rain_duration=globalclimate.rain_duration.values[0],
                        relative_humidity=globalclimate.relative_humidity.values[0],
                        wind_speed=globalclimate.wind_speed.values[0])   
        
        # Run a disease cycle
        #grow(g)
        infect(g, dt) # , position_checker)
        if globalclimate.update_events.values[0]:
            update(g, dt, growth_control_model=controler)
                
        global_rain_intensity = globalclimate.rain.values[0]
        if global_rain_intensity != 0. :
            scene = plot3d(g)
            disperse(g, scene, dispersor, "septoria") 
            wash(g, washor, global_rain_intensity)
    
        # Save variables
        outputs.compute_severity(g)
        # outputs.compute_states(g)

    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    # y_max = max([max(incubating), max(chlorotic), max(necrotic), max(sporulating)])+1
    
    # fig = plt.figure()
    
    # ax1 = fig.add_subplot(4,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), outputs.severity)
    # plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_inc)
    # ylim([0, y_max])
    # ylabel('Incubating surface')
    
    # ax2 = fig.add_subplot(4,1,2)
    # plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_chlo, color='g')
    # ylabel('Chlorotic surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])

    # ax3 = fig.add_subplot(4,1,3)
    # plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_nec, color='r')
    # ylabel('Necrotic surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])

    # ax4 = fig.add_subplot(4,1,4)
    # plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_spo, color='k')
    # plot(pandas.date_range(start_date, end_date, freq='H'), created_surface, color='b', linestyle='--')
    # plot(pandas.date_range(start_date, end_date, freq='H'), sum_surface, color='k', linestyle='--')
    # ylabel('Sporulating surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])
    # xlabel('Simulation time step')
    
    # fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g
