""" Draw the outputs to verify the well functioning of dispersal units. """

# Useful imports
import random as rd
import numpy as np
import pandas
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator, DateFormatter
from alinea.astk.TimeControl import *
from alinea.alep.alep_color import alep_colormap, green_yellow_red

# Imports for weather
from datetime import datetime
from alinea.alep.alep_weather import get_septoria_weather

# Imports for wheat
from alinea.alep.wheat import initialize_stand, find_blade_id, find_leaf_ids
from alinea.astk.plant_interface import grow_canopy
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.disease_outputs import LeafInspector
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence

# Imports for powdery_mildew
from openalea.vpltk import plugin
from alinea.alep.dispersal_emission import PowderyMildewWindEmission
from alinea.alep.dispersal import PowderyMildewWindDispersal

# Useful functions ############################################################
def update_plot(g):
    from alinea.alep.architecture import set_property_on_each_id
    from alinea.alep.disease_outputs import compute_severity_by_leaf
    from alinea.adel.mtg_interpreter import plot3d
    from openalea.plantgl.all import Viewer

    # Compute severity by leaf
    severity_by_leaf = compute_severity_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'severity', severity_by_leaf, label = 'LeafElement')

    # Visualization
    g = alep_colormap(g, 'severity', cmap=green_yellow_red(levels=100),
                      lognorm=False, zero_to_one=False, vmax=100)

    leaves = get_leaves(g, label='LeafElement')
    pos_sen = g.property('position_senescence')
    for leaf in leaves:
        if pos_sen[leaf]==0.:
            g.node(leaf).color = (157, 72, 7)

    scene = plot3d(g)
    Viewer.display(scene)
    return scene

def run_simulation_septoria():
    # Initialization #####################################################
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Choose dates of simulation and initialize the value of date
    start_date = datetime(2000, 10, 1, 1, 00, 00)
    end_date = datetime(2001, 07, 01, 00, 00)
    date = None

    # Read weather and adapt it to septoria (add wetness)
    weather = get_septoria_weather(data_file='meteo01.csv')

    # Initialize a wheat canopy
    g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, 
                                                     width=0.2, sowing_density=150,
                                                     plant_density=150, inter_row=0.12)

    # Initialize the models for septoria
    septoria = plugin_septoria()
    inoculator = RandomInoculation()
    growth_controler = NoPriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
    emitter = SeptoriaRainEmission()
    transporter = Septo3DSplash(reference_surface=domain_area)
    washor = RapillyWashing()

    # Define the schedule of calls for each model
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    weather_timing = TimeControl(delay=1, steps=nb_steps)
    wheat_timing = TimeControl(delay=24, steps=nb_steps, model=wheat, weather=weather, start_date=start_date)
    septo_timing = TimeControl(delay=1, steps=nb_steps)
    timer = TimeControler(weather=weather_timing, wheat=wheat_timing, disease = septo_timing)

    # Call leaf inspectors for target blades (top 3)
    inspectors = {}
    # for rank in range(1,3):
        # inspectors[rank] = LeafInspector(g, blade_id=find_blade_id(g, leaf_rank = rank, only_visible=False))
    inspectors[1] = LeafInspector(g, blade_id=96)
    inspectors[2] = LeafInspector(g, blade_id=88)
    inspectors[3] = LeafInspector(g, blade_id=80)
    inspectors[4] = LeafInspector(g, blade_id=72)
    dates = []
    # Simulation #########################################################
    for t in timer:
        # print(timer.numiter)
        # Update date
        date = (weather.next_date(t['weather'].dt, date) if date!=None else start_date)
        dates.append(date)
        print(date)
        
        # Get weather for date and add it as properties on leaves
        _, data = weather.get_weather(t['weather'].dt, date)
        set_properties(g,label = 'LeafElement',
                        wetness = data.wetness.values[0],
                        temp = data.temperature_air.values[0],
                        rain_intensity = data.rain.values[0],
                        rain_duration = data.rain_duration.values[0],
                        relative_humidity = data.relative_humidity.values[0],
                        wind_speed = data.wind_speed.values[0])

        # Grow wheat canopy
        grow_canopy(g,wheat,t['wheat'])
        update_healthy_area(g, label = 'LeafElement')
        # Note : The position of senescence goes back to its initial value after
        # a while for undetermined reason
        # --> temporary hack for keeping senescence position low when it is over
        positions = g.property('position_senescence')
        are_green = g.property('is_green')
        leaves = get_leaves(g, label = 'LeafElement')
        vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
        positions.update({vid:(0 if positions[vid]==1 and not are_green[vid] else positions[vid]) 
                          for vid in vids})
        
        # Develop disease
        if data.dispersal_event.values[0]==True and timer.numiter >= 1000 and timer.numiter <= 2000:
            # Refill pool of initial inoculum to simulate differed availability
            dispersal_units = generate_stock_du(nb_dus=10, disease=septoria)
            initiate(g, dispersal_units, inoculator)
          
        infect(g, t['disease'].dt, infection_controler, label='LeafElement')
        update(g, t['disease'].dt, growth_controler, sen_model, label='LeafElement')
        for inspector in inspectors.itervalues():
            inspector.compute_severity(g)
            inspector.update_disease_area(g)
            inspector.update_green_area(g)
            inspector.update_area(g)
        if data.dispersal_event.values[0]==True:
            disperse(g, emitter, transporter, "septoria", label='LeafElement')
            wash(g, washor, data.rain.values[0], label='LeafElement')

        # if timer.numiter%24 == 0:
            # update_plot(g)
            
    # Tout stocker dans un dataframe avec les dates en index
    outputs = {}
    for id, inspector in inspectors.iteritems():
        outs = {'severity':inspectors[id].severity,
                'disease_area': inspectors[id].leaf_disease_area,
                'green_area': inspectors[id].leaf_green_area,
                'total_area':inspectors[id].leaf_area}
        outputs[id] = pandas.DataFrame(data=outs, index=dates)
    return outputs

def run_simulation_powdery_mildew():
    # Initialization #####################################################
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Choose dates of simulation and initialize the value of date
    start_date = datetime(2000, 10, 1, 1, 00, 00)
    end_date = datetime(2001, 07, 01, 00, 00)
    date = None

    # Read weather and adapt it to septoria (add wetness)
    weather = get_septoria_weather(data_file='meteo01.csv')

    # Initialize a wheat canopy
    g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, 
                                                     width=0.2, sowing_density=150,
                                                     plant_density=150, inter_row=0.12)

    # Initialize the models for septoria
    diseases=plugin.discover('alep.disease')
    powdery_mildew = diseases['powdery_mildew'].load()
    inoculator = RandomInoculation()
    growth_controler = NoPriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')   
    emitter = PowderyMildewWindEmission()
    transporter = PowderyMildewWindDispersal(cid=0.04, k_beer=0.65, label='LeafElement')
    washor = RapillyWashing()

    # Define the schedule of calls for each model
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    weather_timing = TimeControl(delay=1, steps=nb_steps)
    wheat_timing = TimeControl(delay=24, steps=nb_steps, model=wheat, weather=weather, start_date=start_date)
    mildew_timing = TimeControl(delay=1, steps=nb_steps)
    timer = TimeControler(weather=weather_timing, wheat=wheat_timing, disease=mildew_timing)

    # Call leaf inspectors for target blades (top 3)
    inspectors = {}
    # for rank in range(1,3):
        # inspectors[rank] = LeafInspector(g, blade_id=find_blade_id(g, leaf_rank = rank, only_visible=False))
    inspectors[1] = LeafInspector(g, blade_id=96)
    inspectors[2] = LeafInspector(g, blade_id=88)
    inspectors[3] = LeafInspector(g, blade_id=80)
    inspectors[4] = LeafInspector(g, blade_id=72)
    dates = []
    # Simulation #########################################################
    for t in timer:
        # print(timer.numiter)
        # Update date
        date = (weather.next_date(t['weather'].dt, date) if date!=None else start_date)
        dates.append(date)
        print(date)
        
        # Get weather for date and add it as properties on leaves
        _, data = weather.get_weather(t['weather'].dt, date)
        set_properties(g,label = 'LeafElement',
                        wetness = data.wetness.values[0],
                        temp = data.temperature_air.values[0],
                        rain_intensity = data.rain.values[0],
                        rain_duration = data.rain_duration.values[0],
                        relative_humidity = data.relative_humidity.values[0],
                        wind_speed = data.wind_speed.values[0],
                        wind_direction = (0, 0, 1),
                        age = 1.)

        # Grow wheat canopy
        grow_canopy(g,wheat,t['wheat'])
        update_healthy_area(g, label = 'LeafElement')
        # Note : The position of senescence goes back to its initial value after
        # a while for undetermined reason
        # --> temporary hack for keeping senescence position low when it is over
        positions = g.property('position_senescence')
        are_green = g.property('is_green')
        leaves = get_leaves(g, label = 'LeafElement')
        vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
        positions.update({vid:(0 if positions[vid]==1 and not are_green[vid] else positions[vid]) 
                          for vid in vids})
        
        # Develop disease
        if data.dispersal_event.values[0]==True and timer.numiter >= 1000 and timer.numiter <= 2000:
            # Refill pool of initial inoculum to simulate differed availability
            dispersal_units = generate_stock_du(nb_dus=10, disease=powdery_mildew)
            initiate(g, dispersal_units, inoculator)
        
        du_dict = g.property('dispersal_units')
        for du_list in du_dict.values():
            for du in du_list:
                du.set_position(position=[random.random(), 0])
                
        infect(g, t['disease'].dt, infection_controler, label='LeafElement')
        update(g, t['disease'].dt, growth_controler, sen_model, label='LeafElement')
        for inspector in inspectors.itervalues():
            inspector.compute_severity(g)
            inspector.update_disease_area(g)
            inspector.update_green_area(g)
            inspector.update_area(g)
        disperse(g, emitter, transporter, "powdery_mildew", label='LeafElement')
        if data.dispersal_event.values[0]==True:
            wash(g, washor, data.rain.values[0], label='LeafElement')

        # if timer.numiter%24 == 0:
            # update_plot(g)
            
    # Tout stocker dans un dataframe avec les dates en index
    outputs = {}
    for id, inspector in inspectors.iteritems():
        outs = {'severity':inspectors[id].severity,
                'disease_area': inspectors[id].leaf_disease_area,
                'green_area': inspectors[id].leaf_green_area,
                'total_area':inspectors[id].leaf_area}
        outputs[id] = pandas.DataFrame(data=outs, index=dates)
    return outputs

def draw_outputs(out_septo, out_mildew):
    date_1 = datetime(2001, 2, 1, 1, 00, 00)
    date_2 = datetime(2001, 7, 1, 0, 00, 00)
    date_seq = pandas.date_range(date_1,date_2, freq='H')
    months = MonthLocator(bymonthday=(1,10,20))
    month_fmt = DateFormatter('%b-%d')

    # fig = plt.figure()
    # ax1 = fig.add_subplot(1,2,1)
    # ax1.plot(date_seq, out_septo[1].severity[date_1:date_2], color='b', linewidth=2, label='leaf 1')
    # ax1.plot(date_seq, out_septo[2].severity[date_1:date_2], color='g', linewidth=2, label='leaf 2')
    # ax1.plot(date_seq, out_septo[3].severity[date_1:date_2], color='r', linewidth=2, label='leaf 3')
    # ax1.plot(date_seq, out_septo[4].severity[date_1:date_2], color='k', linewidth=2, label='leaf 4')
    # ax1.set_xticklabels(ax1.get_xticklabels(), fontsize=14, rotation=30, ha='right')
    # ax1.xaxis.set_major_locator(months)
    # ax1.xaxis.set_major_formatter(month_fmt)
    # ax1.set_ylabel('Severity for septoria (in cm2)', fontsize=20)
    # lim = ax1.set_ylim((0,100))
    # h, l = ax1.get_legend_handles_labels()
    # ax1.legend(h,l)

    # ax2 = fig.add_subplot(1,2,2)
    # ax2.plot(date_seq, out_mildew[1].severity[date_1:date_2], color='b', linewidth=2, label='leaf 1')
    # ax2.plot(date_seq, out_mildew[2].severity[date_1:date_2], color='g', linewidth=2, label='leaf 2')
    # ax2.plot(date_seq, out_mildew[3].severity[date_1:date_2], color='r', linewidth=2, label='leaf 3')
    # ax2.plot(date_seq, out_mildew[4].severity[date_1:date_2], color='k', linewidth=2, label='leaf 4')
    # ax2.set_xticklabels(ax2.get_xticklabels(), fontsize=14, rotation=30, ha='right')
    # ax2.xaxis.set_major_locator(months)
    # ax2.xaxis.set_major_formatter(month_fmt)
    # ax2.set_ylabel('Severity for powdery mildew(in cm2)', fontsize=20)
    # lim = ax2.set_ylim((0,100))
    # h, l = ax2.get_legend_handles_labels()
    # ax2.legend(h,l)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    ax1.plot(date_seq, out_septo[1].disease_area[date_1:date_2], color='b', linewidth=2, label='leaf 1')
    ax1.plot(date_seq, out_septo[2].disease_area[date_1:date_2], color='g', linewidth=2, label='leaf 2')
    ax1.plot(date_seq, out_septo[3].disease_area[date_1:date_2], color='r', linewidth=2, label='leaf 3')
    ax1.plot(date_seq, out_septo[4].disease_area[date_1:date_2], color='k', linewidth=2, label='leaf 4')
    ax1.set_xticklabels(ax1.get_xticklabels(), fontsize=14, rotation=30, ha='right')
    ax1.tick_params(axis='y', labelsize=14)
    ax1.xaxis.set_major_locator(months)
    ax1.xaxis.set_major_formatter(month_fmt)
    ax1.set_ylabel('Disease area for septoria (in cm2)', fontsize=20)
    lim = ax1.set_ylim((0,25))
    h, l = ax1.get_legend_handles_labels()
    ax1.legend(h,l, loc=2)
    
    # ax3 = twinx(ax1)
    # ax3.plot(date_seq, out_septo[1].green_area[date_1:date_2], 'b--', linewidth=2, label='leaf 1')
    # ax3.plot(date_seq, out_septo[2].green_area[date_1:date_2], 'g--', linewidth=2, label='leaf 2')
    # ax3.plot(date_seq, out_septo[3].green_area[date_1:date_2], 'r--', linewidth=2, label='leaf 3')
    # ax3.plot(date_seq, out_septo[4].green_area[date_1:date_2], 'k--', linewidth=2, label='leaf 4')
    
    ax2 = fig.add_subplot(1,2,2)
    ax2.plot(date_seq, out_mildew[1].disease_area[date_1:date_2], color='b', linewidth=2, label='leaf 1')
    ax2.plot(date_seq, out_mildew[2].disease_area[date_1:date_2], color='g', linewidth=2, label='leaf 2')
    ax2.plot(date_seq, out_mildew[3].disease_area[date_1:date_2], color='r', linewidth=2, label='leaf 3')
    ax2.plot(date_seq, out_mildew[4].disease_area[date_1:date_2], color='k', linewidth=2, label='leaf 4')
    ax2.set_xticklabels(ax2.get_xticklabels(), fontsize=14, rotation=30, ha='right')
    ax2.tick_params(axis='y', labelsize=14)
    ax2.xaxis.set_major_locator(months)
    ax2.xaxis.set_major_formatter(month_fmt)
    ax2.set_ylabel('Disease area for powdery mildew(in cm2)', fontsize=20)
    lim = ax2.set_ylim((0,25))
    h, l = ax2.get_legend_handles_labels()
    ax2.legend(h,l)

    ax3 = twinx(ax1)
    ax3.plot(date_seq, out_septo[1].total_area[date_1:date_2], 'b--', linewidth=2, label='leaf 1')
    ax3.plot(date_seq, out_septo[2].total_area[date_1:date_2], 'g--', linewidth=2, label='leaf 2')
    ax3.plot(date_seq, out_septo[3].total_area[date_1:date_2], 'r--', linewidth=2, label='leaf 3')
    ax3.plot(date_seq, out_septo[4].total_area[date_1:date_2], 'k--', linewidth=2, label='leaf 4')
    
    show(False)