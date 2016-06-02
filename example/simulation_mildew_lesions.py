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
from alinea.weather.global_weather import *
from alinea.alep.alep_weather import add_wetness

# Imports for vine
from alinea.alep.vine import Vine
from alinea.astk.plant_interface import new_canopy, grow_canopy
from alinea.alep.architecture import (update_healthy_area, add_area_topvine,
                                      set_properties, set_property_on_each_id)

# Imports for powdery mildew
from openalea.vpltk import plugin
from alinea.alep.protocol import *
from alinea.alep.powdery_mildew import *
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.disease_outputs import VineLeafInspector
from alinea.alep.inoculation import InoculationYoungLeaves
from alinea.alep.dispersal_emission import PowderyMildewWindEmission
from alinea.alep.dispersal_transport import PowderyMildewWindDispersal
from alinea.alep.growth_control import GrowthControlVineLeaf
from alinea.alep.infection_control import BiotrophDUProbaModel

def update_plot(g):
    from alinea.alep.architecture import set_property_on_each_id
    from alinea.alep.disease_outputs import compute_severity_by_leaf
    from alinea.adel.mtg_interpreter import plot3d
    from openalea.plantgl.all import Viewer

    # Compute severity by leaf
    severity_by_leaf = compute_severity_by_leaf(g, label = 'lf')
    set_property_on_each_id(g, 'severity', severity_by_leaf, label = 'lf')
                       
    # Visualization
    g = alep_colormap(g, 'severity', cmap=green_yellow_red(levels=100),
                      lognorm=False, zero_to_one=False, vmax=100)
    brown = (100,70,30)
    trunk_ids = [n for n in g if g.label(n).startswith('tronc')]
    for id in trunk_ids:
        trunk = g.node(id)
        trunk.color = brown
    
    scene = plot3d(g)
    Viewer.display(scene)
    return scene

class NoDeposit:
    def disperse(self, g, dispersal_units, time_control = None):
        deposits = {}
        return deposits
    
def run_simulation():
    # Initialization #####################################################
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Choose dates of simulation and initialize the value of date
    start_date = datetime(2001, 03, 1, 1, 00, 00)
    end_date = datetime(2001, 05, 1, 00, 00, 00)
    # end_date = datetime(2001, 07, 01, 00, 00)
    date = None

    # Read weather and adapt it to powdery mildew (add wetness)
    weather = Weather(data_file='meteo02.csv')
    weather = add_wetness(weather)
    
    # Initialize a vine canopy
    vine = Vine()
    g,_ = new_canopy(vine, age = 1)

    # Initialize the models for septoria
    diseases=plugin.discover('alep.disease')
    powdery_mildew = diseases['powdery_mildew'].load()
    inoculator = InoculationYoungLeaves(age_max=5.)
    growth_controler = GrowthControlVineLeaf()
    infection_controler = BiotrophDUProbaModel()
    emitter = PowderyMildewWindEmission()
    # transporter = PowderyMildewWindDispersal()
    transporter = NoDeposit()
    
    # Define the schedule of calls for each model
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    weather_timing = TimeControl(delay=1, steps=nb_steps)
    vine_timing = TimeControl(delay=24, steps=nb_steps)
    mildew_timing = TimeControl(delay=1, steps=nb_steps)
    timer = TimeControler(weather=weather_timing, vine=vine_timing, disease = mildew_timing)

    # Call leaf inspectors for target blades (top 3)
    target_leaf = 6
    inspector = VineLeafInspector(leaf_id=target_leaf)
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
        set_properties(g,label = 'lf',
                        wetness = data.wetness.values[0],
                        temp = data.temperature_air.values[0],
                        rain_intensity = data.rain.values[0],
                        relative_humidity = data.relative_humidity.values[0],
                        wind_speed = data.wind_speed.values[0],
                        wind_direction = (1,0,0))

        # Grow vine canopy
        g,_ = grow_canopy(g, vine, t['vine'])
        add_area_topvine(g)
        update_healthy_area(g, label = 'lf')
        
        # Develop disease
        labels = g.property('label')
        # if timer.numiter < 1500 and labels[target_leaf]=='lf' and timer.numiter%50==0:
        if timer.numiter > 100 and timer.numiter < 170 and labels[target_leaf]=='lf' and timer.numiter%rd.randint(2,5)==0:
        # if timer.numiter < 1500 and timer.numiter%50==0:
            # Refill pool of initial inoculum to simulate differed availability
            dispersal_units = generate_stock_du(nb_dus=rd.randint(1,5), disease=powdery_mildew)
            try:
                g.node(target_leaf).dispersal_units += dispersal_units
            except:
                g.node(target_leaf).dispersal_units = dispersal_units
            # dispersal_units = generate_stock_du(nb_dus=1000, disease=powdery_mildew)
            # initiate(g, dispersal_units, inoculator, label='lf')
            
        infect(g, t['disease'].dt, infection_controler, label='lf')
        update(g, t['disease'].dt, growth_controler, label='lf')
        disperse(g, emitter, transporter, "powdery_mildew", label='lf')
        
        inspector.update_data(g)
            
        if timer.numiter%24 == 0:
            update_plot(g)
        
        outs = {'ratio_latent':inspector.ratio_latent,
                'ratio_spo':inspector.ratio_spo,
                'ratio_empty':inspector.ratio_empty}
        outputs = pandas.DataFrame(data=outs, index=dates)                   
    return outputs

def draw_outputs(outputs):
    date_1 = datetime(2001, 03, 5, 1, 00, 00)
    date_2 = datetime(2001, 04, 10, 00, 00, 00)
    date_seq = pandas.date_range(date_1,date_2, freq='H')
    months = MonthLocator(bymonthday=([1]+range(5,30,5)))
    month_fmt = DateFormatter('%b-%d')

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(date_seq, outputs.ratio_latent[date_1:date_2], color='g', linewidth=2, label='latent')
    ax1.plot(date_seq, outputs.ratio_spo[date_1:date_2], color='r', linewidth=2, label='sporulating')
    ax1.plot(date_seq, outputs.ratio_empty[date_1:date_2], color='k', linewidth=2, label='empty')
    ax1.set_xticklabels(ax1.get_xticklabels(), fontsize=20, rotation=30, ha='right')
    ax1.xaxis.set_major_locator(months)
    ax1.xaxis.set_major_formatter(month_fmt)
    ax1.set_ylabel('% Disease coverage', fontsize=30)   
    lim = ax1.set_ylim((0,max(outputs.ratio_empty)+10))
    ax1.tick_params(axis='y', labelsize=20)
    h, l = ax1.get_legend_handles_labels()
    ax1.legend(h,l, prop={'size':20})

    # fig = plt.figure()
    # ax1 = fig.add_subplot(1,3,1)
    # ax1.plot(date_seq, outputs.latent[date_1:date_2], color='g', linewidth=2)
    # ax1.set_xticklabels(ax1.get_xticklabels(), fontsize=14, rotation=30, ha='right')
    # ax1.xaxis.set_major_locator(months)
    # ax1.xaxis.set_major_formatter(month_fmt)
    # ax1.set_ylabel('% Latent', fontsize=20)
    # lim = ax1.set_ylim((0,100))
        
    # ax2 = fig.add_subplot(1,3,2)
    # ax2.plot(date_seq, outputs.ratio_spo[date_1:date_2], color='r', linewidth=2)
    # ax2.set_xticklabels(ax2.get_xticklabels(), fontsize=14, rotation=30, ha='right')
    # ax2.xaxis.set_major_locator(months)
    # ax2.xaxis.set_major_formatter(month_fmt)
    # ax2.set_ylabel('% Sporulating', fontsize=20)
    # lim = ax2.set_ylim((0,100))

    # ax3 = fig.add_subplot(1,3,3)
    # ax3.plot(date_seq, outputs.ratio_empty[date_1:date_2], color='k', linewidth=2)
    # ax3.set_xticklabels(ax3.get_xticklabels(), fontsize=14, rotation=30, ha='right')
    # ax3.xaxis.set_major_locator(months)
    # ax3.xaxis.set_major_formatter(month_fmt)
    # ax3.set_ylabel('% Empty', fontsize=20)
    # lim = ax3.set_ylim((0,100))
    show(False)
