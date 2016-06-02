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
from alinea.astk.TimeControl import *
import alinea.septo3d
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather
from alinea.alep.alep_weather import wetness_rapilly, is_raining

# Imports for wheat
from alinea.alep.wheat import initialize_stand, find_blade_id, find_leaf_ids
from alinea.astk.plant_interface import grow_canopy
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.disease_outputs import LeafInspector, save_image
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence

# Useful functions ############################################################
from alinea.alep.architecture import set_property_on_each_id
from alinea.alep.disease_outputs import compute_severity_by_leaf
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer
def update_plot(g):
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

def run_simulation():
    # Initialization #####################################################
    # Set the seed of the simulation
    rd.seed(0)
    np.random.seed(0)

    # Read weather and adapt it to septoria (add wetness)
    meteo_path = shared_data(alinea.septo3d, 'meteo98-99.txt')
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    seq = pandas.date_range(start = "1998-10-01 01:00:00",
                            end = "1999-07-01 01:00:00", 
                            freq='H')

    # Initialize a wheat canopy
    g, wheat, domain_area, domain = initialize_stand(age=0., length=0.1, 
                                                    width=0.2, sowing_density=150,
                                                    plant_density=150, inter_row=0.12, 
                                                    seed=8)

    # Initialize the models for septoria
    # septoria = new_septoria(senescence_treshold=senescence_treshold)
    septoria = plugin_septoria()
    inoculator = RandomInoculation()
    growth_controler = NoPriorityGrowthControl()
    infection_controler = BiotrophDUPositionModel()
    sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
    emitter = SeptoriaRainEmission(domain_area=domain_area)
    transporter = Septo3DSplash(reference_surface=domain_area)
    washor = RapillyWashing()

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
    # for rank in range(1,3):
        # inspectors[rank] = LeafInspector(g, blade_id=find_blade_id(g, leaf_rank = rank, only_visible=False))
    inspectors[1] = LeafInspector(g, blade_id=96)
    # inspectors[2] = LeafInspector(g, blade_id=88)
    # inspectors[3] = LeafInspector(g, blade_id=80)
    dates = []
    # Simulation #########################################################
    for i,controls in enumerate(zip(weather_timing, wheat_timing, 
                                    septo_timing, rain_timing)):
        weather_eval, wheat_eval, septo_eval, rain_eval = controls
        
        # Update date
        date = weather_eval.value.index[0]
        dates.append(date)
        
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
            print(date)
            g,_ = grow_canopy(g, wheat, wheat_eval.value)
            # Note : The position of senescence goes back to its initial value after
            # a while for undetermined reason
            # --> temporary hack for keeping senescence position low when it is over
            positions = g.property('position_senescence')
            are_green = g.property('is_green')
            areas = g.property('area')
            senesced_areas = g.property('senesced_area')
            leaves = get_leaves(g, label = 'LeafElement')
            vids = [leaf for leaf in leaves if leaf in g.property('geometry')]
            positions.update({vid:(0 if (positions[vid]==1 and not are_green[vid]) or
                                        (positions[vid]>0 and round(areas[vid],5)==round(senesced_areas[vid],5))
                                        else positions[vid]) for vid in vids})
            
        # Develop disease
        if septo_eval:
            sen_model.find_senescent_lesions(g, label = 'LeafElement')
            update_healthy_area(g, label = 'LeafElement')
            if rain_eval and i <= 500:
                # Refill pool of initial inoculum to simulate differed availability
                dispersal_units = generate_stock_du(nb_dus=rd.randint(0,5), disease=septoria)
                initiate(g, dispersal_units, inoculator)
            infect(g, septo_eval.dt, infection_controler, label='LeafElement')
            update(g, septo_eval.dt, growth_controler, sen_model, label='LeafElement')
        if rain_eval:
            if rain_eval.value.rain.mean()>0:
                g, nb = disperse(g, emitter, transporter, "septoria", label='LeafElement')
                for inspector in inspectors.itervalues():
                    inspector.update_du_variables(g)
                wash(g, washor, rain_eval.value.rain.mean(), label='LeafElement')
                # Save outputs after washing
                infection_controler.control_position(g)
                for inspector in inspectors.itervalues():
                    inspector.update_du_variables(g)
                    inspector.update_green_area(g)
                    inspector.update_healthy_area(g)
            else:
                for inspector in inspectors.itervalues():
                    inspector.nb_dus += [0, 0]
                    inspector.nb_dus_on_green += [0, 0]
                    inspector.nb_dus_on_healthy += [0, 0]
                    inspector.update_green_area(g)
                    inspector.update_healthy_area(g)
        else:
            for inspector in inspectors.itervalues():
                inspector.nb_dus += [0, 0]
                inspector.nb_dus_on_green += [0, 0]
                inspector.nb_dus_on_healthy += [0, 0]
                inspector.update_green_area(g)
                inspector.update_healthy_area(g)
        
        for inspector in inspectors.itervalues():
            inspector.compute_nb_infections(g)
        
        if wheat_eval:
            update_plot(g)
            # scene = plot3d(g)
            # index = i/24
            # if index < 10 :
                # image_name='./images_septo2/image0000%d.png' % index
            # elif index < 100 :
                # image_name='./images_septo2/image000%d.png' % index
            # elif index < 1000 :
                # image_name='./images_septo2/image00%d.png' % index
            # elif index < 10000 :
                # image_name='./images_septo2/image0%d.png' % index
            # else :
                # image_name='./images_septo/image%d.png' % index
            # save_image(scene, image_name=image_name)
            
    # Tout stocker dans un dataframe avec les dates en index
    outs = {'nb_dus':inspectors[1].nb_dus[::2],
            'nb_unwashed_dus':inspectors[1].nb_dus[1::2],
            'nb_dus_on_healthy':inspectors[1].nb_dus_on_healthy[1::2],
            'nb_infections':inspectors[1].nb_infections,
            'green_area':inspectors[1].leaf_green_area,
            'healthy_area':inspectors[1].leaf_healthy_area}
    outputs = pandas.DataFrame(data=outs, index=dates)
    return outputs

def draw_outputs(outputs):    
    date_1 = datetime.datetime(1999, 3, 1, 1, 00, 00)
    date_2 = datetime.datetime(1999, 7, 1, 0, 00, 00)
    date_seq = pandas.date_range(date_1,date_2, freq='H')
    months = MonthLocator(bymonthday=(1,10,20))
    month_fmt = DateFormatter('%b-%d')

    fig = plt.figure()
    # ymax = 300
    ymax = max(outputs.nb_dus[date_1:date_2])+100
    ax1 = fig.add_subplot(2,2,1)
    ax1.set_yscale('symlog')
    ax1.vlines(date_seq, [0], outputs.nb_dus[date_1:date_2], color='b', linewidth=2)
    ax1.set_xticklabels(ax1.get_xticklabels(), fontsize=15, rotation=30, ha='right')
    ax1.xaxis.set_major_locator(months)
    ax1.xaxis.set_major_formatter(month_fmt)
    ax1.set_ylabel('Total intercepted', fontsize=20)
    lim = ax1.set_ylim((0,ymax))
    

    ax2 = fig.add_subplot(2,2,2)
    ax2.set_yscale('symlog')
    ax2.vlines(date_seq, [0], outputs.nb_unwashed_dus[date_1:date_2], color='b', linewidth=2)
    lim = ax2.set_ylim((0,ymax))
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=30, ha='right')
    ax2.xaxis.set_major_locator(months)
    ax2.xaxis.set_major_formatter(month_fmt)
    ax2.set_ylabel('Total unwashed', fontsize=20)

    ax3 = fig.add_subplot(2,2,3)
    ax3.set_yscale('symlog')
    ax3.vlines(date_seq, [0], outputs.nb_dus_on_healthy[date_1:date_2], color='b', linewidth=2)
    ax3.set_xticklabels(ax3.get_xticklabels(), fontsize=15, rotation=30, ha='right')
    ax3.xaxis.set_major_locator(months)
    ax3.xaxis.set_major_formatter(month_fmt)
    lim = ax3.set_ylim((0,ymax))
    ax3.set_ylabel('Total on healthy tissue', fontsize=20)

    ax4 = twinx(ax3)
    ax4.plot(date_seq, outputs.green_area[date_1:date_2], color='g',
                linestyle='--', linewidth=2, label="green area")
    ax4.plot(date_seq, outputs.healthy_area[date_1:date_2], color='b', 
                linestyle='--', linewidth=2, label="healthy area")
    ax4.set_xticklabels(ax4.get_xticklabels(), fontsize=15, rotation=30, ha='right')
    ax4.xaxis.set_major_locator(months)
    ax4.xaxis.set_major_formatter(month_fmt)
    h, l = ax4.get_legend_handles_labels()
    ax4.legend(h,l)
    ax4.set_ylabel('Leaf area (in cm2)', fontsize=20)

    ax5 = fig.add_subplot(2,2,4)
    ax5.set_yscale('symlog')
    ax5.vlines(date_seq, [0], outputs.nb_infections[date_1:date_2], color='b', linewidth=2)
    ax5.set_xticklabels(ax5.get_xticklabels(), fontsize=15, rotation=30, ha='right')
    ax5.xaxis.set_major_locator(months)
    ax5.xaxis.set_major_formatter(month_fmt)
    # lim = ax5.set_ylim((0,20))
    lim = ax3.set_ylim((0,ymax))
    ax5.set_ylabel('Total incubating', fontsize=20)

    show(False)
