"""
Make images of a simulation of wheat/septoria epidemics for a movie.
"""

# Useful imports
import random as rd
import numpy as np
import pandas
from alinea.astk.TimeControl import *

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
from alinea.alep.inoculation import RandomInoculation
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence

# Imports for movie
from alinea.alep.architecture import set_property_on_each_id
from alinea.alep.disease_outputs import compute_severity_by_leaf
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer

# /!\ TEMP /!\ #################################################################
from openalea.plantgl.all import Viewer
import os.path

def save_image(scene, image_name='%s/img%04d.%s', directory='.', index=0, ext='png'):
    '''
    Save an image of a scene in a specific directory

    Parameters
    ----------

        - scene: a PlantGL scene
        - image_name: a string template 
            The format of the string is dir/img5.png
        - directory (optional: ".") the directory where the images are written
        - index: the index of the image
        - ext : the image format

    Example
    -------

        - Movie:
            convert *.png movie.mpeg
            convert *.png movie.gif
            mencoder "mf://*.png" -mf type=png:fps=25 -ovc lavc -o output.avi
            mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq "mf://*.png" -mf type=png:fps=18 -of avi  -o output.avi
            
    '''

    if not image_name:
        image_name='{directory}/img{index:0>4d}.{ext}'
    filename = image_name.format(directory=directory, index=index, ext=ext)
    Viewer.frameGL.saveImage(filename)
    return scene,

# Useful functions ############################################################
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

# Initialization ##############################################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Choose dates of simulation and initialize the value of date
start_date = datetime(2000, 10, 1, 1, 00, 00)
end_date = datetime(2000, 10, 31, 00, 00)
# end_date = datetime(2001, 07, 01, 00, 00)
date = start_date

# Read weather and adapt it to septoria (add wetness)
weather = get_septoria_weather(data_file='meteo01.csv')

# Initialize a wheat canopy
g, wheat, domain_area = initialize_stand(age=0., length=0.1,
                                        width=0.2, sowing_density=150,
                                        plant_density=150, inter_row=0.12)

# Initialize the models for septoria
septoria = plugin_septoria()
inoculator = RandomInoculation()
growth_controler = NoPriorityGrowthControl()
infection_controler = BiotrophDUPositionModel()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
dispersor = Septo3DSplash(reference_surface=domain_area)

# Define the schedule of calls for each model
nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
weather_timing = TimeControl(delay=1, steps=nb_steps)
wheat_timing = TimeControl(delay=24, steps=nb_steps, model=wheat, weather=weather, start_date=start_date)
septo_timing = TimeControl(delay=1, steps=nb_steps)
plot_timing = TimeControl(delay=24, steps=nb_steps)
timer = TimeControler(weather=weather_timing, wheat=wheat_timing, disease = septo_timing, plotting=plot_timing)

# Initialization #####################################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Choose dates of simulation and initialize the value of date
start_date = datetime(2000, 10, 1, 1, 00, 00)
# end_date = datetime(2000, 10, 31, 00, 00)
end_date = datetime(2001, 07, 01, 00, 00)
date = None

# Read weather and adapt it to septoria (add wetness)
weather = get_septoria_weather(data_file='meteo01.csv')

# Initialize a wheat canopy
g, wheat, domain_area = initialize_stand(age=0., length=0.1, 
                                        width=0.2, sowing_density=150,
                                        plant_density=150, inter_row=0.12)

# Initialize the models for septoria
septoria = plugin_septoria()
inoculator = RandomInoculation()
growth_controler = NoPriorityGrowthControl()
infection_controler = BiotrophDUPositionModel()
sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
dispersor = Septo3DSplash(reference_surface=domain_area)
washor = RapillyWashing()

# Define the schedule of calls for each model
nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
weather_timing = TimeControl(delay=1, steps=nb_steps)
wheat_timing = TimeControl(delay=24, steps=nb_steps, model=wheat, weather=weather, start_date=start_date)
septo_timing = TimeControl(delay=1, steps=nb_steps)
timer = TimeControler(weather=weather_timing, wheat=wheat_timing, disease = septo_timing)

# Simulation #########################################################
for t in timer:
    # Update date
    date = (weather.next_date(t['weather'].dt, date) if date!=None else start_date)
    
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
    if data.dispersal_event.values[0]==True and timer.numiter <= 1500:
        # Refill pool of initial inoculum to simulate differed availability
        dispersal_units = generate_stock_du(nb_dus=20, disease=septoria)
        initiate(g, dispersal_units, inoculator)
      
    infect(g, t['disease'].dt, infection_controler, label='LeafElement')
    update(g, t['disease'].dt, growth_controler, sen_model, label='LeafElement')
    if data.dispersal_event.values[0]==True:
        disperse(g, dispersor, "septoria", label='LeafElement')
        wash(g, washor, data.rain.values[0], label='LeafElement')

    if t['wheat'].dt > 0:
        update_plot(g)
        scene = plot3d(g)
        index = timer.numiter/24
        if index < 10 :
            image_name='./images_septo/image00%d.png' % index
        elif index < 100 :
            image_name='./images_septo/image0%d.png' % index
        else :
            image_name='./images_septo/image%d.png' % index
        save_image(scene, image_name=image_name)