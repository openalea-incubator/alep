""" Make images of a simulation of vine/oidium epidemics for a movie.
"""
# Useful imports
import random as rd
import numpy as np
import pandas
from alinea.astk.TimeControl import *

# Imports for weather
from datetime import datetime
from alinea.weather.global_weather import *
from alinea.alep.alep_weather import add_wetness

# Imports for grapevine
from alinea.alep.vine import Vine
from alinea.astk.plant_interface import new_canopy, grow_canopy
from alinea.alep.architecture import (update_healthy_area, add_area_topvine,
                                      set_properties, set_property_on_each_id)

# Imports for powdery mildew
from openalea.vpltk import plugin
from alinea.alep.protocol import *
from alinea.alep.powdery_mildew import *
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal import PowderyMildewWindDispersal
from alinea.alep.growth_control import GrowthControlVineLeaf
from alinea.alep.infection_control import BiotrophDUProbaModel

# Imports for movie
from alinea.alep.disease_outputs import compute_severity_by_leaf
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer

# /!\ TEMP /!\ ################################################################
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
    severity_by_leaf = compute_severity_by_leaf(g, label = 'lf')
    print(severity_by_leaf)
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

# Initialization ##############################################################
# Set the seed of the simulation
rd.seed(0)
np.random.seed(0)

# Choose dates of simulation and initialize the value of date
start_date = datetime(2001, 05, 1, 1, 00, 00)
# end_date = datetime(2001, 06, 1, 00, 00)
end_date = datetime(2001, 07, 01, 00, 00)
date = start_date

# Read weather and adapt it to septoria (add wetness)
weather = Weather(data_file='meteo01.csv')
weather = add_wetness(weather)

# Initialize a vine canopy
vine = Vine()
g,_ = new_canopy(vine, age = 6)

# Initialize the models for powdery mildew
diseases=plugin.discover('alep.disease')
powdery_mildew = diseases['powdery_mildew'].load()
inoculator = RandomInoculation()
growth_controler = GrowthControlVineLeaf()
infection_controler = BiotrophDUProbaModel()
dispersor = PowderyMildewWindDispersal()

# Define the schedule of calls for each model
nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
weather_timing = TimeControl(delay=1, steps=nb_steps)
vine_timing = TimeControl(delay=24, steps=nb_steps)
mildew_timing = TimeControl(delay=1, steps=nb_steps)
plot_timing = TimeControl(delay=24, steps=nb_steps)
timer = TimeControler(weather=weather_timing, vine=vine_timing, disease=mildew_timing, plotting=plot_timing)

# Simulation #########################################################
for t in timer:
    # Update date
    date = (weather.next_date(t['weather'].dt, date) if date!=None else start_date)
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
    # if timer.numiter < 1500 and timer.numiter%50==0:
        # # Refill pool of initial inoculum to simulate differed availability
        # dispersal_units = generate_stock_du(nb_dus=10, disease=powdery_mildew)
        # initiate(g, dispersal_units, inoculator, label='lf')
    
    if timer.numiter < 500 and labels[188]=='lf' and timer.numiter%50==0:
        # Refill pool of initial inoculum to simulate differed availability
        dispersal_units = generate_stock_du(nb_dus=10000, disease=powdery_mildew)
        initiate(g, dispersal_units, inoculator, label='lf')
        # try:
            # g.node(188).dispersal_units += dispersal_units
        # except: 
            # g.node(188).dispersal_units = dispersal_units
    
    from alinea.alep.disease_outputs import count_dispersal_units, count_lesions, count_lesions_by_leaf
    # print('nb_lesions : %d' % count_lesions(g))
    nb = count_lesions_by_leaf(g)
    if len(nb.values())==0:
        print(0)
    else:
        print(max(nb.values()), max(nb.iterkeys(), key=lambda k: nb[k]))  
    infect(g, t['disease'].dt, infection_controler, label='lf')
    update(g, t['disease'].dt, growth_controler, label='lf')
    disperse(g, dispersor, "powdery_mildew", label='lf')

    from alinea.alep.disease_outputs import plot_lesions
    if t['plotting'].dt > 0:
        update_plot(g)
        scene = plot3d(g)
        index = timer.numiter/24
        if index < 10 :
            image_name='./images_oidium/image00%d.png' % index
        elif index < 100 :
            image_name='./images_oidium/image0%d.png' % index
        else :
            image_name='./images_oidium/image%d.png' % index
        save_image(scene, image_name=image_name)