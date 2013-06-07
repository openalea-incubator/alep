"""
Make images of a simulation of wheat/septoria epidemics for a movie.
"""

from alinea.astk.TimeControl import *

from alinea.alep.wheat_examples import adel_mtg, adel_mtg2, adel_one_leaf
from alinea.alep.inoculation import InoculationFirstLeaves
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.dispersal import RandomDispersal
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.septoria import *
from alinea.alep.protocol import *
from alinea.adel.mtg_interpreter import plot3d

from alinea.alep.disease_operation import *
from alinea.alep.disease_outputs import *
from alinea.alep.architecture import *
from alinea.alep.alep_color import alep_colormap, green_white

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
###################################################################################
# Define a plant or canopy
g = adel_mtg2()

# Add missing properties needed for the simulation
set_properties(g,label = 'LeafElement',
               surface=5., healthy_surface=5., position_senescence=None)
set_properties(g,label = 'LeafElement',
                wetness=True,
                temp=22.,
                relative_humidity=85.,
                wind_speed=0.5)

                   
# Create a pool of dispersal units (DU)
diseases=plugin.discover('alep.disease')
septoria = diseases['septoria_exchanging_rings'].load()
nb_du = 10
dispersal_units = generate_stock_du(nb_du, disease=septoria)

# Distribute the DU 
inoculator = InoculationFirstLeaves()
initiate(g, dispersal_units, inoculator)

# Simulation
controler = NoPriorityGrowthControl()
dispersor = Septo3DSplash(reference_surface=1./200)
# dispersor = RandomDispersal()

nsteps = 2000

wheat_timing = TimeControl(delay=1, steps = nsteps)
septo_timing = TimeControl(delay=1, steps = nsteps)
plot_timing = TimeControl(delay=24, steps = nsteps)
timer = TimeControler(wheat = wheat_timing, disease = septo_timing, ploting = plot_timing)

for t in timer:
    print(timer.numiter)
    if timer.numiter>400 and timer.numiter%100==0:
        rain_intensity = 2.
        rain_duration = 1.
    else:
        rain_intensity = 0.
        rain_duration = 0.
    set_properties(g,label = 'LeafElement',
                wetness=True,
                temp=22.,
                relative_humidity=85.,
                rain_intensity = rain_intensity,
                rain_duration = rain_duration,
                wind_speed=0.5)

    infect(g, t['disease'].dt, label='LeafElement')
    update(g, t['disease'].dt, controler, label='LeafElement')

    if rain_intensity>0.:
        disperse(g, dispersor, "septoria", label='LeafElement')  
    
    if t['ploting'].dt > 0:
        print('ploting...')
        plot_lesions(g)
        index = timer.numiter/24
        if index>9:
            index+=900
        scene = plot3d(g)
        save_image(scene, image_name='image%d.png' % index)