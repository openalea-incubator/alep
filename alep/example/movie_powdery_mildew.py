"""
Demonstrate a simulation of vine/oidium epidemics
"""

from alinea.alep.vine import Vine
from alinea.astk.TimeControl import *

from alinea.alep.inoculation import RandomInoculation
from alinea.alep.growth_control import GrowthControlVineLeaf
from alinea.alep.dispersal import RandomDispersal
from alinea.alep.powdery_mildew import *
from alinea.alep.protocol import *
from alinea.adel.mtg_interpreter import plot3d

from alinea.alep.disease_operation import *
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


def update_plot(g):
    # Count lesions by id & add it as MTG property ####################################
    # nb_lesions_by_leaf = count_lesions_by_leaf(g, label = 'lf')
    surface_lesions_by_leaf = count_lesion_surfaces_by_leaf(g, label = 'lf')
    # print(surface_lesions_by_leaf)
    # print(g.property('surface'))
    set_property_on_each_id(g, 'surface_lesions', surface_lesions_by_leaf, label = 'lf')
                       
    # Visualization ###################################################################
    g = alep_colormap(g, 'surface_lesions', cmap=green_white(levels=10), lognorm=False)
    brown = (100,70,30)
    trunk_ids = [n for n in g if g.label(n).startswith('tronc')]
    for id in trunk_ids:
        trunk = g.node(id)
        trunk.color = brown
    scene = plot3d(g)
    Viewer.display(scene)
    return scene

vine = Vine()
g0 = vine.setup_canopy(age=6)
# vine.plot(g0)

# Add missing properties needed for the simulation
# The simulation requires the following properties on leaf elements:
#   - 'surface': total surface of the leaf element (in cm2)
#   - 'healthy_surface': surface of the leaf element without lesion or senescence (in cm2)
#   - 'age': age of the leaf (in decimal_days)
#   - 'position_senescence': position of the senescence on blade axis
set_properties(g0,label = 'lf',
               surface=5., healthy_surface=5., position_senescence=None)

fungus = powdery_mildew()
PowderyMildewDU.fungus=fungus
nb_dus = 100
dispersal_units = ([PowderyMildewDU(nb_spores=1, status="emitted") for i in range(nb_dus)])
inoculator=RandomInoculation()
initiate(g0,dispersal_units,inoculator, label='lf')


controler = GrowthControlVineLeaf()
dispersor = RandomDispersal()
# nsteps = 49
nsteps = 1000

vine_timing = TimeControl(delay=24, steps = nsteps)
mildew_timing = TimeControl(delay =1, steps = nsteps)
plot_timing = TimeControl(delay=24, steps = nsteps)
timer = TimeControler(vine = vine_timing, disease = mildew_timing, ploting = plot_timing)


g=g0
# scene = vine.generate_scene(g)
scene = plot3d(g)

for t in timer:
    print(timer.numiter)
    set_properties_on_new_leaves(g,label = 'lf',
                             surface=5., healthy_surface=5.,
                             position_senescence=None)
    set_properties(g,label = 'lf',
                    wetness=True,
                    temp=22.,
                    rain_intensity=0.,
                    rain_duration=0.,
                    relative_humidity=85.,
                    wind_speed=0.5)
                        
    infect(g, t['disease'].dt, label='lf')
    update(g, t['disease'].dt, controler, label='lf')
    disperse(g, dispersor, "powdery_mildew", label='lf')  
    
    g = vine.grow(g,t['vine'])
    # scene = vine.generate_scene(g)
    if t['ploting'].dt > 0:
        print('ploting...')
        scene = update_plot(g)
        index = timer.numiter/24
        save_image(scene, image_name='image%d.png' % index)