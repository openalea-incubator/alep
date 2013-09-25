""" Demonstration of the dispersal of spores of powdery mildew by wind.

    1. Generate a vine canopy
    2. Choose an initial source of inoculum
        --> Set a lesion with 1000 dispersal units
    3. Emit spores in the direction of the wind
        --> Call the strategy of dispersal 'PowderyMildewWindDispersal'
    4. Each DU deposited produce a lesion with a probability of 0.1
    5. If success, set a lesion
    5. Emit again and repeat n times
"""
# Imports ##########################################################
from alinea.alep.vine import Vine
from random import random, seed
from alinea.alep.disease_operation import (generate_stock_du,
                                           generate_stock_lesions,
                                           generate_lesions_with_emission)
from alinea.alep.disease_outputs import (plot_lesions, 
                                         count_lesions,
                                         count_dispersal_units,
                                         count_lesions_by_leaf)
from alinea.alep.architecture import add_area_topvine, set_property_on_each_id, set_properties
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.alep.dispersal import PowderyMildewWindDispersal
from alinea.alep.protocol import disperse
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer
from openalea.vpltk import plugin

def save_image(scene, image_name='%s/img%04d.%s', directory='.', index=0, ext='png'):
    if not image_name:
        image_name='{directory}/img{index:0>4d}.{ext}'
    filename = image_name.format(directory=directory, index=index, ext=ext)
    Viewer.frameGL.saveImage(filename)
    return scene,
    
# Color leaves with dispersal units ################################
def update_plot(g):
    # Count lesions by id & add it as MTG property 
    nb_lesions_by_leaf = count_lesions_by_leaf(g)
    set_property_on_each_id(g, 'nb_lesions', nb_lesions_by_leaf, label = 'lf')
    
    # Visualization
    g = alep_colormap(g, 'nb_lesions', cmap=green_yellow_red(levels=100), 
                      lognorm=False, zero_to_one=False, vmax=75)
    labels = g.property('label')
    trunk_ids = [k for k,l in labels.iteritems() 
                if l.startswith('tronc') or l.startswith('en')]
    brown = (100,70,30)
    for id in trunk_ids:
        trunk = g.node(id)
        trunk.color = brown
    scene = plot3d(g)
    Viewer.display(scene)
    return scene
    
# Dispersal ########################################################
seed(2)

# Generate a fully developed vine canopy
vine = Vine()
g = vine.setup_canopy(age=50)
add_area_topvine(g, conversion_factor=1000., label = 'lf')

# List all leaves in the scene
labels = g.property('label')
leaves = [k for k,l in labels.iteritems() if l.startswith('lf')]

# Create a lesion and make it emit a stock of dispersal units
diseases=plugin.discover('alep.disease')
powdery_mildew = diseases['powdery_mildew'].load()
nb_lesions = 1
nb_dus = 150
lesions = generate_lesions_with_emission(nb_lesions, nb_dus, powdery_mildew)
  
# Deposit the initial lesion the source leaf
source = 921
g.node(source).lesions = lesions

# Define the wind direction
wind_direction = (1,0,0)
set_properties(g,label = 'lf', wind_direction=wind_direction)

# Call dispersal function
dispersor = PowderyMildewWindDispersal()
disperse(g, dispersor, "powdery_mildew", label='lf')

# Generate a lesion for each deposited DU
dispersal_units = g.property('dispersal_units')

# Display vine with affected leaves
update_plot(g)

# Repeat operation multiple times
for t in range(11):
    for vid, DU_list in dispersal_units.iteritems():
        for DU in DU_list:
            if random()<0.1:
                if not g.node(vid).lesions:
                    g.node(vid).lesions = []
                g.node(vid).lesions += generate_lesions_with_emission(nb_lesions, nb_dus, powdery_mildew)
    print(count_lesions(g))
    disperse(g, dispersor, "powdery_mildew", label='lf')

    # Generate a lesion for each deposited DU
    dispersal_units = g.property('dispersal_units')

    # Display vine with affected leaves
    update_plot(g)
    scene = plot3d(g)
    # save_image(scene, image_name='image%d.png' % t)