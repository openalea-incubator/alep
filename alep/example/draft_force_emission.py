""" Draft illustrating the extension of the class of septoria with a forced emission.
"""
# Imports ##########################################################
from alinea.adel.astk_interface import AdelWheat
from alinea.astk.plant_interface import *
from alinea.alep.protocol import *
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.disease_operation import generate_lesions_with_emission
from alinea.alep.disease_outputs import count_dispersal_units_by_leaf, plot_lesions
from alinea.alep.architecture import set_properties,set_property_on_each_id
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer
from openalea.vpltk import plugin
from time import sleep

# Plot function ####################################################
def update_plot(g, source):      
    # Count dispersal units by id & add it as MTG property 
    nb_dus_by_leaf = count_dispersal_units_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'nb_dus', nb_dus_by_leaf, label = 'LeafElement')

    # Visualization
    g = alep_colormap(g, 'nb_dus', cmap=green_yellow_red(levels=10), lognorm=False, zero_to_one=False)
    blue = (0, 0, 180)
    g.node(source).color = blue
    scene = plot3d(g)
    Viewer.display(scene)
    return scene

# Dispersal ########################################################
# Generate a fully developed wheat plant
wheat = AdelWheat()
g,_ = new_canopy(wheat,age=1000)

# Select the source leaf
source = 50

# Define wind direction
rain_intensity = 3
set_properties(g,label = 'LeafElement', rain_intensity=3)

# Make it emit a stock of dispersal units
diseases = plugin.discover('alep.disease')
septoria = diseases['septoria_exchanging_rings'].load()
g.node(source).lesions  = generate_lesions_with_emission(nb_lesions=1, nb_dus=10000, disease=septoria)
update_plot(g, source)

# Call dispersal function
dispersor = Septo3DSplash()
disperse(g, dispersor, "septoria", label='LeafElement')

# Display vine with affected leaves
sleep(2)
update_plot(g, source)
