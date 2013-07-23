""" Example to use a colormap to visualize how much lesions are on the leaves of a MTG.

"""
# Imports #########################################################################
from alinea.alep.wheat import adel_mtg2
from alinea.adel.mtg_interpreter import *
from alinea.alep.alep_color import alep_colormap, green_yellow_red, green_white
from openalea.plantgl.all import *
# from openalea.mtg import color

from alinea.alep.disease_operation import *
from alinea.alep.disease_outputs import *
from alinea.alep.architecture import *
from alinea.alep.inoculation import RandomInoculation

# MTG generation ##################################################################
g = adel_mtg2()

# Lesions random distribution #####################################################
# distribute_disease(g, fungal_object='lesion', 
                   # nb_objects=500, disease_model='septoria_exchanging_rings', 
                   # initiation_model=RandomInoculation())
distribute_disease(g, fungal_object='lesion', 
                   nb_objects=500, disease_model='powdery_mildew', 
                   initiation_model=RandomInoculation())

# Count lesions by id & add it as MTG property ####################################
nb_lesions_by_leaf = count_lesions_by_leaf(g, label = 'LeafElement')
set_property_on_each_id(g, 'nb_lesions', nb_lesions_by_leaf, label = 'LeafElement')
                   
# Visualization ###################################################################
cmap = 'jet'
g = alep_colormap(g, 'nb_lesions', cmap=green_white(levels=10), lognorm=False)
scene = plot3d(g)
Viewer.display(scene)
    
