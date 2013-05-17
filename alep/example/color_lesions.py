""" Example to use a colormap to visualize how much lesions are on the leaves of a MTG.

"""
# Imports #########################################################################
from alinea.alep.wheat_examples import adel_mtg2
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
from openalea.mtg import color

from alinea.alep.disease_operation import *
from alinea.alep.inoculation import RandomInoculation

# MTG generation ##################################################################
g = adel_mtg2()

# Lesions random distribution #####################################################
distribute_disease(g, fungal_object='lesion', 
                   nb_objects=100, disease_model='powdery_mildew', 
                   initiation_model=RandomInoculation())
                   
# Visualization ###################################################################
cmap = 'YlOrRd'
g = color.colormap(g, 'lesions', cmap=cmap, lognorm=False)