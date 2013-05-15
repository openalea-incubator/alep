""" Tutorial for septoria"""
import random as rd

from alinea.alep.wheat_examples import adel_mtg, adel_mtg2, adel_one_leaf
from alinea.alep.architecture import *
from alinea.alep.disease_operation import *
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.protocol import *

from openalea.vpltk import plugin

# Define a plant or canopy
g = adel_one_leaf()

# Add missing properties needed for the simulation
# The simulation requires the following properties on leaf elements:
#   - 'surface': total surface of the leaf element
#   - 'healthy_surface': surface of the leaf element without lesion or senescence
#   - 'position_senescence': position of the senescence on blade axis
set_properties(g,label = 'LeafElement',
               surface=5, healthy_surface=5, position_senescence=None)
               
# discover the disease implementation for septoriose
diseases=plugin.discover('alep.disease')
#septoria_classes = [kls for kls in diseases if kls.load().name == 'septoria']
# or 
septoria = diseases['septoria_exchanging_rings'].load()

# Create a pool of dispersal units (DU)
nb_du = 2
dispersal_units = generate_stock_du(nb_du, disease=septoria)

# Distribute the DU 
inoculator = RandomInoculation()
initiate(g, dispersal_units, inoculator)



