# Wheat ###############################################################
from alinea.alep.wheat import *

def wheat(wheat_model='adel_mtg2'):
    """ Create a MTG according to given model.
    
    Parameters
    ----------
    wheat_model: model
        Model of wheat among : adel_mtg(), adel_one_leaf(), adel_mtg2()
        
    Returns
    -------
    g: MTG
        Wheat MTG with properties
    """
    return eval(wheat_model + '()')

# Scene ###############################################################
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

def scene_from_g(g):
    """ Generate a scene from a MTG.
    
    Use the method plot3d.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    
    Returns
    -------
    scene: scene
        Scene with the MTG
    """
    return plot3d(g)
    
# Set properties on g #################################################
from alinea.alep.architecture import set_properties

def set_properties_node(g, label='LeafElement', dict=None):
    """ Adapt the function set_properties to receive variable inputs.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG where new properties will be added
    dict: dict
        Dictionnary with variable inputs
    
    Returns
    -------
    g: MTG
        Updated MTG with given properties
    """
    return set_properties(g, label=label, **dict)
    
# Distribute dispersal units on g #####################################
from alinea.alep.disease_operation import *
from alinea.alep.protocol import *
from alinea.alep.inoculation import RandomInoculation
from openalea.vpltk import plugin

def distribute_dispersal_units(g, nb_du=1, disease_model='powdery_mildew'):
    """ Distribute dispersal units on the MTG.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    nb_du: int
        Number of dispersal units to distribute on the mtg
    disease_model: model
        Type of model to compute disease lesion development
        
    Returns
    -------
    g: MTG
        Updated MTG with dispersal units
    """
    # Create a pool of dispersal units (DU)
    diseases=plugin.discover('alep.disease')
    disease = diseases[disease_model].load()
    dispersal_units = generate_stock_du(nb_du=nb_du, disease=disease)
    # Distribute the DU 
    inoculator = RandomInoculation()
    initiate(g, dispersal_units, inoculator)
    return g
    
