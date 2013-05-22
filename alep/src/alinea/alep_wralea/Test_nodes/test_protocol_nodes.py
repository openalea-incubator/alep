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
    
# Operations on fungal objects #######################################
from alinea.alep.disease_operation import distribute_disease
from alinea.alep.inoculation import RandomInoculation

def distribute_dispersal_units(g, nb_objects=1, 
                               disease_model='powdery_mildew',
                               initiation_model=RandomInoculation()):
    """ Use the method 'distribute_disease' to distribute dispersal units on g. 
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    nb_objects: int
        Number of dispersal units or lesions to distribute on the MTG
    disease_model: model
        Type of model to compute disease lesion development
    initiation_model: model
        Model that sets the position of each DU/lesion in stock on g
        Requires a method named 'allocate' (see doc)
        
    Returns
    -------
    g: MTG
        Updated MTG with dispersal units
    """
    distribute_disease(g, fungal_object='dispersal_unit', 
                       nb_objects=nb_objects, 
                       disease_model=disease_model,
                       initiation_model=initiation_model)
                               
def distribute_lesions(g, nb_objects=1, 
                       disease_model='powdery_mildew',
                       initiation_model=RandomInoculation()):
    """ Use the method 'distribute_disease' to distribute lesions on g. 
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    nb_objects: int
        Number of dispersal units or lesions to distribute on the MTG
    disease_model: model
        Type of model to compute disease lesion development
    initiation_model: model
        Model that sets the position of each DU/lesion in stock on g
        Requires a method named 'allocate' (see doc)
        
    Returns
    -------
    g: MTG
        Updated MTG with lesions
    """
    distribute_disease(g, fungal_object='lesion', 
                       nb_objects=nb_objects, 
                       disease_model=disease_model,
                       initiation_model=initiation_model)                            

    
