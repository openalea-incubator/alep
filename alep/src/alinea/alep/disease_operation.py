""" Disease operation utilities

The aim of this module is to provide all the tools needed to manipulate
fungal objects on the MTG.

"""

import random as rd
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

from alinea.alep.protocol import *
from alinea.alep.inoculation import RandomInoculation
from openalea.vpltk import plugin

def generate_stock_du(nb_dus, disease):
    """ Generate a stock of dispersal units.
    
    Parameters
    ----------
    nb_du: int
        Number of dispersal units to create in the stock
    disease: model
        Implementation for a model of fungal disease
        
    Returns
    -------
    dus: list of objects
        List of dispersal units of the given disease
    """
    DU = disease.dispersal_unit()
    return [DU(nb_spores=rd.randint(1,100), status='emitted')
                        for i in range(nb_dus)]
                        
def generate_stock_lesions(nb_lesions, disease):
    """ Generate a stock of lesions.
    
    Parameters
    ----------
    nb_lesions: int
        Number of lesion to create in the stock
    disease: model
        Implementation for a model of fungal disease
        
    Returns
    -------
    lesions: list of objects
        List of lesions of the given disease
    """
    lesion = disease.lesion()
    return [lesion(nb_spores=rd.randint(1,100)) for i in range(nb_lesions)]
    
def generate_lesions_with_emission(nb_lesions, nb_dus, disease):
    """ Generate a stock of lesions.
    
    Parameters
    ----------
    nb_lesions: int
        Number of lesions to create in the stock
    nb_dus: int
        Number of dispersal units emitted by each lesion
    disease: model
        Implementation for a model of fungal disease
        
    Returns
    -------
    lesions: list of objects
        List of lesions of the given disease
    """
    lesions = generate_stock_lesions(nb_lesions, disease)
    for l in lesions:
        l.stock_spores = nb_dus
    LesionKlass = disease.lesion()
    # Overwrite the function of emission of the chosen lesion model
    def simple_emission(lesion, leaf=None):
        """ Set the emission of a lesion to a given number of dispersal units
        independently from its surface, age, status, stock, etc.
        """
        return generate_stock_du(nb_dus, disease)
    LesionKlass.emission = simple_emission
    return lesions

def overwrite_method(LesionKlass, method=None):
    """ Overwrite a specific method in a lesion class.
    
    Parameters
    ----------
    LesionKlass: model
        Class for the model of lesion of a given disease
    method: 
        Method to be overwited
    """
    if method:
        LesionKlass.method = method
    
def distribute_disease(g,
                       fungal_object='lesion', 
                       nb_objects=1, 
                       disease_model='powdery_mildew',
                       initiation_model=RandomInoculation(),
                       label="LeafElement"):
    """ Distribute fungal objects on the MTG.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    fungal_object: str
        Type of fungal object. Choose between : 'dispersal_unit' or 'lesion'
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
        Updated MTG with dispersal units or lesions
    """
    # Create a pool of dispersal units (DU)
    diseases=plugin.discover('alep.disease')
    disease = diseases[disease_model].load()
    if fungal_object=='dispersal_unit':
        objects = generate_stock_du(nb_dus=nb_objects, disease=disease)
    elif fungal_object=='lesion':
        objects = generate_stock_lesions(nb_lesions=nb_objects, disease=disease)
    else:
        raise Exception('fungal object is not valid: choose between ''dispersal_unit'' or ''lesion')
    # Distribute the DU 
    initiate(g, objects, initiation_model, label=label)
    return g
    
def distribute_dispersal_units(g, nb_dus=1, 
                               disease_model='septoria_exchanging_rings',
                               initiation_model=RandomInoculation(),
                               label="LeafElement"):
    """ Distribute dispersal units on the MTG.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    nb_dus: int
        Number of dispersal units to distribute on the MTG
    disease_model: model
        Type of model to compute disease lesion development
    initiation_model: model
        Model that sets the position of each DU/lesion in stock on g
        Requires a method named 'allocate' (see doc)
        
    Returns
    -------
    g: MTG
        Updated MTG with dispersal units or lesions
    """
    distribute_disease(g,
                       fungal_object='dispersal_unit', 
                       nb_objects=nb_dus, 
                       disease_model=disease_model,
                       initiation_model=initiation_model,
                       label=label)
                       
def distribute_lesions(g, nb_lesions=1, 
                       disease_model='septoria_exchanging_rings',
                       initiation_model=RandomInoculation(),
                       label="LeafElement"):
    """ Distribute lesions on the MTG.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    nb_lesions: int
        Number of lesions to distribute on the MTG
    disease_model: model
        Type of model to compute disease lesion development
    initiation_model: model
        Model that sets the position of each DU/lesion in stock on g
        Requires a method named 'allocate' (see doc)
        
    Returns
    -------
    g: MTG
        Updated MTG with dispersal units or lesions
    """
    distribute_disease(g,
                       fungal_object='lesion', 
                       nb_objects=nb_lesions, 
                       disease_model=disease_model,
                       initiation_model=initiation_model,
                       label=label)