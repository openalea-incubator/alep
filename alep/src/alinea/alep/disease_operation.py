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
                        
def generate_stock_lesions(nb_lesions, disease, position=None):
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
    return [lesion(nb_spores=rd.randint(1,100), position=position) for i in range(nb_lesions)]

def generate_lesions_with_emission(nb_lesions, nb_dus, disease):
    """ Generate a stock of lesions that are emitting dispersal units.
    
    Extend the class 'emission' of the lesion in order to force the emission
    of a list of dispersal units, independently from its surface, age, status,
    stock, etc.
    
    Parameters
    ----------
    nb_lesions: int
        Number of lesion to create in the stock
    nb_dus: int
        Number of dispersal units emitted by each lesion
    disease: model
        Implementation for a model of fungal disease
    
    Returns
    -------
    lesions: list of objects
        List of lesions of the given disease that are emitting dispersal units    
    """
    # Override the method emission of the disease
    LesionKlass = disease.lesion()
    def emission_forced(lesion, leaf=None):
        """ Create a list of dispersal units emitted by the lesion.
        
        Force the emission of a given number of dispersal unit by the lesion
        independently from its surface, age, status, stock, etc.        
        """
        # For dispersal, stock of spores needs to be positive
        # TODO : change this condition to a method 'can_emit'
        lesion.is_active = True
        lesion.stock_spores = 1
        DU = disease.dispersal_unit()
        return [DU(nb_spores=rd.randint(1,100), status='emitted',
                    position=lesion.position) for i in range(nb_dus)]
    LesionKlass.emission = emission_forced   
    lesions = [LesionKlass(nb_spores=1) for i in range(nb_lesions)]
    return lesions
    
def override_method(LesionKlass, method=None):
    """ Override a specific method in a lesion class.
    
    Parameters
    ----------
    LesionKlass: model
        Class for the model of lesion of a given disease
    method: 
        Method to be overridden
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
    diseases = plugin.discover('alep.disease')
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