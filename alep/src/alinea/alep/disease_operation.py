""" Disease operation utilities

The aim of this module is to provide all the tools needed to manipulate
fungal objects on the MTG.

"""

import random as rd
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

def generate_stock_du(nb_du, disease):
    """ Generate a stock of dispersal units.
    
    Parameters
    ----------
    disease: model
        Implementation for a model of fungal disease
        
    Returns
    -------
    dus: list of objects
        List of dispersal units of the given disease
    """
    DU = disease.dispersal_unit()
    return [DU(nb_spores=rd.randint(1,100), status='emitted')
                        for i in range(nb_du)]
                        
def generate_stock_lesion(nb_du, disease):
    """ Generate a stock of lesions.
    
    Parameters
    ----------
    disease: model
        Implementation for a model of fungal disease
        
    Returns
    -------
    lesions: list of objects
        List of lesions of the given disease
    """
    lesion = disease.dispersal_unit()
    return [lesion(nb_spores=rd.randint(1,100)) for i in range(nb_du)]
    
def plot_lesions(g):
    """ plot the plant with infected elements in red """
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)