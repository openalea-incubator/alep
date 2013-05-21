""" Adaptation of the colormap function from openalea.mtg to the needs of alep. """
# Imports #########################################################################
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize, LogNorm, LinearSegmentedColormap

def green_yellow_red(levels=10):
    """ Generate a colormap from green to yellow then red.
    
    """
    return LinearSegmentedColormap.from_list(name='green_yellow_red', 
                                             colors =[(0., 1., 0.), 
                                                      (1., 1., 0.), 
                                                      (1, 0., 0.)],
                                             N=levels)

def alep_colormap(g, property_name, cmap='jet',lognorm=True):
    """ Apply a colormap on a given MTG property to compute the 'color' property

    The colormap are thus defined in matplotlib.
    If lognorm is set to True, then the values are normalised on a log scale.
    """
    prop = g.property(property_name)

    keys = prop.keys()
    v = np.array(prop.values())

    if type(cmap) is str:
        try:
            _cmap = cm.get_cmap(cmap)
        except:
            raise Exception('This colormap does not exist')
    else:
        _cmap = cmap
            
    norm = Normalize() if not lognorm else LogNorm() 
    values = norm(v)

    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    g.properties()['color'] = dict(zip(keys,colors))
    return g
