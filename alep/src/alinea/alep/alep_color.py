""" Adaptation of the colormap function from openalea.mtg to the needs of alep. """
# Imports #########################################################################
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize, LogNorm, ListedColormap

def green_yellow_red(levels=10):
    """ Generate a colormap from green to yellow then red.
    
    """
    colors =[(0., 0.5, 0.), 
             (1., 1., 0.), 
             (1, 0., 0.),
             (1, 0., 0.)]
    return ListedColormap(colors=colors, name='green_yellow_red', N=levels)

def green_white(levels=10):
    """ Generate a colormap from green to white.
    
    """
    colors =[(0., 0.5, 0.),
             (1., 1., 1.)]
    return ListedColormap(colors=colors, name='green_white', N=levels)

def green_lightblue_blue(levels=10):
    """ Generate a colormap from green to light blue then dark blue.
    """
    colors =[(0., 1., 0.), 
             (0., 1., 1.), 
             (0., 0., 1.)]
    return ListedColormap(colors=colors, name='green_lightblue_blue', N=levels)
                                             
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
    # values = v

    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    g.properties()['color'] = dict(zip(keys,colors))
    return g

    
def plot_pesticide(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue, lognorm=False):
    """ plot the plant with pesticide doses """
    prop = g.property(property_name)
    keys = prop.keys()
    value = []
    for k, val in prop.iteritems():
        value.append(val[compound_name])
        v = np.array(value)

    if type(cmap) is str:
        try:
            _cmap = cm.get_cmap(cmap())
        except:
            raise Exception('This colormap does not exist')
    else:
        _cmap = cmap()

    green = (0,180,0)

    norm = Normalize(vmin=0, vmax=max(v)) if not lognorm else LogNorm(vmin=0, vmax=max(v)) 
    values = norm(v)

    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    for v in g.vertices(scale=g.max_scale()): 
        n = g.node(v)
        if 'surfacic_doses' in n.properties():
            n.properties()['color'] = dict(zip(keys,colors))
        else : 
            n.color = green

    scene = plot3d(g)
    Viewer.display(scene)
    return g

