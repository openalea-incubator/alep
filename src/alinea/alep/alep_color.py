""" Adaptation of the colormap function from openalea.mtg to the needs of alep. """
# Imports #########################################################################
import numpy as np
from pylab import *
from matplotlib import cm
from matplotlib.colors import Normalize, LogNorm, ListedColormap, LinearSegmentedColormap

def green_yellow_red(levels=10):
    """ Generate a colormap from green to yellow then red.
    
    """
    colors =[(0., 0.5, 0.), 
             (1., 1., 0.), 
             (1, 0., 0.)]
    return LinearSegmentedColormap.from_list(name='green_yellow_red', colors=colors, N=levels)

def green_white(levels=10):
    """ Generate a colormap from green to white.
    
    """
    colors =[(0., 0.5, 0.),
             (1., 1., 1.)]
    return LinearSegmentedColormap.from_list(colors=colors, name='green_white', N=levels)

def green_lightblue_blue(levels=10):
    """ Generate a colormap from green to light blue then dark blue.
    """
    colors =[(0., 1., 0.), 
             (0., 1., 1.), 
             (0., 0., 1.)]
    return LinearSegmentedColormap.from_list(colors=colors, name='green_lightblue_blue', N=levels)
                                             
def alep_colormap(g, property_name, cmap='jet', lognorm=True, zero_to_one=True, vmax=None):
    """ Apply a colormap on a given MTG property to compute the 'color' property

    The colormap are thus defined in matplotlib.
    If lognorm is set to True, then the values are normalised on a log scale.
    """
    prop = g.property(property_name)

    keys = list(prop.keys())
    v = np.array(list(prop.values()))

    if type(cmap) is str:
        try:
            _cmap = cm.get_cmap(cmap)
        except:
            raise Exception('This colormap does not exist')
    else:
        _cmap = cmap
    
    if zero_to_one == True:
        norm = Normalize(vmin=0, vmax=1.) if not lognorm else LogNorm(vmin=0, vmax=1.)
    else:
        if vmax!=None:
            norm = Normalize(vmin=0, vmax=vmax) if not lognorm else LogNorm(vmin=0, vmax=vmax)
        else:
            norm = Normalize(vmin=0, vmax=max(v)) if not lognorm else LogNorm(vmin=0, vmax=max(v))
    values = norm(v)
    # values = v

    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    g.properties()['color'] = dict(list(zip(keys,colors)))
    return g
   
def plot_pesticide(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue, lognorm=False):
    """ plot the plant with pesticide doses """
    prop = g.property(property_name)
    keys = list(prop.keys())
    value = []
    for k, val in prop.items():
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
            n.properties()['color'] = dict(list(zip(keys,colors)))
        else : 
            n.color = green

    scene = plot3d(g)
    Viewer.display(scene)
    return g

def cmap_competition(N=None):
    """ Create a cmap from jet with a finite number of elements."""
    cmap = cm.jet
    if N==None:
        N = cmap.N
    cmaplist = [cmap(i) for i in range(N)]
    shuffle(cmaplist)
    # force the first color entry to be white
    cmaplist[0] = (1.0,1.0,1.0,1.0)
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, N)

    bounds = np.linspace(0,N,N+1)
    norm = mpl.colors.BoundaryNorm(bounds, N)
    return cmap, bounds, norm