""" Architecture utilities

The aim of this module is to provide all the tools needed to adapt a given architecture 
into one that provides all the good parameters and organ name definition.

"""

def leaves(g, leaf_name='LeafElement'):
    return [n for n in g if g.label(n).startswith(leaf_name)]

def default_properties(g,vids,props):
    for name, default in props.iteritems():
        g.add_property(name)
        prop = g.property(name)
        if default is not None:
            prop.update({vid:default for vid in vids})
    return props
    
def set_properties(g,
                   label = 'LeafElement',
                   **kwds):
    """ Give initial values for plant properties of each LeafElement. 
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    surface: float
        Initial surface of each leaf element
    position_senescence: float
        Position of senescence on blade axis
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy
    """
    vids = leaves(g,leaf_name=label)
    default_properties(g,vids,kwds)
        
    return g
