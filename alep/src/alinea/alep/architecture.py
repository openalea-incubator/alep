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

def set_properties_on_new_leaves(g,
                                 label = 'LeafElement',
                                 **kwds):
    """ Give initial properties to newly emerged leaves.
    """
    vids = leaves(g,leaf_name=label)
    for name, default in kwds.iteritems():
        prop = g.property(name)
        new_vids = [n for n in g if g.label(n).startswith(label) if n not in prop]
        if default is not None:
            prop.update({vid:default for vid in new_vids})
    
def set_property_on_leaves(g,
                            property_name,
                            property_dict,
                            label = 'LeafElement'):
    """
    """ 
    vids = [n for n in g if g.label(n).startswith(label)]
    property_dict = dict({(k,property_dict.get(k, 0)) for k in vids})
    g.add_property(property_name)
    prop = g.property(property_name)
    prop.update(property_dict)
    
def set_property_on_each_id(g,
                            property_name,
                            property_dict,
                            label = 'LeafElement'):
    """
    """ 
    vids = [n for n in g]
    property_dict = dict({(k,property_dict.get(k, 0)) for k in vids})
    g.add_property(property_name)
    prop = g.property(property_name)
    prop.update(property_dict)