""" Architecture utilities

The aim of this module is to provide all the tools needed to adapt a given architecture 
into one that provides all the good parameters and organ name definition.

"""

def get_leaves(g, leaf_name='LeafElement'):
    labels = g.property('label')
    return [k for k,l in labels.iteritems() if l.startswith(leaf_name)]

def add_surface_topvine(g, conversion_factor=1000., label='lf'):
    """ Compute the surface of the leaves in topvine in cm2
    and add it as a MTG property.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    conversion_factor: float
        Conversion factor to pass from the default size unit of the MTG to cm2
    label: str
        Label of the part of the MTG concerned by the calculation
    """
    from openalea.plantgl import all as pgl
    
    geometries = g.property('geometry')
    surf = g.property('surface')
    h_surf = g.property('healthy_surface')
    if len(surf)==0:
        g.add_property('surface')
        surf = g.property('surface')
    if len(h_surf)==0:
        g.add_property('healthy_surface')
        h_surf = g.property('healthy_surface')
    new_vids = [n for n in g if g.label(n).startswith(label) if n not in surf]
    # Add surface information to each leaf
    surf.update({vid:pgl.surface(geometries[vid][0])*conversion_factor for vid in new_vids})
    # At start healthy surface equals total surface
    h_surf.update({vid:surf[vid] for vid in new_vids})
    
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
    vids = get_leaves(g,leaf_name=label)
    default_properties(g,vids,kwds)
    return g

def set_properties_on_new_leaves(g,
                                 label = 'LeafElement',
                                 **kwds):
    """ Give initial properties to newly emerged leaves.
    """
    vids = get_leaves(g,leaf_name=label)
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
    property_dict = dict({(k,property_dict.get(k, 0.)) for k in vids})
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
    property_dict = dict({(k,property_dict.get(k, 0.)) for k in vids})
    g.add_property(property_name)
    prop = g.property(property_name)
    prop.update(property_dict)