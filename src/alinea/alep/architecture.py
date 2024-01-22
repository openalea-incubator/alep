""" Architecture utilities

The aim of this module is to provide all the tools needed to adapt a given architecture 
into one that provides all the good parameters and organ name definition.

"""

def get_leaves(g, label='LeafElement'):
    labels = g.property('label')
    return [k for k,l in labels.items() if l.startswith(label)]

def get_total_leaf_area(g, label='LeafElement'):
    leaves = get_leaves(g)
    return sum(g.node(leaf).area for leaf in leaves)
    
def add_area_topvine(g, conversion_factor=1000., label='lf'):
    """ Compute the area of the leaves in topvine in cm2
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
    area = g.property('area')
    g_area = g.property('green_area')
    if len(area)==0:
        g.add_property('area')
        area = g.property('area')
    if len(g_area)==0:
        g.add_property('green_area')
        g_area = g.property('green_area')
    new_vids = [n for n in g if g.label(n).startswith(label) if n not in area]
    # Add area information to each leaf
    area.update({vid:pgl.surface(geometries[vid][0])*conversion_factor for vid in new_vids})
    # At start green area equals total area
    g_area.update({vid:area[vid] for vid in new_vids})
    
def default_properties(g,vids,props):
    for name, default in props.items():
        g.add_property(name)
        prop = g.property(name)
        if default is not None:
            prop.update({vid:default for vid in vids})
    return props
    
def set_properties(g,
                   label = 'LeafElement',
                   **kwds):
    """ Give values for plant properties of each LeafElement. 
    
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
    vids = get_leaves(g, label=label)
    default_properties(g,vids,kwds)
    return g

def update_healthy_area(g, label = 'LeafElement'):
    """ Update values for healthy area to each LeafElement. 
    
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
    healthy_areas: dict([vid:healthy_area])
        Dictionary of healthy area (in cm2) of each vid
    """
    from alinea.alep.disease_outputs import compute_healthy_area_by_leaf
    from alinea.alep.architecture import set_properties
    healthy_areas = g.property('healthy_area')
    if len(healthy_areas)==0:
        g.add_property('healthy_area')
        healthy_areas = g.property('healthy_area')
    healthy_areas.update(compute_healthy_area_by_leaf(g, label))
    return g
    
def set_properties_on_new_leaves(g,
                                 label = 'LeafElement',
                                 **kwds):
    """ Give initial properties to newly emerged leaves.
    """
    vids = get_leaves(g, label=label)
    for name, default in kwds.items():
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