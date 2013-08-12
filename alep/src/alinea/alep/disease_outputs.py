""" Utilities for computing outputs from the disease model.

The aim of this module is to provide all the tools needed to compute
the outputs of the disease models. 
"""

def count_lesions(g):
    """ Count lesions of the mtg.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    nb_lesions: int
        Number of lesions on the MTG
    """
    lesions = g.property('lesions')
    return sum(len(l) for l in lesions.itervalues())

def count_dispersal_units(g):
    """ Count dispersal units of the mtg.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    nb_dispersal_units: int
        Number of dispersal units on the MTG
    """
    dispersal_units = g.property('dispersal_units')
    return sum(len(l) for l in dispersal_units.itervalues())
    
def count_lesions_by_leaf(g):
    """ Count lesions on each leaf of the MTG.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    nb_lesions_by_leaf: dict([id:nb_lesions])
        Number of lesions on each part of the MTG given by the label
    """
    lesions = g.property('lesions')
    return {k:len(v) for k,v in lesions.iteritems()}

def count_dispersal_units_by_leaf(g, label='LeafElement'):
    """ Count dispersal units on each part of the MTG given by the label.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    nb_dispersal_units_by_leaf: dict([id:nb_dispersal_units])
        Number of dispersal units on each part of the MTG given by the label
    """
    dispersal_units = g.property('dispersal_units')
    return {k:len(v) for k,v in dispersal_units.iteritems()}

def plot_lesions(g):
    """ plot the plant with infected elements in red """
    from alinea.adel.mtg_interpreter import plot3d
    from openalea.plantgl.all import Viewer
    
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

def plot_dispersal_units(g):
    """ plot the plant with infected elements in red """
    from alinea.adel.mtg_interpreter import plot3d
    from openalea.plantgl.all import Viewer
    
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'dispersal_units' in n.properties():
            n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)

def compute_total_severity(g):
    """ Compute disease severity on the whole plant.
    
    Severity is the ratio between disease surface and green leaf area (in %).
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    severity: float
        Ratio between disease surface and green leaf area (in %)
    """
    # Initiation
    green_areas = g.property('green_area')
    healthy_areas = g.property('healthy_area')
    
    total_green_area = sum(green_areas.values())
    total_disease_area = total_green_area - sum(healthy_areas.values())

    # Compute ratio, i.e. severity
    return 100 * total_disease_area / total_green_area if total_green_area > 0. else 0.

def compute_lesion_areas_by_leaf(g, label='LeafElement'):
    """ Compute healthy area on each part of the MTG given by the label.
    
    Healthy area is green area (without senescence) minus the surface of lesions.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    lesion_surfaces_by_leaf: dict([id:surface_lesions])
        Surface of the lesions on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    lesions = g.property('lesions')
    return {vid:(sum(l.surface for l in lesions[vid] if l.is_active)
            if vid in lesions.keys() else 0.) for vid in vids} 

def compute_healthy_area_by_leaf(g, label='LeafElement'):
    """ Compute healthy area on each part of the MTG given by the label.
    
    Healthy area is green area (without senescence) minus the surface of lesions.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    lesion_surfaces_by_leaf: dict([id:surface_lesions])
        Surface of the lesions on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    green_areas = g.property('green_area')
    lesion_areas = compute_lesion_areas_by_leaf(g, label)
    return {vid:(green_areas[vid] - lesion_areas[vid] 
            if round(green_areas[vid], 10)>round(lesion_areas[vid], 10) else 0.)
            for vid in vids}
    
def compute_severity_by_leaf(g, label='LeafElement'):
    """ Compute severity of the disease on each part of the MTG given by the label.
    
    Severity is the ratio between disease surface and green leaf area (in %).
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    severity_by_leaf: dict([id:nb_lesions])
        Severity on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    green_areas = g.property('green_area')
    healthy_areas = g.property('healthy_area')
    return {vid:(100*(1-(healthy_areas[vid]/float(green_areas[vid]))) if green_areas[vid]>0. else 0.) for vid in vids}

def compute_total_necrosis(g):
    """ Compute necrosis percentage on the whole plant.
    
    Necrosis percentage ratio between necrotic (and sporulating) disease surface and total area of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    necrosis_percentage: float
        Ratio between necrotic (and sporulating) disease area and total area of leaves (in %)
    """
    # Leaf
    green_areas = g.property('green_area')
    total_green_area = sum(green_areas.values())
    
    # Disease
    lesions = g.property('lesions')
    if lesions:
        lesions = [l for les in lesions.values() for l in les 
                    if (l.status>=l.fungus.NECROTIC and not l.is_senescent)]
        total_necrotic_area = sum(l.surface for l in lesions)
    else:
        total_necrotic_area = 0.

    # Compute ratio, i.e. necrosis percentage
    return 100 * total_necrotic_area / total_green_area if total_green_area > 0. else 0.

######################################################################
class LeafInspector:
    def __init__(self, g, leaf_number=1, label='LeafElement'):
        """ Find the ids of the leaf elements on the chosen blade of the main stem.
        
        First leaf is the upper one. This method computes the id of the 'blade'.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        leaf_number: int
            Number of the chosen leaf
        label: str
            Label of the part of the MTG concerned by the calculation
        """
        labels = g.property('label')
        mets = [n for n in g if g.label(n).startswith('metamer') and g.order(n)==0]
        bids = [v for v,l in labels.iteritems() if l.startswith('blade')]
        blade_id = bids[len(mets)-leaf_number]
        # Find leaf elements on the blade
        self.leaf_elements = [id for id in g.components(blade_id) if labels[id].startswith(label)]
        # Initialize ratios (surfaces in state compared to green surface on leaf)
        self.ratio_inc = []
        self.ratio_chlo = []
        self.ratio_nec = []
        self.ratio_spo = []
        # Initialize total severity
        self.severity = []
        # Initialize necrosis percentage
        self.necrosis = []
    
    def compute_states(self, g):
        """ Compute surface of lesions in chosen state on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        # Compute total area on the blade
        total_green_area = compute_green_area(g)
        
        # Compute disease surface
        lesions = g.property('lesions')
        les = [l for id in self.leaf_elements for l in lesions[id] if not l.is_senescent]
                        
        if les:
            for l in les:
                l.compute_all_surfaces()
            
            surface_inc = sum(l.surface_inc for l in les)
            surface_chlo = sum(l.surface_chlo for l in les)
            surface_nec = sum(l.surface_nec for l in les)
            surface_spo = sum(l.surface_spo for l in les)
        else:
            surface_inc = 0.
            surface_chlo = 0.
            surface_nec = 0.
            surface_spo = 0.

        self.ratio_inc.append(100 * surface_inc / total_green_area if total_green_area>0. else 0.)
        self.ratio_chlo.append(100 * surface_chlo / total_green_area if total_green_area>0. else 0.)
        self.ratio_nec.append(100 * surface_nec / total_green_area if total_green_area>0. else 0.)
        self.ratio_spo.append(100 * surface_spo / total_green_area if total_green_area>0. else 0.)
        assert (round((surface_inc + surface_chlo + surface_nec + surface_spo),4) == 
                round(sum(l.surface for l in les), 4))
    
    def compute_severity(self, g):
        """ Compute severity on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        # Compute green area on the blade
        total_green_area = compute_green_area(g)
        # Compute disease area on the blade
        total_disease_area = total_green_area - compute_healthy_area(g)
        # Compute severity
        self.severity.append(100 * total_disease_area / total_green_area if total_green_area>0. else 0.)       
    
    def compute_healthy_area(self, g):
        """ Compute healthy area on the leaf.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
           
        Returns
        -------
            Leaf healthy area
        """
        healthy_areas = g.property('healthy_area')
        return sum(healthy_areas[id] for id in leaf_elements)
    
    def compute_green_area(self, g):
        """ Compute healthy area on the leaf.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
           
        Returns
        -------
            Leaf healthy area
        """
        green_areas = g.property('healthy_area')
        return sum(green_areas[id] for id in leaf_elements)
    
    def count_dispersal_units(self, g):
        """ count DU of the leaf.
   
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
           
        Returns
        -------
            Number of dus on the leaf
        """
        dispersal_units = g.property('dispersal_units')
        return sum(1 for vid in self.leaf_elements for du in dispersal_units[vid] if du.is_active)
        
    def count_DU_on_healthy(self, g, nb_unwashed):
        """ Count DU of the leaf that are on healthy tissue.
        
        Same calculation as in BiotrophProbaModel. 
        Might not be the actual number of DUs on healthy tissue.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        nb_unwashed: int
            Number of DUs staying after washing
           
        Returns
        -------
            Number of dus on healthy tissue on the leaf
        """
        severity = self.compute_severity(g)/100
        return round(severity * nb_unwashed)