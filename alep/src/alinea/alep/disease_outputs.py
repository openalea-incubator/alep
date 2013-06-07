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
    
def count_lesions_by_leaf(g, label='LeafElement'):
    """ Count lesions on each part of the MTG given by the label.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    nb_lesions_by_leaf: dict([id:nb_lesions])
        Number of lesions on each part of the MTG given by the label
    """
    lesions = g.property('lesions')
    return {k:len(v) for k,v in lesions.iteritems()}

def count_lesion_surfaces_by_leaf(g, label='LeafElement'):
    """ Count the surface of lesions on each part of the MTG given by the label.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    surface_lesions_by_leaf: dict([id:nb_lesions])
        Surface of lesions on each part of the MTG given by the label
    """
    lesions = g.property('lesions')
    return {k:sum(l.surface for l in v) for k,v in lesions.iteritems()}
    
def compute_severity_by_leaf(g, label='LeafElement'):
    """ Compute severity of the disease on each part of the MTG given by the label.
    
    Severity is the ratio between disease surface and total leaf surface.
    
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
    lesions = g.property('lesions')
    surfaces = g.property('surface')
    return {k:(sum(l.surface for l in v)/float(surfaces[k]) if surfaces[k]>0 else 0.) for k,v in lesions.iteritems()}
    
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

def compute_total_severity(g):
    """ Compute disease severity on the whole plant.
    
    Severity is the ratio between disease surface and total surface of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    severity: float
        Ratio between disease surface and total surface of leaves (in %)
    """
    # Initiation
    surfaces = g.property('surface')
    healthy_surfaces = g.property('healthy_surface')
    
    total_leaf_surface = sum(surfaces.values())
    total_disease_surface = total_leaf_surface - sum(healthy_surfaces.values())

    # Compute ratio, i.e. severity
    severity = 100 * total_disease_surface / total_leaf_surface if total_leaf_surface > 0. else 0.
    
    return severity

def compute_total_necrosis(g):
    """ Compute necrosis percentage on the whole plant.
    
    Necrosis percentage ratio between necrotic (and sporulating) disease surface and total surface of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    necrosis_percentage: float
        Ratio between necrotic (and sporulating) disease surface and total surface of leaves (in %)
    """
    # Leaf
    surfaces = g.property('surface')
    total_leaf_surface = sum(surfaces.values())
    
    # Disease
    lesions = g.property('lesions')
    if lesions:
        lesions = [l for les in lesions.values() for l in les if l.status>=l.fungus.NECROTIC]
        total_necrotic_surface = sum(l.surface for l in lesions)
    else:
        total_necrotic_surface = 0.

    # Compute ratio, i.e. severity
    necrosis_percentage = 100 * total_necrotic_surface / total_leaf_surface if total_leaf_surface > 0. else 0.
    
    return necrosis_percentage
    
######################################################################
class LeafInspector:
    def __init__(self, g, leaf_number=1):
        """ Find the ids of the leaf elements on the chosen blade of the main stem.
        
        First leaf is the upper one. This method returns the id of the 'blade'
        carrying 'LeafElements'.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        leaf_number: int
            Number of the chosen leaf
        """
        labels = g.property('label')
        surfaces = g.property('surface')
        mets = [n for n in g if g.label(n).startswith('metamer') and g.order(n)==0]
        bids = [v for v,l in labels.iteritems() if l.startswith('blade')]
        blade_id = bids[len(mets)-leaf_number]
        # Find leaf elements on the blade
        self.leaf_elements = [id for id in g.components(blade_id) if labels[id].startswith('LeafElement')]
        # Initialize ratios (surfaces in state compared to green surface on leaf)
        self.ratio_inc = []
        self.ratio_chlo = []
        self.ratio_nec = []
        self.ratio_spo = []
        # Initialize total severity
        self.severity = []
        # Initialize necrosis percentage
        self.necrosis_percentage = []
    
    def compute_states(self, g):
        """ Compute surface of lesions in chosen state on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        surfaces = g.property('surface')
        leaf_elements = self.leaf_elements
        
        # Compute total surface on the blade
        total_surface = sum(surfaces[id] for id in leaf_elements)
        
        # Compute disease surface
        lesions = g.property('lesions')
        les=[]
        for id in leaf_elements:
            if id in lesions.keys():
                les += lesions[id]

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

        self.ratio_inc.append(100 * surface_inc / total_surface if total_surface>0. else 0.)
        self.ratio_chlo.append(100 * surface_chlo / total_surface if total_surface>0. else 0.)
        self.ratio_nec.append(100 * surface_nec / total_surface if total_surface>0. else 0.)
        self.ratio_spo.append(100 * surface_spo / total_surface if total_surface>0. else 0.)
    
    def compute_severity(self, g):
        """ Compute severity on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        surfaces = g.property('surface')
        leaf_elements = self.leaf_elements
        
        # Compute green surface on the blade
        total_surface = sum(surfaces[id] for id in leaf_elements)
        # Compute disease surface
        lesions = g.property('lesions')
        disease_surface = sum([l.surface for id in leaf_elements for l in lesions[id]])

        self.severity.append(100 * disease_surface / total_surface if total_surface>0. else 0.)
    
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
        leaf_elements = self.leaf_elements
        return sum(1 for id in leaf_elements for d in dispersal_units[id] if d.is_active)
        
    def count_DU_on_green(self, g, nb_unwashed):
        """ Count DU of the leaf that are on green tissue.
        
        Same calculation as in BiotrophProbaModel. 
        Might not be the actual number of DUs on green tissue.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        nb_unwashed: int
            Number of DUs staying after washing
           
        Returns
        -------
            Number of dus on green tissue on the leaf
        """
        leaf_elements = self.leaf_elements
        healthy = self.compute_healthy_surface(g)
        surface = sum(g.node(vid).surface for vid in leaf_elements)
        ratio = healthy / surface if (surface>0. and healthy>0.) else 0.
        
        return round(ratio * nb_unwashed)
    
    def compute_healthy_surface(self, g):
        """ Compute healthy surface on the leaf.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
           
        Returns
        -------
            Leaf healthy surface
        """
        leaf_elements = self.leaf_elements
        return sum(g.node(vid).healthy_surface for vid in leaf_elements)