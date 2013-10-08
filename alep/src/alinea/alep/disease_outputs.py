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
    return sum(len(du) for du in dispersal_units.itervalues())
    
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

def compute_lesion_areas_by_leaf(g, label='LeafElement'):
    """ Compute lesion area on each part of the MTG given by the label.
    
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
    return {vid:(sum(l.surface for l in lesions[vid])
            if vid in lesions.keys() else 0.) for vid in vids} 

def compute_green_lesion_areas_by_leaf(g, label='LeafElement'):
    """ Compute lesion areas on each green part of the MTG given by the label.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    green_lesion_area_by_leaf: dict([id:lesion_area])
        Surface of the lesions on each green part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    lesions = g.property('lesions')
    return {vid:(sum(l.surface for l in lesions[vid] if not l.is_senescent)
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
    healthy_by_leaf: dict([id:healthy_area])
        Healthy area on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    green_areas = g.property('green_area')
    areas = g.property('area')
    lesion_areas = compute_lesion_areas_by_leaf(g, label)
    green_lesion_areas = compute_green_lesion_areas_by_leaf(g, label)
    senescent_lesions = {vid:(lesion_areas[vid]-green_lesion_areas[vid]) for vid in vids}
    green_areas = {vid:min(green_areas[vid], areas[vid]-senescent_lesions[vid]) for vid in vids}
    return {vid:(green_areas[vid] - green_lesion_areas[vid] 
            if round(green_areas[vid], 10)>round(green_lesion_areas[vid], 10) else 0.)
            for vid in vids}
    
def compute_severity_by_leaf(g, label='LeafElement'):
    """ Compute severity of the disease on each part of the MTG given by the label.
    
    Severity is the ratio between disease surface and total leaf area (in %).
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    severity_by_leaf: dict([id:severity])
        Severity on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    total_areas = g.property('area')
    lesion_areas = compute_lesion_areas_by_leaf(g, label)
    return {vid:(100*lesion_areas[vid]/float(total_areas[vid]) if total_areas[vid]>0. else 0.) for vid in vids}

def compute_necrosis_percentage_by_leaf(g, label='LeafElement'):
    """ Compute necrosis percentage on each part of the MTG given by the label.
    
    Necrosis percentage is the ratio between necrotic area and total leaf area.
    A tissue is necrotic if it is covered by a lesion in one of these states:
        - NECROTIC
        - SPORULATING
        - EMPTY
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrosis_by_leaf: dict([id:necrosis_percentage])
        Necrosi percentage on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    total_areas = g.property('area')
    lesions = g.property('lesions')
    necrotic_areas = {}
    for vid in total_areas.iterkeys():
        if vid in lesions.keys():
            necrotic_areas[vid] = sum(lesion.necrotic_area for lesion in lesions[vid])
        else:
            necrotic_areas[vid] = 0.
    return {vid:(100*necrotic_areas[vid]/float(total_areas[vid]) if total_areas[vid]>0. else 0.) for vid in vids}

def compute_necrotic_area_by_leaf(g, label='LeafElement'):
    """ Compute necrosis percentage on each part of the MTG given by the label.
    
    Necrosis percentage is the ratio between necrotic area and total leaf area.
    A tissue is necrotic if it is covered by a lesion in one of these states:
        - NECROTIC
        - SPORULATING
        - EMPTY
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrotic_area_by_leaf: dict([id:necrosis_percentage])
        Necrotic area on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    total_areas = g.property('area')
    lesions = g.property('lesions')
    necrotic_areas = {}
    for vid in total_areas.iterkeys():
        if vid in lesions.keys():
            necrotic_areas[vid] = sum(lesion.necrotic_area for lesion in lesions[vid])
        else:
            necrotic_areas[vid] = 0.
    return necrotic_areas
    
def compute_total_severity(g, label='LeafElement'):
    """ Compute disease severity on the whole plant.
    
    Severity is the ratio between disease surface and green leaf area (in %).
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    severity: float
        Ratio between disease surface and green leaf area (in %)
    """
    from numpy import mean
    severities = compute_severity_by_leaf(g, label=label)
    return mean(severities.values())
    
def compute_total_necrosis_percentage(g, label='LeafElement'):
    """ Compute necrosis percentage on the whole plant.
    
    Necrosis percentage ratio between necrotic (and sporulating) disease surface and total area of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrosis_percentage: float
        Ratio between necrotic (and sporulating) disease area and total area of leaves (in %)
    """
    # Leaf
    # green_areas = g.property('green_area')
    # total_green_area = sum(green_areas.values())
    
    # Disease
    # lesions = g.property('lesions')
    # if lesions:
        # lesions = [l for les in lesions.values() for l in les 
                    # if (l.status>=l.fungus.NECROTIC and not l.is_senescent)]
        # total_necrotic_area = sum(l.surface for l in lesions)
    # else:
        # total_necrotic_area = 0.

    # Compute ratio, i.e. necrosis percentage
    # return 100 * total_necrotic_area / total_green_area if total_green_area > 0. else 0.
    
    from numpy import mean
    nec = compute_necrosis_percentage_by_leaf(g, label=label)
    return mean(nec.values())

def compute_total_necrotic_area(g, label='LeafElement'):
    """ Compute necrosis percentage on the whole plant.
    
    Necrosis percentage ratio between necrotic (and sporulating) disease surface and total area of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrotic_area: float
        Total area of leaves covered by necrotic surfaces of lesions (in cm2)
    """
    from numpy import mean
    nec = compute_necrotic_area_by_leaf(g, label=label)
    return sum(nec.values())
    
######################################################################
from numpy import mean

class LeafInspector:
    def __init__(self, g, blade_id=None, label='LeafElement'):
        """ Find the ids of the leaf elements on the chosen blade.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        blade_id: int
            Index of the blade to be inspected
        label: str
            Label of the part of the MTG concerned by the calculation
        """
        from numpy import mean
        labels = g.property('label')
        self.label = label
        # Find leaf elements on the blade
        self.ids = [id for id in g.components(blade_id) if labels[id].startswith(label)]
        # Initialize leaf properties to save
        self.leaf_area = []
        self.leaf_green_area = []  
        self.leaf_healthy_area = []
        self.leaf_disease_area = []
        # Initialize variables relative to DUs
        self.nb_dus = []
        self.nb_dus_on_green = []
        self.nb_dus_on_healthy = []
        # Initialize variables relative to number of infections
        self.previous_nb_lesions = 0.
        self.nb_infections = []
        # Initialize surfaces in state
        self.surface_inc = []
        self.surface_chlo = []
        self.surface_nec = []
        self.surface_spo = []
        self.surface_total_nec = []
        # Initialize ratios (surfaces in state compared to leaf area)
        self.ratio_inc = []
        self.ratio_chlo = []
        self.ratio_nec = []
        self.ratio_spo = []
        self.ratio_total_nec = []
        # Initialize total severity
        self.severity = []
        # Initialize necrosis percentage
        self.necrosis = []
    
        # Temporary
        self.previous_nb_lesions = 0.
    
    def update_du_variables(self, g):
        """ Save counts of dispersal units.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy 
        """
        dus = 0.
        dus_on_green = 0.
        dus_on_healthy = 0.
        for id in self.ids:
            leaf = g.node(id)
            try:
                dus += len(leaf.dispersal_units)
                dus_on_green = len([du for du in leaf.dispersal_units if du.position[0] < leaf.position_senescence])
                dus_on_healthy = len([du for du in leaf.dispersal_units if du.can_infect_at_position])
            except:
                pass
        self.nb_dus.append(dus)
        self.nb_dus_on_green.append(dus_on_green)
        self.nb_dus_on_healthy.append(dus_on_healthy)
        
    def update_variables(self, g):
        """ Save leaf properties and disease properties.
        
        Save leaf properties and update the computation of severity, 
        necrosis percentage and ratios.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy   
        """
        area = 0.
        green_area = 0.
        disease_area = 0.
        lesion_list = []
        for id in self.ids:
            leaf = g.node(id)
            area += leaf.area
            green_area += leaf.green_area
            try:
                lesion_list += leaf.lesions
            except:
                pass
        
        if len(lesion_list)>0:
            green_lesion_area = sum(l.surface for l in lesion_list if not l.is_senescent)
            disease_area = sum(l.surface for l in lesion_list)
        else:
            green_lesion_area = 0.
            disease_area = 0.
            
        self.leaf_area.append(area)
        self.leaf_green_area.append(green_area)
        self.leaf_healthy_area.append(green_area - green_lesion_area)
        self.leaf_disease_area.append(disease_area)
        self.compute_ratios(g)
        self.compute_severity(g)
        self.compute_necrosis(g)
    
    def update_area(self,g):
        areas = g.property('area')
        self.leaf_area.append(sum([areas[id] for id in self.ids]))
    
    def update_green_area(self, g):
        green_areas = g.property('green_area')
        self.leaf_green_area.append(sum([green_areas[id] for id in self.ids]))
        
    def update_healthy_area(self, g):
        healthy_areas = g.property('healthy_area')
        self.leaf_healthy_area.append(sum([healthy_areas[id] for id in self.ids]))
    
    def update_disease_area(self, g):
        disease_area = compute_lesion_areas_by_leaf(g, label=self.label)
        self.leaf_disease_area.append(sum([disease_area[id] for id in self.ids]))
        
        lesions = g.property('lesions')
        nb_lesions = 0.
        for id in self.ids:
            if id in lesions:
                nb_lesions += len(lesions[id])
        if nb_lesions < self.previous_nb_lesions:
            import pdb
            pdb.set_trace()
        self.previous_nb_lesions = nb_lesions
    
    def compute_nb_infections(self, g):
        """ Compute the number of infections during time step.
        
        Number of infections corresponds to number of new lesions created.        
        """
        lesions = g.property('lesions')
        nb_infections = 0.
        if len(lesions)>0:
            lesions = sum([lesions[id] for id in self.ids if id in lesions.keys()], [])
            if lesions!=0:
                nb_lesions = len(lesions)
                nb_infections = max(0., nb_lesions - self.previous_nb_lesions)
                self.previous_nb_lesions = nb_lesions

        self.nb_infections.append(nb_infections)
            
    def compute_severity(self, g):
        """ Compute severity on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        severities = compute_severity_by_leaf(g, label=self.label)
        self.severity.append(mean([severities[id] for id in self.ids]))
    
    def compute_necrosis(self, g):
        """ Compute necrosis percentage on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        nec = compute_necrosis_percentage_by_leaf(g, label=self.label)
        self.necrosis.append(mean([nec[id] for id in self.ids]))
    
    def compute_ratios(self, g):
        """ Compute surface of lesions in chosen state on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        total_area = sum([g.node(id).area for id in self.ids])
        lesion_list = []
        for id in self.ids:
            leaf = g.node(id)
            try:
                lesion_list += leaf.lesions
            except:
                pass
        
        surface_inc = 0.
        surface_chlo = 0.
        surface_nec = 0.
        surface_spo = 0.       
        if len(lesion_list)>0:
            for l in lesion_list:
                l.compute_all_surfaces()
                surface_inc += l.surface_inc
                surface_chlo += l.surface_chlo
                surface_nec += l.surface_nec
                surface_spo += l.surface_spo

        self.surface_inc.append(surface_inc)
        self.surface_chlo.append(surface_chlo)
        self.surface_nec.append(surface_nec)
        self.surface_spo.append(surface_spo)
        self.surface_total_nec.append(surface_nec+surface_spo)
        self.ratio_inc.append(100. * surface_inc / total_area if total_area>0. else 0.)
        self.ratio_chlo.append(100. * surface_chlo / total_area if total_area>0. else 0.)
        self.ratio_nec.append(100. * surface_nec / total_area if total_area>0. else 0.)
        self.ratio_spo.append(100. * surface_spo / total_area if total_area>0. else 0.)
        self.ratio_total_nec.append(100. * (surface_nec+surface_spo) / total_area if total_area>0. else 0.)
    
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