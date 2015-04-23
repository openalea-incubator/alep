""" Gather different strategies for modeling inoculation of fungus dispersal units on a MTG.

"""

# Imports #########################################################################
import random
import numpy as np
import collections
from alinea.alep.fungal_objects import DispersalUnit, Lesion
from alinea.alep.architecture import get_leaves
from openalea.plantgl import all as pgl

# Random inoculation ##############################################################

class RandomInoculation:
    """ Template class for inoculum allocation that complies with the guidelines of Alep.
    
    A class for a model of dispersal must include a method 'allocate'. In this example,
    dispersal units are randomly distributed.
    
    """
    
    def allocate(self, g, inoculum, label='LeafElement'):
        """ Select randomly elements of the MTG and allocate them a random part of the inoculum.

        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        inoculum: list of dispersal units OR list of lesions
            Source of fungal objects to distribute on the MTG
        label: str
            Label of the part of the MTG concerned by the calculation
            
        Returns
        -------
        None
            Update directly the MTG
        """        
        vids = [n for n in g if g.label(n).startswith(label)]
        areas = g.property('area')
        vids = [vid for vid in vids if vid in g.property('geometry') if areas[vid]>0.]
        n = len(vids)
        
        if n>0:
            for i in inoculum:
                idx = random.randint(0,n-1)
                v = vids[idx]
                leaf = g.node(v)
                # Set a position for i :
                if i.position==None:
                    i.position = [random.random(), 0] # TODO : improve
                
                #  Attach it to the leaf
                if isinstance(i, Lesion):
                    try:
                        leaf.lesions.append(i)
                    except:
                        leaf.lesions = [i]            
                elif isinstance(i, DispersalUnit):
                    i.deposited()
                    try:
                        leaf.dispersal_units.append(i)
                    except:
                        leaf.dispersal_units = [i]

class InoculationYoungLeaves:
    """ Template class for inoculum allocation that complies with the guidelines of Alep.
    
    A class for a model of dispersal must include a method 'allocate'. In this example,
    dispersal units are distributed only on young leaves.
    
    """
    def __init__(self, age_max=5.):
        """ Initialize the model of inoculation.
        
        Parameters
        ----------
        age_max: float
            Maximal age for a leaf to be selected for infection (in days)
        """
        self.age_max=age_max
    
    def allocate(self, g, inoculum, label='LeafElement'):
        """ Select randomly elements of the MTG and allocate them a random part of the inoculum.

        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        inoculum: list of dispersal units OR list of lesions
            Source of fungal objects to distribute on the MTG
        label: str
            Label of the part of the MTG concerned by the calculation
            
        Returns
        -------
        None
            Update directly the MTG
        """        
        vids = [n for n in g if g.label(n).startswith(label) and g.node(n).age<self.age_max]
        areas = g.property('area')
        vids = [vid for vid in vids if vid in g.property('geometry') if areas[vid]>0.]
        n = len(vids)
        
        if n>0:
            for i in inoculum:
                idx = random.randint(0,n-1)
                v = vids[idx]
                leaf = g.node(v)
                # Set a position for i :
                if i.position==None:
                    i.position = [random.random(), 0] # TODO : improve
                
                #  Attach it to the leaf
                if isinstance(i, Lesion):
                    try:
                        leaf.lesions.append(i)
                    except:
                        leaf.lesions = [i]            
                elif isinstance(i, DispersalUnit):
                    i.deposited()
                    try:
                        leaf.dispersal_units.append(i)
                    except:
                        leaf.dispersal_units = [i]
                        
class InoculationLowerLeaves(object):
    """ Template class for inoculum allocation that complies with the guidelines of Alep.
    
    A class for a model of dispersal must include a method 'allocate'. In this example,
    dispersal units are randomly distributed.
    
    """
        
    def __init__(self, max_height=20.):
        self.max_height = max_height
    
    def is_iterable(self, obj):
        """ Test if object is iterable """
        import collections
        return isinstance(obj, collections.Iterable)
    
    def get_leaf_height(self, leaf_geom):
        from openalea.plantgl import all as pgl
        tesselator = pgl.Tesselator()
        bbc = pgl.BBoxComputer(tesselator)
        if self.is_iterable(leaf_geom):
            bbc.process(pgl.Scene(leaf_geom))
        else:
            bbc.process(pgl.Scene([pgl.Shape(leaf_geom)]))
        return bbc.result.getCenter()[2]
        
    def allocate(self, g, inoculum, label='LeafElement'):
        """ Select lower elements of the MTG and allocate them a random part of the inoculum.

        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        inoculum: list of dispersal units OR list of lesions
            Source of fungal objects to distribute on the MTG
        label: str
            Label of the part of the MTG concerned by the calculation
            
        Returns
        -------
        None
            Update directly the MTG
        """
        from alinea.alep.architecture import get_leaves
        
        leaves = get_leaves(g, label='LeafElement')
        geometries = g.property('geometry')
        vids = list(leaf for leaf in leaves if leaf in geometries.iterkeys()
                        and self.get_leaf_height(geometries[leaf])<=self.max_height)
        n = len(vids)
        
        if n>0:
            for i in inoculum:
                idx = random.randint(0,n-1)
                v = vids[idx]
                leaf = g.node(v)
                # Set a position for i :
                i.position = [0, 0] # TODO : improve
                
                #  Attach it to the leaf
                if isinstance(i, Lesion):
                    try:
                        leaf.lesions.append(i)
                    except:
                        leaf.lesions = [i]            
                elif isinstance(i, DispersalUnit):
                    #i.deposited()
                    try:
                        leaf.dispersal_units.append(i)
                    except:
                        leaf.dispersal_units = [i]

class AirborneInoculum:
    """ Represent an airborne innoculum """
    def __init__(self, fungus, mutable = False, 
                 density = 100, domain_area = 1.):
        self.fungus = fungus
        self.mutable = mutable
        self.density = density
        self.domain_area = domain_area

    def emission(self, g):
        return [self.fungus.dispersal_unit(mutable=self.mutable)
                 for i in range(self.domain_area / self.density)]

def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)

class AirborneContamination:
    """ Model of airborne inoculation """
    def __init__(self, fungus, mutable = False, 
                 domain_area = 1., convUnit = 100.,
                 layer_thickness = 1., k_beer = 0.5):
        self.fungus = fungus
        self.mutable = mutable        
        self.domain_area = domain_area 
        self.convUnit = convUnit
        self.layer_thickness = layer_thickness
        self.k_beer = k_beer

    def leaves_in_grid(self, g, label = 'LeafElement'):
        geometries = g.property('geometry')
        centroids = g.property('centroid')
        tesselator = pgl.Tesselator()
        bbc = pgl.BBoxComputer(tesselator)      
        leaves = get_leaves(g, label=label)
        leaves = [l for l in leaves if l in geometries]
        
        # Get centroids        
        def centroid(vid):
            if is_iterable(geometries[vid]):
                bbc.process(pgl.Scene(geometries[vid]))
            else:
                bbc.process(pgl.Scene([pgl.Shape(geometries[vid])]))
            center = bbc.result.getCenter()
            centroids[vid] = center
        
        for vid in leaves:
            centroid(vid)
            
        # Define grid (horizontal layers)
        zs = [c[2] for c in centroids.itervalues()]
        minz = min(zs)
        maxz = max(zs) + self.layer_thickness
        layers = {l:[] for l in np.arange(minz, maxz, self.layer_thickness)}
        
        # Distribute leaves in layers
        for vid, coords in centroids.iteritems():
            z = coords[2]
            ls = layers.keys()
            i_layer = np.where(map(lambda x: x<=z<x+self.layer_thickness 
                                    if z!=maxz - self.layer_thickness
                                    else x<=z<=x+self.layer_thickness , ls))[0]
            if len(i_layer) > 0.:
                layers[ls[i_layer]].append(vid)
        
        self.layers = layers
        
    def contaminate(self, g, 
                    density_dispersal_units = 100., 
                    weather_data = None,
                    label='LeafElement'):
        areas = g.property('area')
        self.leaves_in_grid(g, label=label)
        deposits = {}
        nb_dus = density_dispersal_units * self.domain_area * self.convUnit
        sorted_layers = sorted(self.layers.keys(), reverse = True)
        for layer in sorted_layers:
            if nb_dus > 0.:
                vids = self.layers[layer]
                area_layer = sum([areas[vid] for vid in vids])
    #            beer_factor = 1-np.exp(-self.k_beer*area_layer/(self.domain_area*self.convUnit))
                proba_du_layer = 1-np.exp(-self.k_beer*(area_layer/self.convUnit**2)/self.domain_area)
    #            proba_du_layer = beer_factor / nb_dus
                distribution_by_leaf = np.random.binomial(nb_dus, 
                                                          proba_du_layer,
                                                          len(vids))
#                import pdb
#                pdb.set_trace()                
                deposits.update({lf:distribution_by_leaf[i_lf]
                                 for i_lf, lf in enumerate(vids)})
                nb_dus -= sum(distribution_by_leaf)
        return deposits
            
    def plot_distri_layers(self, g, density_dispersal_units = 100.):
        from alinea.alep.architecture import set_property_on_each_id
        from alinea.alep.alep_color import alep_colormap, green_yellow_red
        from alinea.alep.disease_outputs import plot3d_transparency
        from openalea.plantgl.all import Viewer
        # Compute severity by leaf
        deposits = self.contaminate(g)
        set_property_on_each_id(g, 'nb_dispersal_units', deposits)
    
        # Visualization
        g = alep_colormap(g, 'nb_dispersal_units', cmap=green_yellow_red(), 
                          lognorm=False, zero_to_one=False)

        geometries = g.property('geometry') 
        leaves = get_leaves(g, label='LeafElement')
        leaves = [l for l in leaves if l in geometries] 
        transp = {vid:0. for k,v in self.layers.iteritems() for vid in v}
        set_property_on_each_id(g, 'transparency', transp)                         
        scene = plot3d_transparency(g)
        Viewer.display(scene)