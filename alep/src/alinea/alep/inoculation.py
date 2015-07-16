""" Gather different strategies for modeling inoculation of fungus dispersal units on a MTG.

"""

# Imports #########################################################################
import random
import numpy as np
import collections
from alinea.alep.fungal_objects import DispersalUnit, Lesion, Fungus
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

def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)

class AirborneContamination:
    """ Model of airborne inoculation """
    def __init__(self, fungus, group_dus = False, mutable = False, 
                 domain_area = 1., convUnit = 0.01,
                 layer_thickness = 1., k_beer = 0.5):
        if fungus is not None:
            self.fungus = fungus
        else:
            self.fungus = Fungus()
        self.group_dus = group_dus
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

    def emission(self, g, weather_data = None, density_dispersal_units = 0.):
        return density_dispersal_units * self.domain_area
    
    def contaminate(self, g, nb_dus = 0., weather_data = None,
                    label='LeafElement'):
                        
        def sum_nb(nb_leaves, nb_du):
            if nb_leaves == 1:
                return [nb_du]
            elif nb_du == 0:
                return [0] + sum_nb(nb_leaves-1, nb_du)
            else:
                nb_du_avg = float(nb_du/nb_leaves)
                nb_du_sup = 2.*nb_du_avg
                if nb_du_sup >= 1:
                    nb_on_vid = int(round(max(0, min(nb_du, np.random.normal(nb_du_avg, nb_du_sup)))))
                else:
                    nb_on_vid = 1 if np.random.random()<nb_du_sup else 0
                return [nb_on_vid] + sum_nb(nb_leaves-1, nb_du - nb_on_vid)
        
        areas = g.property('area')
        self.leaves_in_grid(g, label=label)
        deposits = {}
        sorted_layers = sorted(self.layers.keys(), reverse = True)
        for layer in sorted_layers:
            if nb_dus > 0.:              
                vids = self.layers[layer]
                area_layer = sum([areas[vid] for vid in vids])
                proba_du_layer = 1-np.exp(-self.k_beer*(area_layer*self.convUnit**2)/self.domain_area)
                nb_dus_in_layer = np.random.binomial(nb_dus, proba_du_layer)
                if len(vids)>0.:                    
                    distribution_by_leaf = sum_nb(len(vids), nb_dus_in_layer)
                    np.random.shuffle(distribution_by_leaf)
                    deposits.update({lf:distribution_by_leaf[i_lf]
                                     for i_lf, lf in enumerate(vids)
                                     if distribution_by_leaf[i_lf] > 0.})
                    nb_dus -= sum(distribution_by_leaf)
        
        for vid, nb_dus in deposits.iteritems():
            if self.fungus.group_dus==True:
                du = self.fungus.dispersal_unit()
                du.set_nb_dispersal_units(nb_dispersal_units = nb_dus)
                deposits[vid] = [du]
            else:
                dus = []
                for d in nb_dus:
                    du = self.fungus.dispersal_unit()
                    dus.append(du)
                deposits[vid] = dus
        return deposits
            
    def view_distri_layers(self, g, density_dispersal_units = 1000., 
                           vmax = None):
        from alinea.alep.architecture import set_property_on_each_id
        from alinea.alep.alep_color import alep_colormap, green_yellow_red
        from alinea.alep.disease_outputs import plot3d_transparency
        from openalea.plantgl.all import Viewer
        import matplotlib.pyplot as plt
        # Compute severity by leaf
        deposits = {k:sum([du.nb_dispersal_units for du in v]) for k,v in 
                    self.contaminate(g, density_dispersal_units).iteritems()}
        set_property_on_each_id(g, 'nb_dispersal_units', deposits)
    
        # Visualization
        if vmax is None:
            vmax = 2*max(deposits.values())/3
        g = alep_colormap(g, 'nb_dispersal_units', cmap=green_yellow_red(), 
                          lognorm=False, zero_to_one=False)

        d = [np.arange(vmax)]
        fig, ax = plt.subplots()
        ax.imshow(d, cmap = green_yellow_red())

        geometries = g.property('geometry') 
        leaves = get_leaves(g, label='LeafElement')
        leaves = [l for l in leaves if l in geometries] 
        transp = {vid:0. for k,v in self.layers.iteritems() for vid in v}
        set_property_on_each_id(g, 'transparency', transp)   
        for id in g:
            if not id in deposits:
                g.node(id).color = (1,1,1)
                g.node(id).transparency = 0.7
            elif deposits[id]==0.:
                g.node(id).color = (1,1,1)
                g.node(id).transparency = 0.7
            else:
                g.node(id).transparency = 0.                         
        scene = plot3d_transparency(g)
        Viewer.display(scene)
        
    def plot_distri_layers(self, g, density_dispersal_units = 1000.):
        import pandas as pd
        import matplotlib.pyplot as plt
        self.leaves_in_grid(g)
        deposits = {k:sum([du.nb_dispersal_units for du in v]) for k,v in 
                    self.contaminate(g, density_dispersal_units).iteritems()}
        depo_layer = {k:sum([deposits[vid] for vid in v if vid in deposits])
                        for k,v in self.layers.iteritems()}
        df = pd.DataFrame([[k,v] for k, v in depo_layer.iteritems()])
        df = df.sort(0)
        df[2] = df[1]/df[1].sum()
        fig, ax = plt.subplots()
        ax.plot(df[2], df[0], 'k')
        ax.set_ylabel('Height', fontsize = 16)
        ax.set_xlabel('Proportion of deposits', fontsize = 16)
        return df