""" Gather different strategies for modeling inoculation of fungus dispersal units on a MTG.

"""

# Imports #########################################################################
import random
from alinea.alep.fungal_objects import DispersalUnit, Lesion

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
                    i.deposited()
                    try:
                        leaf.dispersal_units.append(i)
                    except:
                        leaf.dispersal_units = [i]