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

class InoculationFirstLeaves:
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
        vids = [12, 13, 23, 24, 34, 35, 178, 179, 311, 312, 444, 445]
        n = len(vids)
            
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