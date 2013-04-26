""" Gather different strategies for modeling inoculation of fungus dispersal units on a MTG.

"""

# Imports #########################################################################
import random

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
        inoculum: list of DUs
            Source of dispersal units to distribute on the MTG
        label: str
            Label of the part of the MTG concerned by the calculation
            
        Returns
        -------
        None
            Update directly the MTG
        """        
        vids = [n for n in g if g.label(n).startswith(label)]
        n = len(vids)
        for vid in vids:
            g.node(vid).dispersal_units = []
            
        for i in inoculum:
            idx = random.randint(0,n-1)
            v = vids[idx]
            leaf = g.node(v)
            if not 'dispersal_units' in leaf.properties():
                leaf.dispersal_units = []
            # Deposit a DU from inoculum on node v of the MTG
            i.deposited()
            i.position = [0, 0] # TODO : improve
            leaf.dispersal_units.append(i)