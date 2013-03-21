""" Gather existing models of dispersal of fungus propagules.

"""

# Imports #########################################################################
import random

# Random dispersal ################################################################

class RandomDispersal:
    """ Template class for a dispersal model that complies with the guidelines of Alep.
    
    A class for a model of dispersal must include a method 'disperse'. In this example,
    dispersal units are randomly distributed.
    
    """
    def disperse(self, scene, dispersal_units):
        """ Example method for dispersal with random distribution.
        
        Parameters
        ----------
        scene :
            Scene containing the simulated system
        dispersal_units : dict
            Dispersal units emitted by the lesions on leaves
            
        Returns
        -------
        deposits : dict
            Dispersal units deposited on new position on leaves
        """
        vids = [geom.id for geom in scene]
        n = len(vids)
        deposits = {}

        for vid, dlist in dispersal_units.iteritems():
            for d in dlist:
                if random.random() < 0.5:
                    idx = random.randint(0,n-1)
                    v = vids[idx]
                    deposits.setdefault(v,[]).append(d)
         
        return deposits