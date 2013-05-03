""" Gather different strategies for checking position of dispersal units on a MTG
    that could prevent them from infecting the leaves.

"""

# Imports #########################################################################
from random import random
from math import sqrt, pi

# Random inoculation ##############################################################
class BiotrophDUProbaModel:
    """ Template class for checking position of dispersal units on a MTG that 
        complies with the guidelines of Alep.
    
    A class for a model of dispersal unit position checking must contain the 
    following method:
        - 'check_position': work at the scale of a single dispersal unit.
        Disable it if it can not infect the tissue.
        
    In this example, a dispersal unit can only infect an healthy tissue. 
    Its property 'position' in its interface does not need to be given.
    The probability for it to be on a healthy tissue is computed with 
    the ratio between leaf healthy and total surfaces.
    """   
    def check_position(self, DU=None, leaf=None):
        """ Disable the dispersal unit if its position prevents it from 
            infecting the tissue.
        
        Parameters
        ----------
        DU: dispersal unit instantiation
            The DU with properties (e.g. status, activity, number of spores, position)
        leaf: Leaf sector node of an MTG 
            A leaf sector with the properties 'healty_surface' and 'surface'
            
        Returns
        -------
        Turn off the activity of the DU if not on healthy tissue
        """
        healthy_surface = leaf.healthy_surface
        surface = leaf.surface
        ratio = healthy_surface / surface if surface>0. else 0.

        if random() > ratio:
            # Then the dispersal unit can not infect
            DU.disable()
            
class BiotrophDUPositionModel:
    """ Template class for checking position of dispersal units on a MTG that 
        complies with the guidelines of Alep.
    
    A class for a model of dispersal unit position checking must contain the 
    following method:
        - 'check_position': work at the scale of a single dispersal unit.
        Disable it if it can not infect the tissue.
        
    In this example, a dispersal unit can only infect an healthy tissue. 
    Its property 'position' in its interface is required.
    """   
    def check_position(self, DU=None, leaf=None):
        """ Disable the dispersal unit if its position prevents it from 
            infecting the tissue.
        
        Parameters
        ----------
        DU: dispersal unit instantiation
            The DU with properties (e.g. status, activity, number of spores, position)
        leaf: Leaf sector node of an MTG 
            A leaf sector with the properties 'healty_surface' and 'surface'
            
        Returns
        -------
        Turn off the activity of the DU if not on healthy tissue
        """
        # Test for senescence:
        if leaf.position_senescence and DU.position[0] >= leaf.position_senescence:
            DU.disable()
            return
            
        # Test for spaces taken by other lesions
        try:
            taken_x_coord = [l.position[0] for l in leaf.lesions]
            taken_y_coord = [l.position[1] for l in leaf.lesions]
            # The hypothesis is made that lesions are circular
            radius = [sqrt(l.surface/pi) for l in leaf.lesions]
        except:
            taken_x_coord = []        
            taken_y_coord = []        
        
        if taken_x_coord:
            for x in taken_x_coord:
                if x-radius<=DU.position[0]<=x+radius:
                    DU.disable()
                    return
            for y in taken_y_coord:
                if y-radius<=DU.position[1]<=y+radius:
                    DU.disable()
                    return
        