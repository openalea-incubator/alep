""" Gather different strategies for checking position of dispersal units on a MTG
    that could prevent them from infecting the leaves.

"""

# Imports #########################################################################
import random
from math import sqrt, pi
from alinea.alep.disease_outputs import compute_severity_by_leaf

# Random inoculation ##############################################################
class BiotrophDUProbaModel:
    """ Template class for checking position of dispersal units on a MTG that 
        complies with the guidelines of Alep.
    
    A class for a model of dispersal unit position checking must contain the 
    following method:
        - 'control_position': work at leaf scale.
        Affect the property 'can_infect_at_position' of the DUs.
        
    In this example, a dispersal unit can only infect an healthy tissue. 
    Its property 'position' in its interface does not need to be given.
    The probability for it to be on a healthy tissue is computed with 
    the ratio between leaf healthy and total areas.
    """   
    def control_position(self, g, label='LeafElement'):
        """ Control if the dispersal units can infect at their current position.
        
        Call the method 'can_not_infect_at_position' of DU interface eventually
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
            'dispersal_units' are stored in the MTG as a property
        label: str
            Label of the part of the MTG concerned by the calculation
            
        Returns
        -------
        g: MTG
            Updated MTG representing the canopy
        """
        from alinea.alep.disease_outputs import compute_severity_by_leaf
        dispersal_units = g.property('dispersal_units')
        severities = compute_severity_by_leaf(g, label=label)
        for vid, du in dispersal_units.iteritems():
            # By leaf element, keep only those which are deposited and active
            du = [d for d in du if d.is_active]
            ratio = severities[vid]/100.
            if du:
                nb_not_on_green = int(ratio*len(du))
                for dispersal_unit in random.sample(du, nb_not_on_green):
                    # Is not on green tissue
                    dispersal_unit.can_not_infect_at_position()

class BiotrophDUPositionModel:
    """ Template class for checking position of dispersal units on a MTG that 
        complies with the guidelines of Alep.
    
    A class for a model of dispersal unit position checking must contain the 
    following method:
        - 'control_position': work at leaf scale.
        Affect the property 'can_infect_at_position' of the DUs.
        
    In this example, a dispersal unit can only infect an healthy tissue. 
    Its property 'position' in its interface is required.
    """   
    def control_position(self, g, label='LeafElement'):
        """ Control if the lesion can infect at its current position.
        
        Call the method 'can_not_infect_at_position' of DU interface eventually
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
            'dispersal_units' are stored in the MTG as a property
            
        Returns
        -------
        g: MTG
            Updated MTG representing the canopy
        """
        dispersal_units = g.property('dispersal_units')
        severities = compute_severity_by_leaf(g)
        for vid, du in dispersal_units.iteritems():
            # By leaf element, keep only those which are deposited and active
            # du = [d for d in du if d.is_active and d.status=="deposited"]
            du = [d for d in du if d.is_active]
            leaf = g.node(vid)
        
            # Test for senescence:
            for DU in du:
                if DU.position[0] >= leaf.position_senescence:
                    DU.can_not_infect_at_position()
                elif random.random()<(severities[vid]/100.):
                    DU.can_not_infect_at_position()
                # # Test for spaces taken by other lesions
                # try:
                    # taken_x_coord = [l.position[0] for l in leaf.lesions]
                    # # The hypothesis is made that lesions are circular
                    # radius = [sqrt(l.surface/pi) for l in leaf.lesions]
                    
                    # for x in taken_x_coord:
                        # if DU.can_infect_at_position() and x-radius<=DU.position[0]<=x+radius:
                            # DU.can_not_infect_at_position()
                            # break                  
                # except:
                    # pass
