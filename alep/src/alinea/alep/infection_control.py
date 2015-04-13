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
    def control(self, g, label='LeafElement'):
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
                    dispersal_unit.disable()

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
    def control(self, g, label='LeafElement'):
        """ Control if the lesion can infect at its current position. If not, disable the dispersal unit.
        
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
        labels = g.property('label')
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        dispersal_units = {k:v for k,v in g.property('dispersal_units').iteritems() if len(v)>0.}
        areas = g.property('area')
        lesions = g.property('lesions')
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith(label)]
            leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], []) 
            les_surf = sum([les.surface for les in leaf_lesions])
            leaf_area = sum([areas[lf] for lf in leaf])
            ratio_les_surface = min(1, round(les_surf,3)/round(leaf_area,3)) if round(leaf_area,3)>0. else 0.
            
            for vid in set(leaf) & set(dispersal_units):
                dus_to_keep = []
                dus = []
                for du in dispersal_units[vid]:
                    if du.is_active:
                        if du.status == 'deposited':
                            dus_to_keep.append(du)
                        elif du.status == 'emitted':
                            dus.append(du)
                            
                if len(dus)>0 and dus[0].fungus.group_dus == True:
                    total_nb_dus = len(sum([du.position for du in dus],[]))
                else:
                    total_nb_dus = len(dus)
                nb_on_lesions = int(total_nb_dus*ratio_les_surface)
                for du in range(nb_on_lesions):
                    random.shuffle(dus)
                    dus[0].position = dus[0].position[1:]
                    if dus[0].nb_dispersal_units==0.:
                        dus[0].disable()
                        dus = dus[1:]
                for du in dus:
                    du.set_status(status = 'deposited')
                dispersal_units[vid] = dus_to_keep + dus