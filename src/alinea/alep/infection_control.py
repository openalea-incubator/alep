""" Gather different strategies for checking position of dispersal units on a MTG
    that could prevent them from infecting the leaves.

"""

# Imports #########################################################################
import numpy as np

import collections        
def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)

# Random inoculation ##############################################################
class BiotrophDUProbaModel:
    """ Template class for checking position of dispersal units on a MTG that 
        complies with the guidelines of Alep.
    
    A class for a model of dispersal unit position checking must contain the 
    following method:
        - 'control': work at leaf scale.
        Affect the property 'can_infect_at_position' of the DUs.
        
    In this example, a dispersal unit can only infect an healthy tissue. 
    Its property 'position' in its interface does not need to be given.
    The probability for it to be on a healthy tissue is computed with 
    the ratio between leaf healthy and total areas.
    """
    def __init__(self, max_lesion_density=150, 
                 age_infection=False, fungus='septoria'):
        self.max_lesion_density = max_lesion_density
        self.age_infection=age_infection
        self.fungus=fungus

    def control(self, g, label='LeafElement'):
        """ Control if the dispersal units can infect at their current position.
        
        Call the method 'can_not_infect_at_position' of DU interface eventually
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
            'dispersal_units' are stored in the MTG as a property
            'lesions' are stored in the MTG as a property
            'green_area' is the green area of leaf elements
            'area' is the total area of leaf elements
        label: str
            Label of the part of the MTG concerned by the calculation
            
        Returns
        -------
        g: MTG
            Updated MTG representing the canopy
        """
        labels = g.property('label')
        bids = (v for v,l in labels.items() if l.startswith('blade'))
        DUs = g.property('dispersal_units')
        dispersal_units = {k:v for k,v in DUs.items() if len(v)>0.}
        areas = g.property('area')
        green_areas = g.property('green_area')
        lesions = g.property('lesions')
        
        # Temp
        ages = g.property('age')
#        from alinea.adel.newmtg import adel_labels
#        a_labs = adel_labels(g)
        if not is_iterable(self.fungus):
            fungi = [self.fungus]
        else: 
            fungi = self.fungus        
        
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith(label)]
            leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], [])
            les_surf = sum([les.surface for les in leaf_lesions])
            leaf_area = sum([areas[lf] for lf in leaf])
            leaf_green_area = sum([green_areas[lf] for lf in leaf])
            ratio_les_surface = min(1, round(les_surf,3)/round(leaf_area,3)) if round(leaf_area,3)>0. else 0.
            ratio_green = min(1, round(leaf_green_area,3)/round(leaf_area,3)) if round(leaf_area,3)>0. else 0.
            lesion_density_by_fungus = {}
            for fungus in fungi:
                nb_lesions = sum([les.nb_lesions for les in leaf_lesions if les.fungus.name==fungus])
                lesion_density_by_fungus[fungus] = nb_lesions/leaf_area if leaf_area>0 else 0.


            if round(ratio_green*(1-ratio_les_surface), 10) == 0.:
                for vid in set(leaf) & set(dispersal_units):
                    DUs[vid] = []
                    dispersal_units[vid] = []
            else:
#                if blade==42 and g.node(leaf[0]).dd_sequence[0]>384.:
#                    import pdb
#                    pdb.set_trace()
#                
                for vid in set(leaf) & set(dispersal_units):
                    all_deposit = []
                    for fungus in fungi:
                        dus_to_keep = []
                        dus = []
                        for du in dispersal_units[vid]:
                            if du.is_active and du.fungus.name==fungus:
                                if du.status == 'deposited':
                                    dus_to_keep.append(du)
                                elif du.status == 'emitted':
                                    dus.append(du)
                            
                        if len(dus)>0.:
                            group_dus = dus[0].fungus.group_dus
                            if group_dus == True:
                                total_nb_dus = dus[0].nb_dispersal_units
                            else:
                                total_nb_dus = len(dus)
#                                nb_on_healthy = int(total_nb_dus*ratio_green*(1-ratio_les_surface))
                            nb_on_healthy = int(total_nb_dus*ratio_green)
                            if lesion_density_by_fungus[fungus]>=self.max_lesion_density:
                                nb_on_healthy = 0.
                            
                            # Temp
                            age_lim = 50.
                            if self.age_infection and blade in ages and ages[blade]<age_lim:
                                nb_on_healthy *= 0.
                            #
                                
                            if nb_on_healthy > 0:
                                if group_dus == True:
                                    dus[0].nb_dispersal_units = nb_on_healthy
                                    dus[0].set_status(status = 'deposited')
                                    dus = [dus[0]]
                                else:
                                    dus = [du for i, du in enumerate(dus) 
                                            if i <= nb_on_healthy]
                                    for du in dus:
                                        du.set_status(status = 'deposited')
                                all_deposit += dus_to_keep + dus
                        else:
                            all_deposit += dus_to_keep
                    DUs[vid] = all_deposit
                    dispersal_units[vid] = all_deposit
#                    if len(all_deposit)>0:
#                        import pdb
#                        pdb.set_trace()