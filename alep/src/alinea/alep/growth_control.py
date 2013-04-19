""" Gather different strategies for models that coordinate growth of lesion on 
    leaves according to surface available (i.e. competition).

"""

# Imports #########################################################################

# With no priority between lesions ################################################

class NoPriorityGrowthControl:
    """ Template class for a model of competition between lesions for leaf surface.
    
    A class for a model of growth control must include a method 'control'. In this
    example, no priority is defined between lesions in different states. On the 
    same leaf, the older lesions grow before the others only because of the way the 
    calculation is performed.
    
    """   
    def control(self, g, label='LeafElement'):
        """ Example to limit lesion growth to the healthy surface on leaves.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        label: str
            Label of the part of the MTG concerned by the calculation
        
        Returns
        -------
        None
            Update directly the MTG
        """
        lesions = g.property('lesions')
        surfaces = g.property('surface')
        healthy_surfaces = g.property('healthy_surface')
        labels = g.property('label')

        # Select all the leaves
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith(label)]
            leaf_surface = sum(surfaces[lf] for lf in leaf)
            leaf_healthy_surface = sum(healthy_surfaces[lf] for lf in leaf)
            
            leaf_lesions = [l for lf in leaf for l in lesions.get(lf,[])]
            total_demand = sum(l.growth_demand for l in leaf_lesions)
            
            if total_demand > leaf_healthy_surface:
                for l in leaf_lesions:
                    growth_offer = leaf_healthy_surface * l.growth_demand / total_demand
                    l.control_growth(growth_offer=growth_offer)
                for lf in leaf:
                    healthy_surfaces[lf] = 0.
            else:
                for l in leaf_lesions:
                    growth_offer = l.growth_demand
                    l.control_growth(growth_offer=growth_offer)
                for lf in leaf:
                    gd = sum(l.growth_demand for l in lesions.get(lf,[]))
                    healthy_surfaces[lf] -= gd
                    
