""" Gather different strategies for models that coordinate growth of lesion on 
    leaves according to area available (i.e. competition).

"""

# Imports #########################################################################

# With no priority between lesions ################################################
class NoPriorityGrowthControl:
    """ Template class for a model of competition between lesions for leaf area.
    
    A class for a model of growth control must include a method 'control'. In this
    example, no priority is defined between lesions in different states. On the 
    same leaf, the older lesions grow before the others only because of the way the 
    calculation is performed.
    
    """   
    def control(self, g, label='LeafElement'):
        """ Example to limit lesion growth to the healthy area on leaves.
        
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
            
        ..todo:: 
        - Check for modularity
        - Assert that leaf surface is not negative
        
        """       
        lesions = {k:v for k,v in g.property('lesions').iteritems() if len(v)>0.}
        areas = g.property('area')
        green_areas = g.property('green_area')
        senesced_areas = g.property('senesced_area')
        labels = g.property('label')
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith(label)]
            leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], [])
            les_surf = sum([les.surface for les in leaf_lesions])
            leaf_area = sum([areas[lf] for lf in leaf])
            leaf_green_area = sum([green_areas[lf] for lf in leaf])
            leaf_senesced_area = sum([senesced_areas[lf] for lf in leaf])
            ratio_green = min(1., leaf_green_area/leaf_area) if leaf_area>0. else 0.
            green_lesion_area = les_surf * ratio_green if leaf_senesced_area > les_surf else les_surf - leaf_senesced_area
            leaf_healthy_area = leaf_area - (leaf_senesced_area + green_lesion_area)
            leaf_healthy_area = max(0., round(leaf_healthy_area, 10))   
            total_demand = sum(l.growth_demand for l in leaf_lesions)
            if total_demand > leaf_healthy_area:
                for l in leaf_lesions:
                    growth_offer = round(leaf_healthy_area * l.growth_demand / total_demand, 14)
                    l.control_growth(growth_offer=growth_offer)
            else:
                for l in leaf_lesions:
                    growth_offer = l.growth_demand
                    l.control_growth(growth_offer = l.growth_demand)

class PriorityGrowthControl:
    """ 
    """   
    def control(self, g, label='LeafElement'):
        """ 
        """
        lesions = {k:v for k,v in g.property('lesions').iteritems() if len(v)>0.}
        areas = g.property('area')
        green_areas = g.property('green_area')
        senesced_areas = g.property('senesced_area')
        labels = g.property('label')

        # Select all the leaves
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith(label)]
            leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], []) 
            les_surf = sum([les.surface for les in leaf_lesions])
            leaf_area = sum([areas[lf] for lf in leaf])
            leaf_green_area = sum([green_areas[lf] for lf in leaf])
            leaf_senesced_area = sum([senesced_areas[lf] for lf in leaf])
            ratio_green = min(1., leaf_green_area/leaf_area) if leaf_area>0. else 0.
            green_lesion_area = les_surf * ratio_green if leaf_senesced_area > les_surf else les_surf - leaf_senesced_area
            leaf_healthy_area = leaf_area - (leaf_senesced_area + green_lesion_area)
            leaf_healthy_area = max(0., round(leaf_healthy_area, 10))   
            total_demand = sum(l.growth_demand for l in leaf_lesions)
            if total_demand > leaf_healthy_area:
                prior_lesions = [l for l in leaf_lesions if l.status>=l.fungus.CHLOROTIC]
                non_prior_lesions = [l for l in leaf_lesions if l.status<l.fungus.CHLOROTIC]
                prior_demand = sum(l.growth_demand for l in prior_lesions)
                if prior_demand > leaf_healthy_area:
                    for l in non_prior_lesions:
                        l.control_growth(growth_offer=0.)
                    for l in prior_lesions:
                        growth_offer = round(leaf_healthy_area * l.growth_demand / prior_demand, 14)
                        l.control_growth(growth_offer=growth_offer)
                else:
                    for l in prior_lesions:
                        growth_offer = l.growth_demand
                        l.control_growth(growth_offer=growth_offer)
                    non_prior_demand = sum(l.growth_demand for l in non_prior_lesions)
                    assert non_prior_demand >= (leaf_healthy_area-prior_demand)
                    for l in non_prior_lesions:
                        growth_offer = round((leaf_healthy_area-prior_demand) * l.growth_demand / non_prior_demand, 14)
                        l.control_growth(growth_offer=growth_offer)
            else:
                for l in leaf_lesions:
                    growth_offer = l.growth_demand
                    l.control_growth(growth_offer=growth_offer)

class GrowthControlVineLeaf:
    """ Class for growth control used when the phyto-element is a vine leaf.
    
    """   
    def control(self, g, label='lf'):
        """ ELimit lesion growth to the healthy area on vine leaves.
        
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
        healthy_areas = g.property('healthy_area')
        labels = g.property('label')
        
        # Select all the leaves
        vids = (v for v,l in labels.iteritems() if l.startswith(label))
        for leaf in vids:
            try:
                leaf_healthy_area = healthy_areas[leaf]
            except:
                raise NameError('Set healthy area on the MTG before simulation' '\n' 
                                'See alinea.alep.architecture > set_healthy_area')
            
            leaf_lesions = [l for l in lesions.get(leaf,[]) if l.is_active]
            total_demand = sum(l.growth_demand for l in leaf_lesions)
            
            if total_demand > leaf_healthy_area:
                for l in leaf_lesions:
                    growth_offer = leaf_healthy_area * l.growth_demand / total_demand
                    l.control_growth(growth_offer=growth_offer)
            else:
                for l in leaf_lesions:
                    growth_offer = l.growth_demand
                    l.control_growth(growth_offer=growth_offer)