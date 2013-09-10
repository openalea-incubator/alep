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
        lesions = g.property('lesions')
        labels = g.property('label')
        healthy_areas = g.property('healthy_area')           
        
        # Select all the leaves
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            try:
                leaf = [vid for vid in g.components(blade) if labels[vid].startswith(label)]
                leaf_healthy_area = sum(healthy_areas[lf] for lf in leaf)
            except:
                raise NameError('Set healthy area on the MTG before simulation' '\n' 
                                'See alinea.alep.architecture > set_healthy_area')
                
            leaf_lesions = [l for lf in leaf for l in lesions.get(lf,[]) if l.is_active]
            total_demand = sum(l.growth_demand for l in leaf_lesions)
            
            if total_demand > leaf_healthy_area:
                for l in leaf_lesions:
                    growth_offer = leaf_healthy_area * l.growth_demand / total_demand
                    l.control_growth(growth_offer=growth_offer)
                # for lf in leaf:
                    # # Update healthy area
                    # healthy_areas[lf] = 0.
            else:
                for l in leaf_lesions:
                    growth_offer = l.growth_demand
                    l.control_growth(growth_offer=growth_offer)
                # for lf in leaf:
                    # gd = sum(l.growth_demand for l in lesions.get(lf,[]) if l.is_active)
                    # # Update healthy area
                    # healthy_areas[lf] -= gd
                    # /!\ WARNING /!\
                    # This method leads to local healthy surfaces < 0 for leaf elements.
                    # But the global healthy area on the entire leaf stays >= 0.
                    # TODO : if the surface of a phyto-element is < 0, report the loss
                    # to its neighbour ?
        
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
                # healthy_areas[leaf] = 0.
            else:
                for l in leaf_lesions:
                    growth_offer = l.growth_demand
                    l.control_growth(growth_offer=growth_offer)

                gd = sum(l.growth_demand for l in lesions.get(leaf,[]) if l.growth_is_active)
                # healthy_areas[leaf] -= gd