""" Gather different strategies for models that coordinate growth of lesion on 
    leaves according to area available (i.e. competition).

"""

# Imports #########################################################################
import numpy as np
import random as rd

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
        geom = g.property('geometry')
        labels = g.property('label')
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) 
                    if labels[vid].startswith(label) and vid in geom]
            if len(leaf) > 0.:
                leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], [])
                les_surf = sum([les.surface for les in leaf_lesions])
                leaf_area = sum([areas[lf] for lf in leaf])
                leaf_green_area = sum([green_areas[lf] for lf in leaf])
                leaf_senesced_area = sum([senesced_areas[lf] for lf in leaf])
                ratio_green = min(1., leaf_green_area/leaf_area) if leaf_area>0. else 0.
                green_lesion_area = les_surf * ratio_green if leaf_senesced_area > les_surf else les_surf - leaf_senesced_area
                leaf_healthy_area = leaf_area - (leaf_senesced_area + green_lesion_area)
                leaf_healthy_area = max(0., round(leaf_healthy_area, 10))
                # TEMP
#                leaf_healthy_area *= (1 - les_surf/leaf_area)
                #

                total_demand = sum(l.growth_demand for l in leaf_lesions)
                # formule demande en entree avec taille des lesions moyenne
                if total_demand > leaf_healthy_area:
                    for l in leaf_lesions:
                        growth_offer = round(leaf_healthy_area * l.growth_demand / total_demand, 14)
                        l.control_growth(growth_offer=growth_offer)
#                        l.control_growth(growth_offer=growth_offer*(0.88 - green_lesion_area/leaf_green_area))
                else:
                    for l in leaf_lesions:
                        growth_offer = l.growth_demand
                        l.control_growth(growth_offer = growth_offer)
#                        l.control_growth(growth_offer=growth_offer*(0.88 - green_lesion_area/leaf_green_area))

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
        geom = g.property('geometry')

        # Select all the leaves
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) 
                    if labels[vid].startswith(label) and vid in geom]
            if len(leaf) > 0.:
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
        """ Limit lesion growth to the healthy area on vine leaves.
        
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

# Geometric competition with circular lesions #################################
class GeometricCircleCompetition:
    """ Model of competition between circular lesions for leaf area.
    
        In this model, a Poisson law is used to simulate the congestion between
        circular lesions on leaves in healthy part of the leaf
    
    """
    def true_area_impacted(self, nb_lesions, available_area, mean_lesion_size):
        if available_area>0.:
            return available_area*(1-np.minimum(1.,np.exp(-nb_lesions*mean_lesion_size/available_area)))
        else:
            return 0.
        
    def control(self, g, label='LeafElement'):
        """ Limit lesion growth to healthy area on leaves and simulate 
            congestion between circular lesions.
        """       
        lesions = {k:v for k,v in g.property('lesions').iteritems() if len(v)>0.}
        areas = g.property('area')
        green_areas = g.property('green_area')
        geom = g.property('geometry')
        labels = g.property('label')
        potential_lesion_areas = g.property('potential_lesion_area')
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) 
                    if labels[vid].startswith(label) and vid in geom]

            # Save potential lesion area (same on all vids of leaf)            
            pots = list(set(potential_lesion_areas.keys()) & set(leaf))
            if len(pots)>0:
                potential_lesion_area = potential_lesion_areas[pots[0]]
            else:
                potential_lesion_area = 0.
                
            if len(leaf) > 0:
                leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], [])
                nb_lesions = sum([les.nb_lesions for les in leaf_lesions])
                if nb_lesions>0.:
                    les_surf = sum([les.surface for les in leaf_lesions])
                    total_demand = sum(l.growth_demand for l in leaf_lesions)
                    leaf_area = sum([areas[lf] for lf in leaf])
                    leaf_green_area = sum([green_areas[lf] for lf in leaf])
                    ratio_green = min(1., leaf_green_area/leaf_area) if leaf_area>0. else 0.
                    nb_les_on_green = float(nb_lesions) * ratio_green
                   
                    potential_lesion_area += total_demand
                    true_area = self.true_area_impacted(nb_les_on_green, 
                                                        leaf_green_area, 
                                                        potential_lesion_area*ratio_green/nb_les_on_green)

                    offer = true_area - les_surf

                    for l in leaf_lesions:
                        growth_offer = l.growth_demand * offer/total_demand if total_demand>0. else 0.
                        l.control_growth(growth_offer = growth_offer)
                
                potential_lesion_areas.update({vid:potential_lesion_area for vid in leaf})
                
# Growth control between 2 diseases ###########################################
import random as rd
class SeptoRustCompetition:
    """ Class that control growth of lesions of brown rust and septoria in
        competition. 
    """
    def __init__(self, SeptoModel, RustModel):
        self.SeptoModel = SeptoModel
        self.RustModel = RustModel
        
    def control(self, g, label = 'LeafElement'):
        lesions = {k:v for k,v in g.property('lesions').iteritems() if len(v)>0.}
        areas = g.property('area')
        green_areas = g.property('green_area')
        labels = g.property('label')
        geom = g.property('geometry')
        # Select all the leaves
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) 
                    if labels[vid].startswith(label) and vid in geom]
            if len(leaf) > 0.:
                area = sum([areas[lf] for lf in leaf])
                green_area = sum([green_areas[lf] for lf in leaf])
                leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], [])
                
                # Separate lesions in two groups according to priority level
                prio_les = []
                other_les = []
                for l in leaf_lesions:
                    if isinstance(l, self.SeptoModel) and l.status >= l.fungus.CHLOROTIC:
                        prio_les.append(l)
                    else:
                        other_les.append(l)

                # Grow prioritary lesions over other ones
                if len(prio_les)>0:
                    s_full = sum([l.surface for l in prio_les])
                    offer_prio = max(green_area - s_full*green_area/area, 0)
                    demand_prio = sum(l.growth_demand for l in prio_les)
                    ratio = demand_prio/offer_prio if offer_prio > 0. else 0.
                    if ratio < 1 and ratio > 0:
                        # Surface is not limiting for prioritary lesions
                        demand_others = 0.
                        s_others = 0.
                        # Remove some other lesions that are overgrown
                        nb_dead = int(round(len(other_les)*ratio))
                        deads = rd.sample(other_les, nb_dead)
                        import pdb
                        pdb.set_trace()                        
                        for l in deads:
                            l.disappear()
                            other_les.remove(l)
                        for l in other_les:
                            demand_others += l.growth_demand
                            s_others += l.surface
                        
#                        for l in other_les:
#                            if np.random.random() < ratio:
#                                l.disappear()
#                                other_les.remove(l)
#                            else:
#                                demand_others += l.growth_demand
#                                s_others += l.surface
                        # All prioritary lesions can grow as expected
                        for l in prio_les:
                            l.control_growth(l.growth_demand)
                        offer_others = green_area - s_others*green_area/area - demand_prio
                    else:
                        # Space is limiting already for prioritary lesions
                        # Grow over all others and share space left
                        for l in other_les:
                            l.disappear()
                        offer_each = offer_prio / len(prio_les)
                        for l in prio_les:
                            l.control_growth(offer_each)
                        other_les = []
                else:                   
                    demand_others = sum(l.growth_demand for l in other_les)
                    s_others = sum([l.surface for l in other_les])
                    offer_others = green_area - s_others*green_area/area
                
                # Distribute surface left for other lesions with lower
                # level of priority
                if len(other_les)>0.:                
                    if demand_others < offer_others:
                        for l in other_les:
                            l.control_growth(l.growth_demand)
                    else:
                        offer_each = offer_others / len(other_les)
                        for l in other_les:
                            l.control_growth(offer_each)
#                        
#                import pdb
#                pdb.set_trace()
#                
                