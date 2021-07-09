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
        lesions = {k:v for k,v in g.property('lesions').items() if len(v)>0.}
        areas = g.property('area')
        green_areas = g.property('green_area')
        senesced_areas = g.property('senesced_area')
        geom = g.property('geometry')
        labels = g.property('label')
        bids = (v for v,l in labels.items() if l.startswith('blade'))
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
                    for l in leaf_lesions:
                        growth_offer = round(leaf_healthy_area * l.growth_demand / total_demand, 14)
                        l.control_growth(growth_offer=growth_offer)
                else:
                    for l in leaf_lesions:
                        growth_offer = l.growth_demand
                        l.control_growth(growth_offer = growth_offer)

class PriorityGrowthControl:
    """ 
    """   
    def control(self, g, label='LeafElement'):
        """ 
        """
        lesions = {k:v for k,v in g.property('lesions').items() if len(v)>0.}
        areas = g.property('area')
        green_areas = g.property('green_area')
        senesced_areas = g.property('senesced_area')
        labels = g.property('label')
        geom = g.property('geometry')

        # Select all the leaves
        bids = (v for v,l in labels.items() if l.startswith('blade'))
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
        vids = (v for v,l in labels.items() if l.startswith(label))
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
class GeometricPoissonCompetition:
    """ Model of competition between lesions for leaf area.
    
        In this model, a Poisson law is used to simulate the congestion between
        circular lesions on leaves in healthy part of the leaf
    
    """    
    def true_area_impacted(self, nb_lesions, available_area, mean_lesion_size):
        if available_area>0.:
            return available_area*(1-np.minimum(1.,np.exp(-nb_lesions*mean_lesion_size/available_area)))
        else:
            return 0.
        
    def manage_senescence_border(self, leaf, r, lesions, 
                                 senesced_lengths, lengths):
        s_ls = [senesced_lengths[lf] for lf in leaf]
        s_l = sum(s_ls)
        t_ls = [lengths[lf] for lf in leaf]
        t_l = sum(t_ls)
        new_length = t_l - r*(t_l - s_l)
        idx = -1
        while round(new_length,10) > 0:
            n_l = min(t_ls[idx],new_length)
            s_ls[idx] = n_l
            new_length -= n_l
            idx -= 1
        senesced_lengths.update({lf:s_ls[i] for i,lf in enumerate(leaf)})
        for lf in leaf:
            if lf in lesions:
                for l in lesions[lf]:
                    l.senescence_response(senesced_length=senesced_lengths[lf])
    
    def control(self, g, label='LeafElement'):
        """ Limit lesion growth to healthy area on leaves and simulate 
            congestion between circular lesions.
        """       
        lesions = {k:v for k,v in g.property('lesions').items() if len(v)>0.}
        green_areas = g.property('green_area')
        geom = g.property('geometry')
        labels = g.property('label')
        lengths = g.property('length')
        areas = g.property('area')
        senesced_lengths = g.property('senesced_length')
        bids = (v for v,l in labels.items() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) 
                    if labels[vid].startswith(label) and vid in geom]
            if len(leaf) > 0:
                leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], [])
                nb_lesions = sum([les.nb_lesions_non_sen for les in leaf_lesions])
                if nb_lesions>0.:
                    les_surf = 0.
                    pot_les_surf = 0.
                    les_surf_non_sen = 0.
                    total_demand = 0.
                    for les in leaf_lesions:
                        les_surf += les.surface
                        pot_les_surf += les.potential_surface
                        les_surf_non_sen += les.surface_non_senescent
                        total_demand += les.growth_demand
                    leaf_green_area = sum([green_areas[lf] for lf in leaf])
                    leaf_area = sum([areas[lf] for lf in leaf])
                    if round(les_surf_non_sen,16) < round(leaf_green_area, 16):
                        true_area = self.true_area_impacted(nb_lesions, 
                                                            leaf_green_area, 
                                                            pot_les_surf/nb_lesions)

                        offer = min(leaf_area-les_surf, true_area - les_surf_non_sen)
                    else:
                        offer = 0
                        r = leaf_green_area / les_surf_non_sen
                        if round(r,14) < 1:
                            self.manage_senescence_border(leaf, r, lesions,
                                                          senesced_lengths,
                                                          lengths)

                    if offer > 0:
                        for l in leaf_lesions:
                            growth_offer = l.growth_demand * offer/total_demand if total_demand>0. else 0.
                            l.control_growth(growth_offer = growth_offer)
                    else:
                        for l in leaf_lesions:
                            growth_offer = offer*les.surface/les_surf
                            l.control_growth(growth_offer = growth_offer)

# Growth control between 2 diseases ###########################################
class SeptoRustCompetition:
    """ Class that control growth of lesions of brown rust and septoria in
        competition. 
    """     
    def __init__(self):
        self.disease_area=0.
    
    def true_area_impacted(self, nb_lesions, available_area, mean_lesion_size):
        if available_area>0.:
            return available_area*(1-np.minimum(1.,np.exp(-nb_lesions*mean_lesion_size/available_area)))
        else:
            return 0.
    
    def manage_senescence_border(self, leaf, r, lesions, 
                                 senesced_lengths, lengths):
        s_ls = [senesced_lengths[lf] for lf in leaf]
        s_l = sum(s_ls)
        t_ls = [lengths[lf] for lf in leaf]
        t_l = sum(t_ls)
        new_length = t_l - r*(t_l - s_l)
        idx = -1
        while round(new_length,10) > 0:
            n_l = min(t_ls[idx],new_length)
            s_ls[idx] = n_l
            new_length -= n_l
            idx -= 1
        senesced_lengths.update({lf:s_ls[i] for i,lf in enumerate(leaf)})
        for lf in leaf:
            if lf in lesions:
                for l in lesions[lf]:
                    l.senescence_response(senesced_length=senesced_lengths[lf])
    
    def control(self, g, label = 'LeafElement'):
        """ Model of competition. 
        
        TODO : Review: lots of limit cases, difficult to understand. Aaaaargh
        """
        lesions = {k:v for k,v in g.property('lesions').items() if len(v)>0.}
        areas = g.property('area')
        green_areas = g.property('green_area')
        labels = g.property('label')
        geom = g.property('geometry')
        lengths = g.property('length')
        areas = g.property('area')
        senesced_lengths = g.property('senesced_length')
        bids = (v for v,l in labels.items() if l.startswith('blade'))
        
        # Temp
#        from alinea.adel.newmtg import adel_labels
#        a_labs = adel_labels(g)
#        ages = g.property('age')
        
        for blade in bids:
            leaf = [vid for vid in g.components(blade) 
                    if labels[vid].startswith(label) and vid in geom]
            if len(leaf) > 0:
                leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], [])
                leaf_area = sum([areas[lf] for lf in leaf])
                leaf_green_area = sum([green_areas[lf] for lf in leaf])
                s_prio = 0.
                s_non_prio = 0.
                s_prio_alive = 0.
                s_non_prio_alive = 0.
                s_prio_non_sen = 0.
                s_non_prio_non_sen = 0.
                s_pot_prio = 0.
                s_pot_non_prio = 0.
                demand_prio = 0.
                demand_non_prio = 0.
                offer_prio = 0.
                offer_non_prio = 0.
                prio_les = []
                non_prio_les = []
                nb_prio = 0.
                nb_non_prio = 0.
                nb_rust = 0.
                s_rust_non_sen = 0.
                for l in leaf_lesions:
                    if l.fungus.name=='septoria' and l.status >= l.fungus.CHLOROTIC:
                        s_prio += l.surface
                        s_prio_alive += l.surface_alive
                        s_prio_non_sen += l.surface_non_senescent
                        s_pot_prio += l.potential_surface
                        prio_les.append(l)
                        demand_prio += l.growth_demand
                        nb_prio += l.nb_lesions_non_sen
                    else:
                        s_non_prio += l.surface
                        s_non_prio_alive += l.surface_alive
                        s_non_prio_non_sen += l.surface_non_senescent
                        s_pot_non_prio += l.potential_surface
                        non_prio_les.append(l)
                        demand_non_prio += l.growth_demand
                        nb_non_prio += l.nb_lesions_non_sen
                        if l.fungus.name=='brown_rust':
                            nb_rust += l.nb_lesions_non_sen
                            s_rust_non_sen += l.surface_non_senescent
                
                if nb_prio>0:
                    if round(s_prio_non_sen, 16) < round(leaf_green_area, 16):
                        true_area_prio = self.true_area_impacted(nb_prio,
                                                                 leaf_green_area,
                                                                 s_pot_prio/nb_prio)  
                        if leaf_green_area==leaf_area:
                            offer_prio = min(leaf_area-s_prio_non_sen, true_area_prio - s_prio_non_sen)
                        elif leaf_area > s_prio:
                            offer_prio = min(leaf_area-s_prio, true_area_prio - s_prio_non_sen)
                        else:
                            offer_prio = 0.                            
                    else:
                        offer_prio = 0.
                        r = leaf_green_area / s_prio_non_sen if s_non_prio_non_sen>0 else 0.
                        if round(r,14) < 1:
                            self.manage_senescence_border(leaf, r, lesions,
                                                          senesced_lengths,
                                                          lengths)

                    offer_prio = max(0., min(offer_prio, demand_prio))
                    for l in prio_les:
                        offer_lesion = l.growth_demand*offer_prio/demand_prio if demand_prio>0. else 0.
                        l.control_growth(offer_lesion)
                        # Temp : Add limiting factor on sporulating capacity of septoria
                        # print "l.377 growth_control modif"
                        # if nb_rust>0:
                            # l.sporulating_capacity = s_rust_non_sen/leaf_green_area if leaf_green_area>0. else 0.

                if nb_non_prio>0:
                    nb_les = nb_prio + nb_non_prio
                    s_pot = s_pot_non_prio + s_pot_prio
                    if round(s_non_prio_non_sen, 16) < round(leaf_green_area, 16):
                        true_area_non_prio = self.true_area_impacted(nb_les,
                                                                     leaf_green_area, 
                                                                     s_pot/nb_les)
                        new_true_area_prio = s_prio_alive+offer_prio
                        offer_non_prio = min(leaf_area-s_non_prio,
                                             true_area_non_prio-s_non_prio_non_sen-new_true_area_prio)
                    else:
                        offer_non_prio = 0.
                        r = leaf_green_area / s_non_prio_non_sen if s_non_prio_non_sen>0 else 0.
                        if round(r,14) < 1:
                            self.manage_senescence_border(leaf, r, lesions,
                                                          senesced_lengths,
                                                          lengths)

                    if leaf_green_area<leaf_area:
                        max_offer = max(0, leaf_green_area - s_non_prio - s_prio_alive - offer_prio)
                        offer_non_prio = min(offer_non_prio, max_offer)
                    
                    if offer_non_prio>0:
                        offers = 0.
                        if offer_non_prio > demand_non_prio:
                            offer_non_prio = demand_non_prio
                        for l in non_prio_les:
                            offer_lesion = l.growth_demand * offer_non_prio/demand_non_prio if demand_non_prio>0. else 0.
                            offers += offer_lesion
                            l.control_growth(offer_lesion)
                    else:
                        offer_non_prio = max(offer_non_prio, -s_non_prio_non_sen)
                        offer_non_prio = max(offer_non_prio, -offer_prio)
                        for l in non_prio_les:
                            offer_lesion = offer_non_prio*l.surface_non_senescent/s_non_prio_non_sen if s_non_prio_non_sen>0. else 0.                        
                            l.control_growth(growth_offer=offer_lesion)
