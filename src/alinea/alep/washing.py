""" Gather different strategies for modeling the washing of fungus dispersal units 
    on leaves by rain.

"""

# Imports #########################################################################
import numpy

# Rapilly washing #################################################################

class RapillyWashing:
    """ Template class for a model of washing that complies with the guidelines of Alep.
    
    A class for a model of dispersal must include a method 'compute_washing_rate'.
    In this example, the rate is computed as in Rapilly and Jolivet, 1976 :
    'Construction d'un modele (episept) permettant la simulation d'une epidemie de 
    Septoria nodorum sur ble'
    
    """       
    def compute_washing_rate(self, g, global_rain_intensity, label='LeafElement'):
        """ Compute the washing rate on each leaf element.
        
        The washing rate is a function of rain duration and intensity.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil).
            Leaves must know the 'rain_intensity' and the 'rain_duration'
        global_rain_intensity: float
            Rain intensity over the canopy to trigger washing
        label: str
            Label of the part of the MTG concerned by the calculation
            
        Returns
        -------
        None
            Update directly the MTG
        """
        if global_rain_intensity > 0.:
            vids = [n for n in g if g.label(n).startswith(label)]
            for v in vids : 
                leaf = g.node(v)
                rain_int = leaf.rain_intensity
                rain_dur = leaf.rain_duration
                # healthy_area = leaf.healthy_area
                # leaf.washing_rate = max(0, min(1, rain_int/(healthy_area+rain_int)*rain_dur))
                area = leaf.area
                if leaf.geometry!=None:
                    leaf.washing_rate = max(0, min(1, rain_int/(area+rain_int)*rain_dur))           
                else:
                    leaf.washing_rate = 0.
            
            
    # TODO : Check if the model of Rapilly is properly implemented above.