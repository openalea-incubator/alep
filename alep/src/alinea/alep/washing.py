""" Gather existing models of washing of fungus dispersal units on leaves by rain.

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
    def __init__(self):
        self.rain_duration = 0.
        self.rain_intensity = {}
        
    def compute_washing_rate(self, g, global_rain_intensity, label='LeafElement'):
        """ Compute the washing rate on each leaf element.
        
        The washing rate is a function of rain duration and intensity.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
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
            self.rain_duration += 1
            
            vids = [n for n in g if g.label(n).startswith(label)]
            for v in vids : 
                leaf = g.node(v)
                if not v in self.rain_intensity:
                    self.rain_intensity[v] = []
                self.rain_intensity[v].append(leaf.rain_intensity)
                leaf.washing_rate = 0.
        else:
            vids = [n for n in g if g.label(n).startswith(label)]
            for v in vids : 
                leaf = g.node(v)
                if v in self.rain_intensity.keys():
                    mean_rain_intensity = numpy.mean(self.rain_intensity[v])
                    leaf.washing_rate = max(0,min(1, mean_rain_intensity / (leaf.healthy_surface + mean_rain_intensity)*self.rain_duration))
                    self.rain_intensity.pop(v)
                else:
                    leaf.washing_rate = 0.
            self.rain_duration = 0.
            
    # TODO : Check if the model of Rapilly is properly implemented above.