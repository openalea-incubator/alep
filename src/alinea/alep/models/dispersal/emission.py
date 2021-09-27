""" Gather different strategies for modeling emission of fungus propagules, 
    first step of dispersal.
"""

from alinea.alep.dispersal import Emission, get_sporulating_lesions

# Imports #########################################################################

#from alinea.alep.powdery_mildew import PowderyMildewDU
from math import exp

# Rain emission ##########################################################
class RainEmission(Emission):
    """ Template class for a model of emission of dispersal units by rain 
        that complies with the guidelines of Alep.
    
    Emission is the first step of dispersal, before transport. A class for a model 
    of emission must include a method 'get_dispersal_units'. In this example, the 
    number of dispersal units emitted is calculated according to the equations of
    Rapilly and Jolivet (1976) and interacts with specific septoria lesions.
    """
    
    def __init__(self, pFA:float=6.19e7, pEclin:float=0.36, Imin:float=0.5):
        
        """Initialize the model with fixed parameters.

        Parameters
        ----------
        pFA : float, 
            proportionality between rain intensity and number of splash droplets emited (m-2), by default 6.19e7
        pEclin : float
            proportion of emitted infectious droplets that contain enough spores and that do not evaporate, by default 0.36
        Imin : float, optional
            minimal rain intensity, by default 0.5
        """        
        super(RainEmission, self).__init__()
        self.pFA = pFA
        self.pEclin = pEclin
        self.Imin = Imin
        
    def get_dispersal_units(self, lesions:dict, exposed_sporulating_areas:dict, fungus_name:str=None, rain_intensity:float=1) -> dict:
        """Compute emission of dispersal units by rain

        Parameters
        ----------
        lesions : dict
            MTG property containing lesions {vid: [lesion,...], ...}
        exposed_sporulating_areas : dict
            MTG property with sporulating areas exposed to rain. Areas should be in square meter
            this quantity can be determined for a vid using a projection models. {vid: [area,...],...}
        fungus_name : str, optional
            name of fungus, by default None
        rain_intensity : int, optional
            Rain intensity= ratio of rain quantity by rain duration, by default 1

        Returns
        -------
        dict
            dispersal_units : dict([leaf_id, list of dispersal units])
            Dispersal units emitted by source lesions
        """        
        
        if rain_intensity <= self.Imin:
            return {}
        
        # Get lesions
        sporulating_lesions = get_sporulating_lesions(lesions, fungus_name)
        assert all([k in exposed_sporulating_areas for k in sporulating_lesions])
        
        DU = {k:[] for k in sporulating_lesions}
        for vid, l in sporulating_lesions.items():
            for i,lesion in enumerate(l):
                nb_du = self.pFA * rain_intensity * self.pEclin * exposed_sporulating_areas[vid][i]
                if nb_du > 0:
                    DU[vid].append(lesion.emission(nb_DU=nb_du))
        return DU
 
 
# Septoria rain emission ##########################################################
class PowderyMildewWindEmission:
    """ Template class for a model of emission of dispersal units by wind 
        that complies with the guidelines of Alep.
    
    Emission is the first step of dispersal, before transport. A class for a model 
    of emission must include a method 'get_dispersal_units'. In this example, the 
    number of dispersal units emitted is calculated according to the equations of
    Willocquet and Clerjeau (1998).
    """
    
    def __init__(self):
        """ Initialize the model with fixed parameters.
        """
        pass
        
    def get_dispersal_units(self, g, fungus_name="powdery_mildew", label='lf',
                            b = -5.8, r = 0.41):
        """ Compute emission of dispersal units by rain splash on wheat.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        fungus_name: str
            Name of the fungus
                    
        Returns
        -------
        dispersal_units : dict([leaf_id, list of dispersal units])
            Dispersal units emitted by leaf.
        """
        b = -5.8
        r = 0.41
        
        # Get lesions
        lesions = {k:[l for l in les if l.fungus.name is fungus_name and l.is_sporulating()] 
                    for k, les in g.property('lesions').items()} 
        
        DU = {}
        for vid, l in lesions.items():
            leaf = g.node(vid)
            wind_speed = leaf.wind_speed
            for lesion in l:
                emissions = []
                stock = lesion.stock_spores
                if stock > 0.:
                    # Dispersal rate
                    dispersal_rate = exp(r*wind_speed+b) / (1+exp(r*wind_speed+b))
                    # Number of spores emitted : One only spore by DU
                    nb_DU_emitted = int(dispersal_rate * stock)
                    
                    if nb_DU_emitted > 0 :
                        # Return emissions
                        PowderyMildewDU.fungus = lesion.fungus
                        emissions = [PowderyMildewDU(nb_spores=1, status='emitted',
                                                position=lesion.position) for i in range(nb_DU_emitted)]
                        # Update stock of spores                       
                        lesion.reduce_stock(nb_spores_emitted = nb_DU_emitted)

                        if vid not in DU:
                            DU[vid] = []
                        DU[vid] += emissions
        return DU