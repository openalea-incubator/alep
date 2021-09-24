""" (deprecated) dispersal strategies for septo3d (c implementation) model
"""

from alinea.alep.architecture import get_total_leaf_area
from alinea.alep.septoria import SeptoriaDU

# Septoria rain emission ##########################################################
class SeptoriaRainEmission:
    """ Template class for a model of emission of dispersal units by rain 
        that complies with the guidelines of Alep.
    
    Emission is the first step of dispersal, before transport. A class for a model 
    of emission must include a method 'get_dispersal_units'. In this example, the 
    number of dispersal units emitted is calculated according to the equations of
    Rapilly and Jolivet (1976) and interacts with specific septoria lesions.
    """
    
    def __init__(self, domain_area=None):
        """ Initialize the model with fixed parameters.
        
        Parameters
        ----------
        domain_area: float
            Domain area of the canopy stand
        """
        self.domain_area = domain_area
        
    def get_dispersal_units(self, g, fungus_name="septoria", 
                            label='LeafElement', weather_data=None,
                            domain_area = None, k_wheat = 0.65):
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
        if domain_area is None:
            domain_area = self.domain_area
        
        # Compute total leaf area
        total_area = get_total_leaf_area(g, label=label)
        
        # Compute intercept with Beer Lambert
        intercept = 1 - exp(-k_wheat*total_area*10**-4/self.domain_area)
        # intercept = 1 - exp(-k_wheat*total_area*10**-4)
        
        # Get lesions
        les = {k:[l for l in v if l.fungus.name is fungus_name and l.is_sporulating()] 
                    for k, v in g.property('lesions').items()} 
        
        # Compute total sporulating area
        total_spo = sum([l.surface_spo for v in list(les.values()) for l in v])
        tot_fraction_spo = total_spo/total_area if total_area>0. else 0.
        
        DU = {}
        for vid, l in les.items():
            for lesion in l:
                # Compute number of dispersal units emitted by lesion
                leaf = g.node(vid)
                total_DU_leaf = 0.36 * 6.19e7 * intercept * tot_fraction_spo *\
                                leaf.rain_intensity * domain_area
                
                if lesion.is_stock_available(leaf):
                    initial_stock = lesion.stock_spores
                    stock_available = int(lesion.stock_spores*2/3.)
                    contribution = lesion.surface_spo/total_spo if total_spo>0. else 0.
                    nb_DU_lesion = int(contribution * total_DU_leaf)
                    
                    # Distribute spores into dispersal units
                    nb_spores_by_DU = 10
                    nb_DU_lesion = min(nb_DU_lesion, stock_available/nb_spores_by_DU)
                    
                    SeptoriaDU.fungus = lesion.fungus
                    emissions = [SeptoriaDU(nb_spores = nb_spores_by_DU, status='emitted', 
                                 position=lesion.position) for i in range(nb_DU_lesion)]

                    # Update stock of spores of the lesion
                    nb_spores_emitted = nb_DU_lesion*nb_spores_by_DU
                    lesion.reduce_stock(nb_spores_emitted)
                    
                    # Update empty surface on lesion
                    lesion.update_empty_surface(nb_spores_emitted, initial_stock)
                    
                    if vid not in DU:
                        DU[vid] = []
                    DU[vid] += emissions

class BenchSeptoriaRainEmission:
    """ Template class for a model of emission of dispersal units by rain 
        that complies with the guidelines of Alep only for benchmark test.
    """
    
    def __init__(self, domain_area=None):
        """ Initialize the model with fixed parameters.
        
        Parameters
        ----------
        domain_area: float
            Domain area of the canopy stand
        """
        self.domain_area = domain_area
        
    def get_dispersal_units(self, g, fungus_name="septoria", label='LeafElement'):
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
        k_wheat = 0.65
        
        # Compute total leaf area
        total_area = get_total_leaf_area(g, label=label)
        
        # Compute intercept with Beer Lambert
        intercept = 1 - exp(-k_wheat*total_area*10**-4/self.domain_area)
        # intercept = 1 - exp(-k_wheat*total_area*10**-4)
        
        # Get lesions
        les = {k:[l for l in v if l.fungus.name is fungus_name and l.is_sporulating()] 
                    for k, v in g.property('lesions').items()} 
        
        # Compute total sporulating area
        total_spo = sum([l.surface_spo for v in list(les.values()) for l in v])
        # tot_fraction_spo = total_spo/(total_area/self.domain_area) if (total_area/self.domain_area)>0. else 0.
        tot_fraction_spo = total_spo/total_area if total_area>0. else 0.

        DU = {}
        for vid, l in les.items():
            for lesion in l:
                # Compute number of dispersal units emitted by lesion
                leaf = g.node(vid)
                total_DU_leaf = 0.36 * 6.19e7 * intercept * tot_fraction_spo * leaf.rain_intensity * self.domain_area

                initial_stock = lesion.stock_spores
                stock_available = int(lesion.stock_spores*2/3.)
                contribution = lesion.surface_spo/total_spo if total_spo>0. else 0.
                nb_DU_lesion = int(contribution * total_DU_leaf)
                
                # Distribute spores into dispersal units
                nb_spores_by_DU = 10
                nb_DU_lesion = min(nb_DU_lesion, stock_available/nb_spores_by_DU)
                
                SeptoriaDU.fungus = lesion.fungus
                emissions = [SeptoriaDU(nb_spores = nb_spores_by_DU, status='emitted', 
                             position=lesion.position) for i in range(nb_DU_lesion)]

                # Update stock of spores of the lesion
                nb_spores_emitted = nb_DU_lesion*nb_spores_by_DU
                lesion.reduce_stock(nb_spores_emitted)
                
                # Update empty surface on lesion
                lesion.update_empty_surface(nb_spores_emitted, initial_stock)
                
                if vid not in DU:
                    DU[vid] = []
                DU[vid] += emissions
        return DU
    