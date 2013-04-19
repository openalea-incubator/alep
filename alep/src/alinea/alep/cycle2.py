# -*- coding: latin1 -*- 

##
##

""" Define the interface for the fungal classes of dispersal units, lesions,
    rings and parameters.

"""
# Dispersal unit ##################################################################
class DispersalUnit(object):
    """ Generic class for a dispersal unit.
    
    """
    def __init__(self, position=None, nb_spores=None, status=None):
        """ Initialize the dispersal unit.
        
        Parameters
        ----------
        position: non defined
            Position of the dispersal unit on the phyto-element
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        status: str
            'emitted' or 'deposited'
        
        Returns
        -------
            None
        """
        self.position = position
        self.nb_spores = nb_spores
        self.status = status
        self.is_active = True
    
    def disable(self):
        """ disable a dispersal unit.
        
        """
        self.is_active = False
    
    def deposited(self):
        """ Change the status of the spore to 'deposited'.
        
        """
        self.status = 'deposited'
    
    def set_position(self, position=None):
        """ Set the position of the DU to position given in argument.
        
        Parameters
        ----------
        position: Non defined
            Position of the DU.
        """
        self.position = position
        
    def create_lesion(self, leaf=None):
        """ Create a new lesion of fungus and disable dispersal unit.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        
        Returns
        -------
            None
        
        """
        les = self.fungus(position=self.position, nb_spores=self.nb_spores)
        try:
            leaf.lesions.append(les)
        except:
            leaf.lesions = [les]
        self.disable()
             
# Lesions #################################################################################

#X class LesionFactory(object):
#X     """
#X     """
#X 
#X     def __init__(self, fungus):
#X         """ Initialize the lesion. 
#X         
#X         :Parameters:
#X           - `fungus` (function): returns a class of specific parameters for 
#X           the chosen fungus (e.g. 'septoria()' or 'powderyMildew()').
#X         """
#X         self.fungus = fungus
#X         
#X     def instantiate(self, nb_spores, position) :
#X         """ instantiate a Lesion
#X      
#X         :Parameters:
#X           - `nb_spores` (int): 
#X         """
#X         l = Lesion(self.fungus, nb_spores, position)
#X         return l
#X         
#X     def instantiate_at_stage(self, nb_spores, position) :
#X         """ force the instantiation of a Lesion at a given stage"""
#X         l = Lesion(self.fungus, nb_spores, position)
#X         #to do : deal with spores
#X         return l
#X 

class Lesion(object):
    """ Define a lesion interface.

    To implement a lesion, you need to implement the following methods:
        - update()
        - growth_control()
        - stock_spores()
        - emission()
        - disable()        
    
    And the lesion have to answer to a set of queries:
        - is_dead
        - surface
        - status
        - age
        - age_physio
        
    ..todo:: improve header and doc
    """
    fungus = None
    def __init__(self, position=None, nb_spores=None):
        """ Initialize the lesion. 
        
        Parameters
        ----------
        position: non defined
            Position of the dispersal unit on the phyto-element
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        
        """
        self.is_active = True
        self.nb_spores = nb_spores
        self.position = position
        self.stock_spores = []
        self.emissions = []
    
    def disable(self):
        """ Disable the lesion.
        
        Set the activity of the lesion to False and its growth demand to 0.
        
        Parameters
        ----------
            None
        """
        self.is_active = False
        self.growth_demand = 0.
    
# Rings ##################################################################################
class Ring(object):
    """ Ring of Lesion at a given age.
    """

# Fungus Parameters (e.g. .ini): config of the fungus ####################################
class Parameters(object):
    def write(self, filename):
        pass
    
    @staticmethod
    def read(filename):
        pass
