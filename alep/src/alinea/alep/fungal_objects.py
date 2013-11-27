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
    def __init__(self, nb_spores=None, position=None, status=None):
        """ Initialize the dispersal unit.
        
        Parameters
        ----------
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        position: non defined
            Position of the dispersal unit on the phyto-element
        status: str
            'emitted' or 'deposited'
        
        Returns
        -------
            None
        """
        self.nb_spores = nb_spores
        self.position = position
        self.status = status
        self.is_active = True
        self.can_infect_at_position = None
    
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
    
    def infection_impossible_at_position(self):
        """ Turn the property 'can_infect_at_position' to False.
        """
        self.can_infect_at_position = False
    
    def infection_possible_at_position(self):
        """ Turn the property 'can_infect_at_position' to True.
        """
        self.can_infect_at_position = True
    
    def create_lesion(self, leaf=None):
        """ Create a new lesion of fungus and disable dispersal unit.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. area, green area, healthy area,
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
        position: not defined
            Position of the dispersal unit on the phyto-element
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        
        """
        # Total activity of the lesion (growth and ageing)
        self.is_active = True
        # Growth activity of the lesion
        self.growth_is_active = True
        # Growth demand of the lesion (cm2)
        self.growth_demand = None
        # Is the lesion on senescent tissue
        self.is_senescent = False
        # Number of spores forming the lesion
        self.nb_spores = nb_spores
        # Position of the center of the lesion
        self.position = position
        # List of dispersal units emitted by the lesion
        self.emissions = []
    
    def disable_growth(self):
        """ Shut down lesion growth activity (turn it to False)
        
        Parameters
        ----------
            None
        """
        self.growth_is_active = False
        self.growth_demand = 0.
    
    def disable(self):
        """ Disable all activities of the lesion.
        
        Set the activity of the lesion to False and its growth demand to 0.
        
        Parameters
        ----------
            None
        """
        self.is_active = False
        self.growth_demand = 0.
                
class BiotrophicLesion(Lesion):
    """ Define the interface for a biotrophic lesion.
    """
    
    def __init__(self, position=None, nb_spores=None):
        super(BiotrophicLesion, self).__init__(nb_spores=nb_spores, position=position)
        # Position of senescence the time step before :
        # Used to compute the time left for growth before senescence occur
        self.old_position_senescence = None

    def update(self, leaf=None):
        if not self.is_senescent:
            self.update_age(leaf=leaf)
        else: 
            self.update_age_before_senescence(leaf=leaf)
            
    def senescence_response(self):
        self.disable()

class HemibiotrophicLesion(Lesion):
    """ Define the interface for an hemibiotrophic lesion.
    """
    def __init__(self, position=None, nb_spores=None):
        super(BiotrophicLesion, self).__init__(nb_spores=nb_spores, position=position)
        # Position of senescence the time step before :
        # Used to compute the time left for growth before senescence occur
        self.old_position_senescence = None

    def update(self, leaf=None):
        if not self.is_senescent:
            self.update_age(leaf=leaf)
        else: 
            self.update_age_before_senescence(leaf=leaf)
            
    def senescence_response(self):
        self.disable_growth()
        self.update_surface_dead()
        self.complete_update()
        
class NecrotrophicLesion(Lesion):
    """ Define the interface for a necrotrophic lesion.
    """
    def __init__(self, position=None, nb_spores=None):
        super(NecrotrophicLesion, self).__init__(nb_spores=nb_spores, position=position)
    
    def update(self, leaf=None):
        self.update_age(leaf=leaf)
    
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
