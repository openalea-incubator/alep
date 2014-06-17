# -*- coding: latin1 -*- 
import copy
##
##

""" Define the interface for the fungal classes of dispersal units, lesions,
    rings and parameters.

"""

# Dispersal unit ##################################################################
class DispersalUnit(object):
    """ Generic class for a dispersal unit.
    
    """
    fungus = None
    def __init__(self, nb_spores=None, position=None, mutable=False):
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
        self.is_active = True
        self.can_infect_at_position = None
        if mutable:
            self.fungus = copy.copy(self.__class__.fungus)
    
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
    
    def set_can_infect(self, can_infect=True):
        """ Turn the property 'can_infect_at_position' to True.
        """
        self.can_infect_at_position = can_infect
    
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
        les = self.fungus.lesion()
        if leaf is None:
            self.disable()
            return les
        else:
            try:
                leaf.lesions.append(les)
            except:
                leaf.lesions = [les]
            self.disable()
    
    @property
    def nb_dispersal_units(self):
        if self.position is None:
            return None
        else:
            return len(self.position)
             
# Lesions #################################################################################
# /!\ TODO : Dans le cas ou on veut des parametres variables par lesion, il faut ecraser la variable
# de classe 'fungus' par une variable d'instance. A mettre en option de init de lesion.

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
    def __init__(self, position=None, nb_spores=None, mutable=False):
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
        if mutable:
            self.fungus = copy.copy(self.__class__.fungus)
        
    def update(self, dt, leaf):
        pass
    
    def emission(self, emission_rate = 1e4):
        """ This method simulates the biological regulation of spore emission. Not generic
        
        Parameters to this function generaly are expressed as a climatic-dependant potential rate of emmission
        """
        nb_dus = int(emission_rate * self.fungus.length_unit**2)
        
        return create_dispersal_units(nb_dus)
    
    def senescence_response(self, **kwds):
        pass
        
    def control_growth(self, growth_offer=0.):
        pass
        
    def become_senescent(self):
        """ The lesion will become senescent during this time step.
        """
        self.is_senescent = True    
    
    def create_dispersal_units(self, nb_dus=1):
        """ Generic method to create new dispersal units.
        """
        return [self.fungus.dispersal_unit() for i in range(nb_dus)]
        
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
        
    def set_position(self, position):
        self.position = position
        

# Fungus Parameters (e.g. .ini): config of the fungus ####################################
#class Parameters(object):
#    def write(self, filename):
#        pass
    
#    @staticmethod
#    def read(filename):
#        pass
        

 
# class _Fungus(object):
    # def __init__(self, name='base_fungus', Lesion=Lesion, DispersalUnit = DispersalUnit, parameters = {}):
        # self.name = name
        # self.Lesion = Lesion
        # self.DispersalUnit = DispersalUnit
        # for k,v in parameters.iteritems():
            # exec "self.%s = %s"%(k,v)
            
class Fungus(object):
    
    def __init__(self, name='template', Lesion=Lesion, DispersalUnit = DispersalUnit, parameters = {}, length_unit = 0.01):
        self.name = name
        self.length_unit = length_unit
        self.Lesion_class = Lesion
        self.DispersalUnit_class = DispersalUnit
        self.parameter_names = parameters.keys()
        self.__dict__.update(parameters)
    
#    @classmethod
    def parameters(self):
        return {k:getattr(self,k) for k in self.parameter_names}
    
#    @classmethod
    def dispersal_unit(self, mutable=False, **kwds):
        self.__dict__.update(kwds) #should theses new attributes be added to parameters ?
        self.DispersalUnit_class.fungus = self
        instance = self.DispersalUnit_class(mutable=mutable)
        return instance
        
#    @classmethod
    def lesion(self, mutable=False, **kwds):
        self.__dict__.update(kwds)
        self.Lesion_class.fungus = self
        instance = self.Lesion_class(mutable=mutable)
        return instance

## TEMP ##
class Parameters(object):
    pass