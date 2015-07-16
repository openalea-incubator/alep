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
    # fungus = None
    def __init__(self, mutable=False):
        """ Initialize the dispersal unit.
        
        Parameters
        ----------
        
        """
        self.mutable = mutable
        self.is_active = True
        self.status = 'emitted'
        self.nb_dispersal_units = 1
        if mutable:
            self.fungus = copy.copy(self.__class__.fungus)
    
    def disable(self):
        """ disable a dispersal unit.
        
        """
        self.is_active = False
    
    
    def set_position(self, position=None):
        """ Set the position of the DU to position given in argument.
        
        Parameters
        ----------
        position: (x,y)
            Position of the DU.
        """
        self.position = position
        
    def set_nb_dispersal_units(self, nb_dispersal_units=1):
        """ Set the number of DUs in cohort to number given in argument.
        
        Parameters
        ----------
        nb_dispersal_units: int
            Number of dispersal_units in the cohort
        """
        self.nb_dispersal_units = nb_dispersal_units
    
    def infect(self, leaf=None):
        self.create_lesion(leaf)
    
    def create_lesion(self, leaf=None, **kwds):
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
        les = self.fungus.lesion(mutable = self.mutable)
        les.__dict__.update(kwds)
        les.set_position(self.position)
        if leaf is None:
            self.disable()
            return les
        else:
            try:
                leaf.lesions.append(les)
            except:
                leaf.lesions = [les]
            self.disable()
             
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
    
    And the lesion have to a provide the following attributes:
        - is_dead
        - surface
        - status
        - age
        - age_physio
        
    ..todo:: improve header and doc
    """
    fungus = None
    def __init__(self, mutable=False):
        """ Initialize the lesion. 
        
        Parameters
        ----------
        position: not defined
            Position of the dispersal unit on the phyto-element
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        
        """
        self.mutable = mutable
        # Total activity of the lesion (growth and ageing)
        self.is_active = True
        # Growth activity of the lesion
        self.growth_is_active = True
        # Growth demand of the lesion (cm2)
        self.growth_demand = None
        # Is the lesion on senescent tissue
        self.is_senescent = False
        # Position of the lesion
        self.position = None
        # Capacity to differ from other lesions of same Fungus
        if mutable:
            self.fungus = copy.copy(self.__class__.fungus)
        
    def update(self, dt, leaf, growth_rate=1e-4):
        pass
    
    def emission(self):
        pass
    
    def senescence_response(self, **kwds):
        pass
        
    def control_growth(self, growth_offer=0.):
        pass        
        
    def disable_growth(self):
        """ Shut down lesion growth activity (turn it to False)
        
        Parameters
        ----------
            None
        """
        self.growth_is_active = False
        self.growth_demand = 0.
        if round(self.surface, 16)==0.:
            self.disable()
    
    def disable(self):
        """ Disable all activities of the lesion.
        
        Set the activity of the lesion to False and its growth demand to 0.
        
        Parameters
        ----------
            None
        """
        self.growth_is_active = False
        self.is_active = False
        self.growth_demand = 0.
        
    def set_position(self, position):
        self.position = position
    
    def become_senescent(self):
        self.is_senescent = True


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
    
    def __init__(self, name='template', Lesion=Lesion,
                 DispersalUnit = DispersalUnit, 
                 parameters = {'name':'template', 'group_dus':'False'},
                 length_unit = 0.01):
        self.name = name
        self.length_unit = length_unit
        self.Lesion_class = Lesion
        self.DispersalUnit_class = DispersalUnit
        self.parameter_names = parameters.keys()
        self.__dict__.update(parameters)
    
#    @classmethod
    def parameters(self, **kwds):
        self.__dict__.update(kwds)
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