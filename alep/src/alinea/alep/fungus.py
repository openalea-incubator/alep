# -*- coding: latin1 -*- 
import copy
##
##

""" Define the interface for the abstract classes of dispersal units and lesions.

"""

# Dispersal unit ###################################################################################
class DispersalUnit(object):
    """ Generic class for a dispersal unit (DU) (or cohort of DUs).
    
    Contains the methods common to all dispersal units in the framework. 
    
    To implement a dispersal unit for a specific fungus, you can override the following method:
        - infect()       
    """
    def __init__(self, mutable=False):
        """ Initialize the dispersal unit (DU).
        
        :Parameters:
         - 'mutable' (bool) - True if each DU has its own parameters (intra-population variability),
         False if all DU of the same fungus share the same parameters.
         
        :Attributes:
         - 'is_active' (bool) - Activity of the DU (bool): if False, DU is dead.
         - 'status' (str) - 'emitted' or 'deposited': can be used to distinguish DUs to disperse and 
         DUs that are already deposited on leaves.
         - 'nb_dispersal_units' (int) - A DispersalUnit object can represent a cohort of DU that are
         deposited on the leaf at the same date. 
        """
        self.is_active = True 
        self.status = 'emitted'
        self.nb_dispersal_units = 1

        # Capacity to differ from other lesions of same Fungus
        self.mutable = mutable
        if mutable:
            self.fungus = copy.copy(self.__class__.fungus)
        
    def infect(self, dt=1, leaf=None, **kwds):
        """ Compute the success of infection by the DU.
        
        To be overridden specifically by fungus type. By default, create a new lesion.
        
        :Parameters:
         - 'dt' (int): Time step of the simulation (in hours)
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
         - **kwds : optional arguments depending on fungus specifications
        """
        self.create_lesion(leaf)
    
    def create_lesion(self, nb_lesions = 1, leaf=None, **kwds):
        """ Create new lesions of same fungus type as DU, and disable self if no DU left in cohort.
        
        :Parameters:
         - 'nb_lesions' (int): Number of lesions to create from DU (or cohort of DU)
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
         - **kwds : optional arguments depending on fungus specifications
        """
        if nb_lesions>0:
            les = self.fungus.lesion(mutable = self.mutable)
            les.__dict__.update(kwds)
            les.set_position(self.position)
            self.set_nb_dispersal_units(max(0, self.nb_dispersal_units - nb_lesions))
            if leaf is None:
                self.disable()
                return les
            else:
                try:
                    leaf.lesions.append(les)
                except:
                    leaf.lesions = [les]
                if self.nb_dispersal_units == 0.:
                    self.disable()
                    return
        
    def disable(self):
        """ Disable the dispersal unit.
        """
        self.is_active = False

    def set_position(self, position=None):
        """ Set the position of the DU to position given in argument.
        
        :Parameters:
         - 'position' (x,y) - Position of the DU: coordinates on the leaf axis.
        """
        self.position = position

    def set_status(self, status = 'deposited'):
        """ Set the status of the DU to given argument.
        
        :Parameters:
         - 'status' (str) - 'emitted' or 'deposited': can be used to distinguish DUs to disperse and 
         DUs that are already deposited on leaves.
        """
        self.status = status
        
    def set_nb_dispersal_units(self, nb_dispersal_units=1):
        """ Set the number of DUs in cohort to number given in argument.
        
        :Parameters:
         - 'nb_dispersal_units' (int) - Number of dispersal_units in the cohort
        """
        self.nb_dispersal_units = nb_dispersal_units
        
        
# Lesions ##########################################################################################
class Lesion(object):
    """ Generic class for a lesion (or cohort of lesion).
    
    Contains the methods common to all dispersal units in the framework. 
    
    To implement a dispersal unit for a specific fungus, you can override the following methods:
        - update()
        - control_growth()
        - emission()
        - senescence_response()
    """
    fungus = None
    def __init__(self, mutable=False, nb_lesions=1):
        """ Initialize the lesion. 
        
        :Parameters:
         - 'mutable' (bool) - True if each lesion has its own parameters (intra-population 
         variability), False if all lesions of the same fungus share the same parameters.
         - 'nb_lesions' (int) - A Lesion object can represent a cohort of lesions that are
         deposited on the leaf at the same date. 
         
        :Attributes:
         - 'is_active' (bool) - Activity of the lesion (growth and ageing): if False, lesion is dead.
         - 'growth_is_active' (bool) - Growing activity of the lesion: if False, the growth is 
         stopped but the lesion is still alive and can continue the infection cycle.
         - 'is_senescent' (bool) - True if the lesion is on naturally senescent parts of the leaf.
         - 'growth_demand' (float) - Potential surface increase of the lesion during
         the time step to come.  
        """
        self.is_active = True
        self.growth_is_active = True
        self.is_senescent = False
        self.growth_demand = 0.
        self.nb_lesions = nb_lesions

        # Capacity to differ from other lesions of same Fungus
        self.mutable = mutable
        if mutable:
            self.fungus = copy.copy(self.__class__.fungus)
        
    def update(self, dt, leaf, **kwds):
        """ Update the lesion: compute growth demand and ageing during time step.
        
        Lesion growth is computed in two separate times : 
            - the 'update' method calculates a potential growth of the lesion according to climatic 
            conditions
            - this growth is effective or reduced according to the area available on the leaf
            
        ..seealso:: :method: Lesion.control_growth()
        
        To be overridden specifically by fungus type. By default, pass.
        
        :Parameters:
         - 'dt' (float) - length of time step
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
         - **kwds (dict): optional arguments depending on fungus specifications
        """
        pass
        
    def control_growth(self, growth_offer=0.):
        """ Increase lesion surface according to area available on leaf for this specific 
        lesion ('growth_offer')
        
        The growth offer is specific of each Lesion object and represents the area that will be 
        colonized by the lesion during the time step. It is necessarily inferior or equal to growth 
        demand (i.e. potential growth). It is computed by an external model that coordinates all the 
        lesions on the same leaf.
        
        ..seealso:: :method: Lesion.update()
        
        To be overridden specifically by fungus type. By default, pass.
        
        :Parameters:
         - growth_offer (float): area available on leaf for lesion growth.
        """
        pass
    
    def emission(self, **kwds):
        """ Compute number of dispersal units emitted according to climatic conditions
        
        To be overridden specifically by fungus type. By default, return 1 dispersal unit.
        
        :Parameters:
         - **kwds (dict): optional arguments depending on fungus specifications
         
        :Returns:
         - 'nb_dispersal_units' (int): number of dispersal units emitted
        """
        return 1
    
    def senescence_response(self, **kwds):
        """ Modification of lesion variables in response to leaf natural senescence
        
        Not mandatory for every fungus type, but can be overriden. By default, pass.
        
        :Parameters:
         - **kwds (dict): optional arguments depending on fungus specifications
        """
        pass    
        
    def disable_growth(self):
        """ Shut down lesion growth activity (turn it to False)
        
        Set the growth activity of the lesion to False and its growth demand to 0. If lesion surface
        is inferior to computing precision and growth is disabled, then kill the lesion.
        """
        self.growth_is_active = False
        self.growth_demand = 0.
        if round(self.surface, 16)==0.:
            self.is_active = False
    
    def disable(self):
        """ Disable all activities of the lesion.
        
        Set the activity of the lesion to False and its growth demand to 0.
        """
        self.disable_growth()
        self.is_active = False
        
    def set_position(self, position=None):
        """ Set the position of the lesion to position given in argument.
        
        :Parameters:
         - 'position' (x,y) - Position of the lesion: coordinates on the leaf axis.
        """
        self.position = position
            
    def become_senescent(self):
        """ Set the senescence marker to True.
        """
        self.is_senescent = True

# Composition to define a fungus type ##############################################################
class Fungus(object):
    """ Defines a fungus type by combining a lesion type, a dispersal unit type and specific 
        parameters
    """
    def __init__(self, Lesion=Lesion,
                 DispersalUnit = DispersalUnit, 
                 parameters = {'name':'template', 'group_dus':'False'},
                 length_unit = 0.01):
        """ Initialize the fungus.
        
        :Parameters:
         - 'Lesion' (Lesion class): Lesion class for a specific fungus.
         - 'DispersalUnit' (Dispersal class): Dispersal unit class for a specific fungus.
         - 'parameters' (dict): dict of fungus specific parameters, mandatory:
            * 'name' (str): name of fungal species
            * 'group_dus' (bool): True if the model operates with cohorts, False if the model 
            operates with strictly individual lesions and DUs.
         - 'length_unit' (float): Unit of conversion for dimensions (uses meters as reference:
            thus 0.01 is centimeter)
        """
        self.length_unit = length_unit
        self.Lesion_class = Lesion
        self.DispersalUnit_class = DispersalUnit
        self.parameter_names = parameters.keys()
        self.__dict__.update(parameters)

    def update_parameters(self, **kwds):
        """ Get and update fungus parameters with parameters in kwds.
        
        :Parameters:
         - 'kwds' (dict): keys and values for new parameters
        """
        self.__dict__.update(kwds)
        if len(kwds)>0:
            self.parameter_names += [k for k in kwds.iterkeys()]
        
    def parameters(self, **kwds):
        """ Get parameters of the fungus.
        
        :Parameters:
         - 'kwds' (dict): keys and values for eventual new parameters
        """
        self.update_parameters(**kwds)
        return {k:getattr(self,k) for k in self.parameter_names}
    
    def dispersal_unit(self, mutable=False, **kwds):
        """ Create a dispersal unit instance of the fungus.
        
        :Parameters:
         - 'kwds' (dict): keys and values for eventual new parameters
        """
        self.update_parameters(**kwds)
        self.DispersalUnit_class.fungus = self
        instance = self.DispersalUnit_class(mutable=mutable)
        return instance

    def lesion(self, mutable=False, **kwds):
        """ Create a lesion instance of the fungus.
        
        :Parameters:
         - 'kwds' (dict): keys and values for eventual new parameters
        """
        self.update_parameters(**kwds)
        self.Lesion_class.fungus = self
        instance = self.Lesion_class(mutable=mutable)
        return instance
