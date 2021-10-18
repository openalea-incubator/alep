# -*- coding: latin1 -*-
import copy
##
##

""" Define the interface for the abstract classes of dispersal units and lesions.
"""


class DotDict(dict):
    """
    a dictionary that supports dot notation
    as well as dictionary access notation
    usage: d = DotDict() or d = DotDict({'val1':'first'})
    set attributes: d.val2 = 'second' or d['val2'] = 'second'
    get attributes: d.val2 or d['val2']
    """
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, dct = None):
        if dct is None:
            dct = {}
        for key, value in dct.items():
            if hasattr(value, 'keys'):
                value = DotDict(value)
            self[key] = value

# Dispersal unit ###################################################################################
class DispersalUnit(object):
    """ Generic class for a dispersal unit (DU) (or cohort of DUs).

    Contains the methods common to all dispersal units in the framework.

    To implement a dispersal unit for a specific fungus, you can override the following method:
        - infect()
    """
    fungus = None

    def __init__(self, nb_dispersal_units=1, **mutations):
        """ Initialize the dispersal unit (DU).

        :Parameters:
         - 'mutable' (bool) - True if each DU has its own parameters (intra-population variability),
         False if all DU of the same fungus share the same parameters.

        :Attributes:
         - 'nb_dispersal_units' (int) - A DispersalUnit object can represent a cohort of DU that are
         deposited on the leaf at the same date.
        """
        self.is_active = True
        self.nb_dispersal_units = nb_dispersal_units
        self.mutations = mutations

    def parameters(self):
        """ Get parameters of parent fungus, updated with mutations
        """
        pars = {}
        if self.fungus:
            pars.update(self.fungus.parameters)
        pars.update(self.mutations)
        return DotDict(pars)

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

    To implement a lesion for a specific fungus, you can override the following methods:
        - update()
        - control_growth()
        - emission()
        - senescence_response()
    """

    # fungus is inherited by Fungus
    fungus = None

    def __init__(self, nb_lesions=1, **mutations):
        """ Initialize a lesion.

        :Parameters:
         - 'nb_lesions' (int) - A Lesion object can represent a cohort of lesions that are
         deposited on the leaf at the same date.

        :Attributes:
         - 'is_active' (bool) - Activity of the lesion (growth and ageing): if False, lesion is dead.
         - 'growth_is_active' (bool) - Growing activity of the lesion: if False, the growth is
         stopped but the lesion is still alive and can continue the infection cycle.
         - 'is_senescent' (bool) - True if the lesion is on naturally senescent parts of the leaf.
         - 'growth_demand' (float) - Potential surface increase of the lesion during
         the time step to come.
         - 'is_sporulating': whether the lesion is emiting spores or not
        """
        self.is_active = True
        self.growth_is_active = True

        self.nb_lesions = nb_lesions

        # to check
        self.is_senescent = False
        self.growth_demand = 0.
        self.position = None

        self.mutations = mutations

    def parameters(self):
        """ Get parameters of parent fungus, updated with mutations
        """
        pars = {}
        if self.fungus:
            pars.update(self.fungus.parameters)
        pars.update(self.mutations)
        return DotDict(pars)

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

# interface with dispersal models
#################################

    def emission(self, emission_demand=1, **mutations):
        """ return dispersal units the lesion can emmit as a function of an emmission demand originated from a dispersal model.
            This method is NOT intended to specify the emission demand, but how the lesion accommodates to such a demand,
            and do required housekeeping

        To be overridden specifically by fungus type.

        :Parameters:
         - emission_demand (int): number of dispersal units to be produced according to dispersal emission model
         - mutations : not used yet

        :Returns:
        - a list of dispersal units. If all dispersal units share the same parameters (no mutations), the list resume
        to a single element, with nb_dispersal_units appropriately set.
        """
        if self.fungus is None:
            raise TypeError("fungus undefined : lesion should be instantiated with fungus method for lesion.emission "
                            "to properly work")
        transmitted_mutations = copy.deepcopy(self.mutations)
        transmitted_mutations.update(mutations)
        return self.fungus.dispersal_unit(emission_demand, **transmitted_mutations)

    def is_sporulating(self):
        """Inform the dispersal model about the current ability of the lesion to emit dispersal units """
        return True

# end dispersal interface
#################################

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
        if round(self.surface, 16) == 0.:
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
    def __init__(self, lesion=Lesion, dispersal_unit=DispersalUnit, name='basic', length_unit=0.01, **parameters):
        """ Initialize the fungus.

        :Parameters:
         - 'Lesion' (Lesion class): Lesion class for a specific fungus.
         - 'DispersalUnit' (Dispersal class): Dispersal unit class for a specific fungus.
            operates with strictly individual lesions and DUs.
         - 'length_unit' (float): Unit of conversion for dimensions (uses meters as reference:
            thus 0.01 is centimeter)
        """
        self.length_unit = length_unit
        self.name = name
        self.Lesion = lesion
        self.DispersalUnit = dispersal_unit
        if parameters is None:
            parameters = {}
        self.parameters = DotDict(parameters)

    def dispersal_unit(self, nb_dispersal_units=1, **mutations):
        """ Create a dispersal unit instance of the fungus.

        :Parameters:
         - 'kwds' (dict): keys and values for eventual new parameters
        """
        self.DispersalUnit.fungus = self
        instance = self.DispersalUnit(nb_dispersal_units=nb_dispersal_units, **mutations)
        return instance

    def lesion(self, nb_lesions=1, **mutations):
        """ Create a lesion instance of the fungus.

        :Parameters:
         - 'kwds' (dict): keys and values for eventual new parameters
        """
        self.Lesion.fungus = self
        instance = self.Lesion(nb_lesions=nb_lesions, **mutations)
        return instance
