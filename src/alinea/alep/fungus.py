# -*- coding: latin1 -*-
import copy
##
##
from alinea.alep.alep_objects import FungalObject, DotDict

""" Define the interface for the abstract classes of dispersal units and lesions.
"""


# Dispersal unit ###################################################################################
class DispersalUnit(FungalObject):
    """ Generic class for a dispersal unit (DU) (or cohort of DUs).

    Contains the methods common to all dispersal units in the framework.

    To specialise a dispersal unit for a specific fungus, you can provide a infection model
    """

    def setup(self):
        """initialise the infection model"""
        pass

    def update(self, environment=None):
        """callled evr stop in the tyime sequence, alllow to update with meteo"""
        pass

    def infection_response(self):
        """Called every infection event. should return the number of lesion to be created" and update number of bjects if needed"""
        mutations = {}
        new = self.n_objects
        self.n_objects -= new
        return new, mutations


class Lesion(FungalObject):
    """ Generic class for a lesion (or cohort of lesion).

    Contains the methods common to all dispersal units in the framework.

    To implement a lesion for a specific fungus, you can override the following methods:
        - update()
        - control_growth()
        - emission()
        - senescence_response()
    """

    growth_demand = None

    def setup(self):
        """initialise the infection model"""

    def update(self, environment=None):
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

    def growth_response(self, **ressources):
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
        self.growth_demand = None
        pass

    def emission_response(self, emission_demand, **ressources):
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
        mutations = {}
        return emission_demand, mutations

    def is_sporulating(self):
        """Inform the dispersal model about the current ability of the lesion to emit dispersal units """
        return True

    def is_growing(self):
        return True


# Composition to define a fungus type ##############################################################
class Fungus(object):
    """ Defines a fungus type by combining a lesion type, a dispersal unit type and specific
        parameters
    """
    def __init__(self, lesion=Lesion, dispersal_unit=DispersalUnit, name='basic', **parameters):
        """ Initialize the fungus.

        :Parameters:
         - 'Lesion' (Lesion class): Lesion class for a specific fungus.
         - 'DispersalUnit' (Dispersal class): Dispersal unit class for a specific fungus.
            operates with strictly individual lesions and DUs.
         - 'length_unit' (float): Unit of conversion for dimensions (uses meters as reference:
            thus 0.01 is centimeter)
        """
        self.Lesion = lesion
        self.DispersalUnit = dispersal_unit
        parameters.update({'name':  name})
        self._parameters = parameters

    def parameters(self):
        return DotDict(self._parameters)

    def dispersal_unit(self, nb_dispersal_units=1, **mutations):
        """ Create a dispersal unit instance of the fungus.

        :Parameters:
         - 'kwds' (dict): keys and values for eventual new parameters
        """
        instance = self.DispersalUnit(kind='DispersalUnit', nb_objects=nb_dispersal_units, **mutations)
        instance.fungus = self
        return instance

    def lesion(self, nb_lesions=1, **mutations):
        """ Create a lesion instance of the fungus.

        :Parameters:
         - 'kwds' (dict): keys and values for eventual new parameters
        """
        instance = self.Lesion(kind='Lesion', nb_objects=nb_lesions, **mutations)
        instance.fungus = self
        return instance
