# -*- coding: latin1 -*- 

##
##

""" Description Lesion

"""
from random import random, randint
from math import (exp, pi, floor, ceil)

# Dispersal units #########################################################################

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
        self.active = True
    
    def disable(self):
        """ disable a dispersal unit.
        
        """
        self.active = False
    
    def deposited(self):
        """ Change the status of the spore to 'deposited'.
        
        """
        self.status = 'deposited'
        
    def create_lesion(self, leaf):
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
        les = self.fungus(nb_spores = self.nb_spores, position = self.position)
        if not 'lesions' in leaf.properties():
            leaf.lesions = []
        leaf.lesions.append(les)
        self.disable()
  
class SeptoriaDU(DispersalUnit):
    """ Define a dispersal unit specific of septoria.
    
    """
    fungus = None
    def __init__(self, position=None, nb_spores=None, status=None):
        """ Initialize the dispersal unit of septoria.
        
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
        super(SeptoriaDU, self).__init__(position=position, nb_spores=nb_spores, status=status)
        self.cumul_wetness = 0.
    
    def infect(self, dt, leaf, **kwds):
        """ Compute infection by the dispersal unit of Septoria.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)

        Returns
        -------
            None
        """
        leaf_wet = leaf.wetness # (boolean): True if the leaf sector is wet during this time step.
        temp = leaf.temp # (float) : mean temperature on the leaf sector during the time step (in degree).
        healthy_surface = leaf.healthy_surface # (float) : healthy surface (=with no lesion) on the leaf sector during the time step (in cm^2).
        
        if healthy_surface > 0. :
            # TODO : Right way to do this ?
            if self.nb_spores == 0.:
                self.disable()
           
            else:
                if self.status == 'deposited':
                    # TODO: design a new equation : see Magarey (2005)
                    if leaf_wet:
                        self.cumul_wetness += 1
                    elif self.cumul_wetness > 0: 
                        assert not leaf_wet
                        self.cumul_wetness = 0.
                        # TODO : find a way to reduce inoculum if wet then dry. 
                        # Following lines are a hack - NO biological meaning
                        if proba(self.fungus.loss_rate):
                            self.disable()
                    else:
                        assert not leaf_wet
                        assert self.cumul_wetness == 0.
                    
                    if (self.fungus.temp_min <= temp <= self.fungus.temp_max) and self.cumul_wetness >= self.fungus.wd_min :
                        # TODO : create a function of the number of spores            
                        spores_factor = self.nb_spores / self.nb_spores # always equals 1 for now
                        if proba(spores_factor):
                            self.create_lesion(leaf)
                    elif self.cumul_wetness == 0 :
                        # TODO : Proba conditionnelle doit se cumuler.
                        if proba(self.fungus.loss_rate): 
                            self.disable()
        else:
            self.disable()

class PowderyMildewDU(DispersalUnit):
    """ Define a dispersal unit specific of powdery mildew.
    
    """
    fungus = None
    def __init__(self, position=None, nb_spores=None, status=None):
        """ Initialize the dispersal unit of powdery mildew.
        
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
        super(SeptoriaDU, self).__init__( position=position, nb_spores=nb_spores, status=status)
        
    def infect(self, dt, leaf, **kwds):
        """ Compute infection by the dispersal unit of powdery mildew.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)

        Returns
        -------
            None  
        """
        # External variables
        leaf_wet = leaf.wetness
        temp = leaf.temp
        relative_humidity = leaf.relative_humidity
        
        # Raw parameters for the calculation
        temp_min = self.fungus.temp_min_for_infection
        temp_max = self.fungus.temp_max_for_infection
        m = self.fungus.m_for_infection
        n = self.fungus.n_for_infection
        max_infection_rate = self.fungus.max_infection_rate
        decay_rate = self.fungus.decay_rate
        a_RH_effect = self.fungus.a_RH_effect
        b_RH_effect = self.fungus.b_RH_effect
        RH_opt = self.fungus.RH_opt_for_infection
        c_wetness_effect = self.fungus.c_wetness_effect
        d_wetness_effect = self.fungus.d_wetness_effect
       
        # TODO : Right way to do this ?
        if self.nb_spores == 0.:
            self.disable()
        
        else:
            if self.status == 'deposited':
                # Temperature factor
                temp_norm_function_for_infection = temp_norm_function(temp, temp_min, temp_max, m, n)
                temp_factor = max_infection_rate * temp_norm_function_for_infection * exp(-decay_rate * leaf.age)
                
                # Relative humidity factor
                RH_factor = min(1., a_RH_effect * relative_humidity + b_RH_effect)
                
                # Wetness factor
                if relative_humidity >= RH_opt or leaf_wet:
                    wetness_factor = min(1., c_wetness_effect - d_wetness_effect * temp)
                else:
                    wetness_factor = 1.
                
                # Spores factor
                # TODO : create a function of the number of spores            
                spores_factor = nb_spores / nb_spores # always equals 1 for now
                
                # Infection rate 
                infection_rate = temp_factor * RH_factor * wetness_factor * spores_factor
                if proba(infection_rate):
                    self.create_lesion(leaf)
                else:
                    self.disable()
                    # TODO : review because too extreme.        
       
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
        
        Returns
        -------
            None
        """
        self.active = True
        self.nb_spores = nb_spores
        self.position = position
        self.rings = []
        self.stock_du = []
        self.emissions = []

class Septoria(Lesion):
    """ Define a lesion specific of septoria.
    """
    
    def __init__(self, nb_spores=None, position=None):
        """ Initialize the lesion of septoria. 
        
        Parameters
        ----------
        position: non defined
            Position of the dispersal unit on the phyto-element
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        
        Returns
        -------
            None

        """
        super(Septoria, self).__init__(nb_spores=nb_spores, position=position)
        ring = SeptoriaRing(lesion = self, status = self.fungus.IN_FORMATION, dt=1.)
        self.rings.append(ring)
        self.surface_dead = 0.
    
    def update(self, dt, leaf, **kwds):
        """ Update the status of the lesion and create a new growth ring if needed.
                
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)

        Returns
        -------
            None
        """
        assert self.active

        # remove non active rings
        self.surface_dead += sum(ring.surface for ring in self.rings if not ring.active)
        self.rings = [ring for ring in self.rings if ring.active] 

        if not self.rings:
            self.disable()
            return

        for ring in self.rings:
            ring.update(dt=dt, leaf=leaf, lesion=self, **kwds)

        # Manage the ring in formation and create a new ring when needed
        if len(self.rings) == 1 and self.surface_dead == 0.:
            Dt = self.fungus.degree_days_to_chlorosis
        else:
            Dt = self.fungus.Dt

        if self.can_form_new_ring(Dt):
            # if self.rings[-1].surface > 10.**-6:
            remaining_age = self.rings[-1].age_dday - Dt
            self.rings[-1].status += 1
            new_ring = SeptoriaRing(lesion = self, status = self.fungus.IN_FORMATION, dt=1.)
            self.rings.append(new_ring)
            self.rings[-1].age_dday += remaining_age
            # else:
                # print('desactivation!!')
                # self.rings[-1].disable()

        # Compute emissions
        self.emission(leaf)

    def can_form_new_ring(self, Dt):
        """ Check if the lesion can form a new ring.
        
        A new ring is created when the physiologic age is reached.
        
        Parameters
        ----------
        Dt : int
            Time step in degree days to create a new ring
        
        Returns
        -------
            None
        """
        
        if self.surface < self.fungus.Smax:
            if self.rings[-1].age_dday < Dt :
                return False
            else:
                return True
        else:
            return False
        
        # Keep it in mind --> To be integrated later
        # lesions = leaf.lesions
        # healthy_surface = leaf.healthy_surface
        # incubating_surface = sum([les.surface for les in lesions if les.status == self.INCUBATING])
        # green_surface = healthy_surface + incubating_surface     

        # if( leaf.temp < self.fungus.basis_for_dday or
            # (self.status == self.fungus.IN_FORMATION and self.age_dday == 0.) or
            # green_surface == 0. ):
            # return False
        # else:
            # return True

        #return ring.can_form_new_ring(leaf, lesion)
        
    def emission(self, leaf, **kwds):
        """ Create a list of dispersal units emitted by the entire lesion.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)

        Returns
        -------
            None
        """
        for ring in self.rings:
            emissions = ring.emission(leaf, lesion=self)
            if emissions:
                self.emissions += emissions
    
    def growth_control(self, growth_offer = 0.):
        """ Reduce surface of the last ring up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Surface available on the leaf for the ring to grow (cm2)
        
        Returns
        -------
            None
        """
        if self.rings:
            self.rings[-1].growth_control(growth_offer)
    
    def disable(self):
        """ Set the activity of the lesion to False.
        
        """
        self.active = False
        self.growth_demand = 0.
    
    def is_dead(self):
        """ Update the status of all the rings to 'DEAD' if the lesion is dead.

        """
        return all(ring.is_dead() for ring in self.rings)
  
    @property
    def surface(self):
        """ Compute the surface of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        surface: float
            Surface of the whole lesion (cm2)
        """
        surf = self.surface_dead + sum(ring.surface for ring in self.rings)
        return surf
    
    @property
    def growth_demand(self):
        """ Compute the growth demand of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        growth_demand: float
            Potential growth of the lesion if surface non limiting
        """
        if self.rings:
            return self.rings[-1].growth_demand
    
    @property
    def status(self):
        """ Compute the status of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        status: int
            Status of the lesion
        """
        if self.rings:
            return self.rings[0].status
    
    @property
    def age(self):
        """ Compute the age of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        age_dday: float
            Age of the lesion in hours
        """
        if self.rings:
            return self.rings[0].age
            
    @property
    def age_dday(self):
        """ Compute the thermal age of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        age_dday: float
            Age of the lesion in degree days
        """
        if self.rings:
            return self.rings[0].age_dday

    @status.setter
    def status(self, value):
        """ Set the status of the lesion to the chosen value.
        
        Parameters
        ----------
        value : int
            Chosen value to set lesion status
            
        Returns
        -------
            None
        """
        if self.rings:
            self.rings[0].status = value


class PowderyMildew(Lesion):
    """ 
    """
    
    def __init__(self, nb_spores=None):
        """ Initialize the lesion. 
        
        Parameters
        ----------
        nb_spores: int
            Number of spores into the DU
            
        Returns
        -------
            None
        """
        super(PowderyMildew, self).__init__(nb_spores=nb_spores)
        self.status = self.fungus.LATENT
        ring = PowderyMildewRing(dt=1, lesion=self, status=self.status)
        self.rings.append(ring)
        self.cumul_wetness = 0.
        self.rings.append(ring)
    
    def update(self, dt, leaf, **kwds):
        """ Update the status of the lesion and create a new growth ring if needed.
                
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        
        Returns
        -------
            None
        """        
        if self.status < self.fungus.DEAD:
            # Update the status of each ring
            for ring in self.rings:
                ring.update(dt, leaf, self, **kwds)
                    
    def is_dead(self):
        """ Update the status of all the rings to 'DEAD' if the lesion is dead.
        
        """
        return all(ring.is_dead() for ring in self.rings)
  
    @property
    def surface(self):
        """ Compute the surface of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        surface: float
            Surface of the whole lesion (cm2)
        """
        surf = sum(ring.surface for ring in self.rings)
        return surf

    @property
    def status(self):
        """ Compute the status of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        status: int
            Status of the lesion
        """
        if self.rings:
            return self.rings[0].status
    
    @property
    def age(self):
        """ Compute the age of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        age_dday: float
            Age of the lesion in hours
        """
        if self.rings:
            return self.rings[0].age
            
    @property
    def age_dday(self):
        """ Compute the thermal age of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        age_dday: float
            Age of the lesion in degree days
        """
        if self.rings:
            return self.rings[0].age_dday

    @status.setter
    def status(self, value):
        """ Set the status of the lesion to the chosen value.
        
        Parameters
        ----------
        value : int
            Chosen value to set lesion status
            
        Returns
        -------
            None
        """
        if self.rings:
            self.rings[0].status = value
                      
 ####################################################################################################
def proba(p):
    """ Compute the occurence of an event according to p.

    Parameters
    ----------
    p : float
        Probability of the event in [0,1]
    
    Returns
    -------
    True or False
    """
    return random() < p

class Ring(object):
    """ Ring of Lesion at a given age.
    """

class SeptoriaRing(Ring):
    """ Ring of Lesion of Septoria at a given age.
    """

    def __init__(self, lesion, status, dt=1.):
        """ Initialize each new ring of septoria. 
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        status : int
            Status of the ring when initiated
            
        Returns
        -------
            None 
        """
        super(SeptoriaRing, self).__init__()
        self.status = status
        self.active = True
        self.growth_demand = 0.
        
        if not lesion.rings:
            self.surface = lesion.fungus.epsilon # TOPO : see later with growth_demand
        else:
            self.surface = 0.
        
        self.age = 0.
        self.age_dday = 0.
        self.cumul_rain_event = 0.
        self.rain_before = False
        self.stock_du = []
        

    def is_in_formation(self, fungus):
        """ Can keep growing!!! """
        ok = self.status in (fungus.INCUBATING, fungus.IN_FORMATION)
        if not ok:
            ok = self.is_chlorotic(fungus) and self.age_dday < fungus.Dt
        return ok

    def is_incubating(self, fungus):
        return self.status == fungus.INCUBATING

    def is_chlorotic(self, fungus):
        return self.status == fungus.CHLOROTIC

    def is_necrotic(self, fungus):
        return self.status == fungus.NECROTIC

    def is_sporulating(self, fungus):
        return self.status == fungus.SPORULATING

    def is_empty(self, fungus):
        return self.status == fungus.EMPTY

    def is_dead(self, fungus):
        return self.status == fungus.DEAD

    def update(self, dt, leaf, lesion=None):
        """ Update the status of the ring.
        
        * Cumulate the age of the ring.
        * Assign leaf data to the ring in order to access it in the methods.
        * Call the property 'stage' depending on the current status of the ring.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
            
        Returns
        -------
            None
          
        """
        self.leaf = leaf
        
        # MANAGE TIME SCALE
        self.age += dt
        ddday = (leaf.temp - lesion.fungus.basis_for_dday)/(24./dt)
        self.age_dday += ddday
        # TODO : What if dt > 24 ?
        
        self.stage(dt=dt, ddday=ddday, lesion=lesion)
          
    def in_formation(self, ddday, lesion=None, **kwds):
        """ Cumulate surface while Dt in the state is not reached.        
        
        Parameters
        ----------
        ddday: float
            Number of degree days in 'dt'
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
            
        Returns
        -------
            None
        """     
        fungus = lesion.fungus
        leaf = self.leaf

        assert self.is_in_formation(fungus)

        if self.is_incubating(fungus): # First ring appears in incubation
            self.growth_demand = fungus.Smin * ddday / fungus.degree_days_to_chlorosis
        else:
            size_before_Smax = fungus.Smax - lesion.surface
            self.growth_demand = min(size_before_Smax, fungus.growth_rate * ddday)
        
        #self.surface += self.growth_demand
    
    def growth_control(self, growth_offer = 0.):
        """ Reduce surface of the last ring up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Surface available on the leaf for the ring to grow (cm2)
        
        Returns
        -------
            None
        """
        self.surface += growth_offer
        hs = self.leaf.healthy_surface
        self.leaf.healthy_surface = hs-growth_offer if hs > growth_offer else 0. 

    def chlorotic(self, lesion=None, **kwds):
        """ Set the status of the ring to NECROTIC when needed.
        
        Each ring entering in the CHLOROTIC stage must wait 110 DD 
        to be NECROTIC.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        
        Returns
        -------
            None
        """
        fungus = lesion.fungus
        assert self.is_chlorotic(fungus)
                
        if self.age_dday >= (fungus.degree_days_to_chlorosis +
                             fungus.degree_days_to_necrosis):
            self.status = fungus.NECROTIC 
            
    def necrotic(self, lesion=None, **kwds):
        """ Set the status of the ring to SPORULATING when needed.
        
        Each ring entering in the CHLOROTIC stage must wait ??? DD 
        to be SPORULATING.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        
        Returns
        -------
            None
        """
        fungus = lesion.fungus

        assert self.is_necrotic(fungus)

        if self.age_dday >= (fungus.degree_days_to_chlorosis + 
                             fungus.degree_days_to_necrosis + 
                             fungus.degree_days_to_sporulation):
            self.status = fungus.SPORULATING # The ring begins the sporulation.

    def sporulating(self, lesion=None, **kwds):
        """ Compute the number of rain events on the ring, 
        and update the status of the ring when needed.
        
        A sporulating ring bears fructifications containing dispersal units.
        These dispersal units are spread by the rain if the relative humidity
        is greater than or equal to 85%. It is assumed that the ring is EMPTY
        after 3 separate rain events.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        
        Returns
        -------
            None

        .. Todo:: Enhance the code to parametrize the code with a function.
        """
        fungus = lesion.fungus
        leaf = self.leaf

        assert self.is_sporulating(fungus)
        
        # (float) : rain intensity on the leaf sector during the time step (in mm/h).
        rain_intensity = leaf.rain_intensity 
        # (float) : relative humidity on the leaf sector during the time step (in %).
        relative_humidity = leaf.relative_humidity 
        
        if (rain_intensity > 0. and 
            relative_humidity >= fungus.rh_min and 
            (not self.rain_before)):
            self.cumul_rain_event += 1
        
        if self.cumul_rain_event >= fungus.rain_events_to_empty:
            self.empty(lesion) # The ring is empty.
        
        # Fill the stock of dispersal units according to production rate and surface of the ring
        if not self.stock_du:
            production = self.surface * fungus.production_rate
            nb_du_produced = int(floor(production))
            if proba(1 - (ceil(production) - production)):
                nb_du_produced += 1
            
            self.stock_du = [SeptoriaDU(nb_spores = randint(1,100), status='emitted')
                                for i in range(nb_du_produced)]
                
    def emission(self, leaf, lesion=None, **kwds):
        """ Create a list of dispersal units emitted by the ring.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        
        Returns
        -------
        emissions: list of DUs
            DUs emitted by the ring
        
        .. Todo:: Implement a real formalism.
        """
        emissions = []
        if self.stock_du and self.status == lesion.fungus.SPORULATING and leaf.rain_intensity > 0. and leaf.relative_humidity >= lesion.fungus.rh_min and not self.rain_before:
            nb_du = len(self.stock_du)
            nb_emitted = min(nb_du, int(floor((leaf.rain_intensity / 10) * nb_du)))
            for i in range(nb_emitted):
                idx = randint(0, nb_du-1)
                emissions.append(self.stock_du.pop(idx))
                nb_du -= 1
                if nb_du == 0 :
                    self.status == lesion.fungus.EMPTY

        if leaf.rain_intensity == 0:
            self.rain_before = False
        else:
            self.rain_before = True
        
        return emissions
        
    def empty(self, lesion=None, **kwds):
        """ Disable the 'EMPTY' ring.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        
        Returns
        -------
            None
        """
        self.status = lesion.fungus.EMPTY
        self.disable()
        
    def dead(self, lesion=None, **kwds):
        """ Disable the 'DEAD' ring.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
            
        Returns
        -------
            None
        """
        self.status = lesion.fungus.DEAD
        self.disable()
    
    def disable(self, **kwds):
        """ Set the activity of the ring to 'False'.
        
        Parameters
        ----------
            None
            
        Returns
        -------
            None
        """
        self.active = False
        
    # @property
    def stage(self, dt, ddday, lesion=None):
        """ Orient the ring to toward the proper function according to its status.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        ddday: int 
            Number of degree days in 'dt'
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
            
        Returns
        -------
            Call the method corresponding to ring status  
        """
        if self.status == lesion.fungus.IN_FORMATION:
            return self.in_formation(ddday=ddday, lesion=lesion)
        elif self.status == lesion.fungus.CHLOROTIC:
            return self.chlorotic(lesion=lesion)
        elif self.status == lesion.fungus.NECROTIC:
            return self.necrotic(lesion=lesion)
        elif self.status == lesion.fungus.SPORULATING:
            return self.sporulating(lesion=lesion)
        elif self.status == lesion.fungus.EMPTY:
            return self.empty(lesion=lesion)
        else:
            return self.dead(lesion=lesion)
       
class PowderyMildewRing(Ring):
    """ Ring of Lesion of PowderyMildew at a given age.
    
    """
    
    def __init__(self, status, dt=1.):
        """ Initialize each new ring. 
        
        Parameters
        ----------
        status: int
            Status of the lesion carrying the ring
        
        Returns
        -------
            None
        """
        super(PowderyMildewRing, self).__init__()
        self.status = status
        self.surface = 0.     
        self.age = 0.
        self.age_dday = 0.
        self.latency_progress = 0.
        self.sporulation_progress = 0.
        self.nb_du_dispersed = 0.
        self.nb_du_on_ring -= 0.
    
    def update(self, dt, leaf, lesion=None):
        """ Update the surface and the status of the ring.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
            
        Returns
        -------
            None
        """
        self.leaf = leaf
        
        # TODO : MANAGE TIME SCALE
        self.age += dt
        ddday = (leaf.temp - lesion.fungus.basis_for_dday)/(24/dt)
        self.age_dday += ddday
        # TODO : What if dt > 24 ?
                   
        # Update the surface of the entire lesion
        self.growth(dt, leaf, lesion, **kwds)
        # Update the status of the entire lesion
        self.stage(dt=dt, ddday=ddday, lesion=lesion)
    
    def growth(self, dt, leaf, lesion=None, **kwds):
        """ Compute growth of the ring.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)

        Returns
        -------
            None
        """
        
        # External variables
        temp = leaf.temp
        healthy_surface = leaf.healthy_surface
        lesions = leaf.lesions
        nb_lesions = len(lesions)
        
        # Raw parameters for the calculation
        diameter_max = lesion.fungus.diameter_max
        diameter_min = lesion.fungus.diameter_min
        leaf_age_effect = lesion.fungus.leaf_age_effect
        temp_min = lesion.fungus.temp_min_for_growth
        temp_max = lesion.fungus.temp_max_for_growth
        m = lesion.fungus.m_for_growth
        n = lesion.fungus.n_for_growth
        growth_rate = lesion.fungus.growth_rate
        half_growth_time = lesion.fungus.half_growth_time

        # Maximum diameter of the lesion according to leaf age
        kmax = diameter_min + (diameter_max - diameter_min) * exp(-leaf_age_effect * leaf.age)
        
        # Ring surface according to temperature
        free_space = healthy_surface / nb_lesions
        temp_norm_function_for_growth = temp_norm_function(temp, temp_min, temp_max, m, n)
        self.surface += min(free_space, 
                           (pi/4) * (kmax * temp_norm_function_for_growth * 
                           growth_rate * exp(growth_rate * (half_growth_time - self.age)) /
                           ((1. + exp(growth_rate * (half_growth_time - self.age)))**2))**2)
    
    def latent(self, leaf, **kwds):
        """ Set the status of the lesion to SPORULATING when needed.
               
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
            
        Returns
        -------
            None
        """
        assert(self.status == lesion.fungus.LATENT)
        
        # External variables
        leaf = self.leaf
        temp = leaf.temp
       
        # Raw parameters for the calculation
        temp_min = lesion.fungus.temp_min_for_latency
        temp_max = lesion.fungus.temp_max_for_latency
        m = lesion.fungus.m_for_latency
        n = lesion.fungus.n_for_latency
        min_latency_duration = lesion.fungus.min_latency_duration
        
        # Latency progress
        temp_norm_function_for_latency = temp_norm_function(temp, temp_min, temp_max, m, n)
        self.latency_progress += temp_norm_function_for_latency / min_latency_duration
        
        # Status update
        if self.latency_progress >= 1.:
            self.status = self.fungus.SPORULATING
     
    def sporulating(self, leaf, **kwds):
        """ Compute the number of dispersal events on the lesion, 
        and update the status of the lesion to EMPTY when needed.
               
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
            
        Returns
        -------
            None 
        """
        assert(self.status == lesion.fungus.SPORULATING)
        
        # External variables
        dispersal_rate = random() # TODO : Ajouter lien vers module de dispersion !!!!!!
        temp = leaf.temp
        
        # Raw parameters for the calculation
        a_for_sporulation = lesion.fungus.a_for_sporulation # TODO : diviser par 24 par rapport au Matlab
        b_for_sporulation = lesion.fungus.b_for_sporulation
        temp_min = lesion.fungus.temp_min_for_sporulation
        temp_max = lesion.fungus.temp_max_for_sporulation
        m = lesion.fungus.m_for_sporulation
        n = lesion.fungus.n_for_sporulation
        nb_du_max = lesion.fungus.nb_du_max
        treshold_nb_du_to_empty = lesion.fungus.treshold_nb_du_to_empty
        
        # Number of passed dispersal events
        if self.nb_du_dispersed > 0.:
            self.cumul_dispersal_events += 1
            
        # Progress of the sporulating period
        if temp_max > temp > temp_min:
            self.sporulation_progress += a_for_sporulation * exp(b_for_sporulation * temp)
        
        # Number of dispersal units produced by surface unit by time step
        if self.sporulation_progress < 1.:
            temp_norm_function_for_sporulation = temp_norm_function(temp, temp_min, temp_max, m, n)
            self.nb_du_on_ring += temp_norm_function_for_sporulation * nb_du_max * self.surface
        
        self.nb_du_dispersed = self.nb_du_on_ring * dispersal_rate
        self.nb_du_on_ring -= self.nb_du_dispersed

        # Status update
        if self.nb_du_on_ring < treshold_nb_du_to_empty:
            self.status = lesion.fungus.EMPTY
        
        # TODO : Calculs a revoir... + Integrer stochasticite.
        
    def empty(self, *args, **kwds):
        """ Assert if the status of the ring is 'EMPTY'.
        
        """
        assert(self.status == lesion.fungus.EMPTY)
    
    def dead(self, *args, **kwds):
        """ Assert if the status of the ring is 'DEAD'.

        """
        assert(self.status == lesion.fungus.DEAD)
     
    @property
    def stage(self):
        if self.status == lesion.fungus.LATENT:
            return self.latent
        elif self.status == lesion.fungus.SPORULATING:
            return self.sporulating
        elif self.status == lesion.fungus.EMPTY:
            return self.empty
        else:
            return
    
    pass

def temp_norm_function(temp_mean, temp_min, temp_max, m, n):
    """ Compute the normalized temperature function.
    
    Compute the normalized temperature function as in Calonnec et al., 2008
    
    Parameters
    ----------
    temp_mean: float
        Mean temperature during the time step (°C)
    temp_min: float
        Minimum temperature during the time step (°C)
    temp_max: float
        Maximum temperature during the time step (°C)
    m: float
        Shape parameter
    n: float
        Shape parameter
    
    Returns
    -------
    temp_norm_function: float
        Normalized temperature function
    """
    # Calculation of the normalized temperature
    temp_norm =  min(max(0.,(temp_mean-temp_min)/(temp_max-temp_min)), 1.);
    
    # Calculation of the normalized temperature function
    temp_norm_function = (temp_norm**n)*((1-temp_norm)**m)*((m+n)**(m+n))/(n**n*m**m);
    return temp_norm_function
    
##############################################################################
# Fungus Parameters (e.g. .ini): config of the fungus

class Parameters(object):
    def write(self, filename):
        pass
    
    @staticmethod
    def read(filename):
        pass     
    
class SeptoriaParameters(Parameters):
    def __init__(self,
                 IN_FORMATION = 0,
                 CHLOROTIC = 1,
                 NECROTIC = 2,
                 SPORULATING = 3,
                 EMPTY = 4,
                 DEAD = 5,
                 INCUBATING = 6,
                 Dt = 20,
                 basis_for_dday = -2.,
                 temp_min = 10.,
                 temp_max = 30.,
                 wd_min = 10.,
                 loss_rate = 1./120,
                 degree_days_to_chlorosis = 220,
                 degree_days_to_necrosis = 110,
                 degree_days_to_sporulation = 20,
                 epsilon = 0.001,
                 Smin = 0.03,
                 Smax = 0.3,
                 growth_rate = 0.0006,
                 rh_min = 85.,
                 rain_events_to_empty = 3,
                 production_rate = 10*0.36*6.19*10**3,
                 *args, **kwds):
        """ Parameters for septoria.
        
        Parameters
        ----------
        Dt: int
            Time step in degree days to create a new ring
        basis_for_dday: float
            Basis temperature for the accumulation of degree days (°C)
        temp_min: float
            Minimal temperature for infection
        temp_max: float
            Maximal temperature for infection
        wd_min: float
            Minimal wetness duration for infection
        loss_rate: float
            Loss rate of dispersal units in 1 hour
        degree_days_to_chlorosis: float
            Thermal time between emergence and chlorosis
            (i.e. incubation for the first rings)
        degree_days_to_necrosis: float
            Thermal time between chlorosis and necrosis
            (i.e. incubation for the first rings)
        degree_days_to_sporulation: float
            Thermal time between necrosis and sporulation
        epsilon: float
            Initial size of incubating lesion (cm2)
        Smin: float
            Initial size of chlorotic lesion (cm2)
        Smax: float
            Lesion maximum size (cm2)
        growth_rate: float
            Lesion growth rate (cm2.dday-1)
        rh_min: float
            Minimal relative humidity for sporulation
        rain_events_to_empty: int
            Number of rain events to empty a sporulating ring
        """
        self.name = "Septoria"
        self.IN_FORMATION = IN_FORMATION
        self.CHLOROTIC = CHLOROTIC
        self.NECROTIC = NECROTIC
        self.SPORULATING = SPORULATING
        self.EMPTY = EMPTY
        self.DEAD = DEAD
        self.INCUBATING = INCUBATING
        self.Dt = Dt
        self.basis_for_dday = basis_for_dday
        self.temp_min = temp_min
        self.temp_max = temp_max
        self.wd_min = wd_min
        self.loss_rate = loss_rate
        self.degree_days_to_chlorosis = degree_days_to_chlorosis
        self.degree_days_to_necrosis = degree_days_to_necrosis
        # TODO : Find value for parameters !!
        self.degree_days_to_sporulation = degree_days_to_sporulation
        self.epsilon = epsilon
        self.Smin = Smin
        self.Smax = Smax
        self.growth_rate = growth_rate
        self.rh_min = rh_min
        self.rain_events_to_empty = rain_events_to_empty
        self.production_rate = production_rate
        # TODO : Improve this parameter. 
        # Rapilly would say : RI = 10 mm/h * fDU = 0.36 * pDr = 6.19e3 spores produced by cm2.

    def __call__(self,nb_spores = None, position = None):
        if Septoria.fungus is None:
            Septoria.fungus = self
        if SeptoriaDU.fungus is None:
            SeptoriaDU.fungus = self
        return Septoria(nb_spores = None, position = None)

def septoria(**kwds):
    return SeptoriaParameters(**kwds)
    
class PowderyMildewParameters(Parameters):
    def __init__(self,
                 LATENT = 0,
                 SPORULATING = 1,
                 EMPTY = 2,
                 DEAD = 3,
                 temp_min_for_infection = 5.,
                 temp_max_for_infection = 33.,
                 m_for_infection = 0.338,
                 n_for_infection = 1.055,
                 max_infection_rate = 0.53,
                 decay_rate = 0.147/24,
                 a_RH_effect = 0.0023,
                 b_RH_effect = 0.8068,
                 RH_opt_for_infection = 85,
                 c_wetness_effect = 1.155,
                 d_wetness_effect = 0.014,
                 diameter_max = 5.,
                 diameter_min = 0.2, 
                 leaf_age_effect = 0.08,
                 temp_min_for_growth = 5.,
                 temp_max_for_growth = 33.,
                 m_for_growth = 0.27,
                 n_for_growth = 1.24,
                 growth_rate = 0.2 / 24,
                 half_growth_time = 13. * 24, 
                 temp_min_for_latency = 5.,
                 temp_max_for_latency = 33.,
                 m_for_latency = 0.27,
                 n_for_latency = 1.24,
                 min_latency_duration = 6. * 24,
                 
                 a_for_sporulation = 0.0227 / 24,
                 b_for_sporulation = 0.0762,
                 temp_min_for_sporulation = 15.,
                 temp_max_for_sporulation = 33.,
                 m_for_sporulation = 1.,
                 n_for_sporulation = 1.,
                 nb_du_max = 100,
                 treshold_nb_du_to_empty = 1,
                 
                 *args, **kwds):
        """ Parameters for septoria.
        
        Parameters
        ----------
        temp_min_for_infection: float
            Minimal temperature for infection (°C)
        temp_max_for_infection: float
            Maximal temperature for infection (°C)
        m_for_infection: float
            Shape parameter for the calculation of the normalized temperature for infection
        n_for_infection: float
            Shape parameter for the calculation of the normalized temperature for infection
        max_infection_rate: float
            Maximal infection rate ([0;1])
        decay_rate: float
            Leaf susceptibility decay with age ([0;1])
        a_RH_effect: float
            Shape parameter for the calculation of the effect of relative humidity on infection
        b_RH_effect: float
            Shape parameter for the calculation of the effect of relative humidity on infection
        RH_opt_for_infection: float
            Optimal relative humidity for infection (%)
        c_wetness_effect: float
            Shape parameter for the calculation of the effect of leaf wetness on infection
        d_wetness_effect: float
            Shape parameter for the calculation of the effect of leaf wetness on infection
        diameter_max: float
            Maximum diameter of the lesion (cm)
        diameter_min: float
            Minimum diameter of the lesion (cm)
        leaf_age_effect: float
            Rate of colony growth with the leaves age ([0;1])
        temp_min_for_growth: float
            Minimal temperature for growth (°C)
        temp_max_for_growth: float
            Maximal temperature for growth (°C)
        m_for_growth: float
            Shape parameter for the calculation of the normalized temperature for growth
        n_for_growth: float
            Shape parameter for the calculation of the normalized temperature for growth
        half_growth_time: float
            Time before 50% of colony growth (hours)
        temp_min_for_latency: float
            Minimal temperature for latency (°C)
        temp_max_for_latency: float
            Maximal temperature for latency (°C)
        m_for_latency: float
            Shape parameter for the calculation of the normalized temperature for latency
        n_for_latency: float
            Shape parameter for the calculation of the normalized temperature for latency
        min_latency_duration: float
            Minimum latency duration (hours)
        a_for_sporulation: float
            Shape parameter for the calculation of the progress 
            of the sporulation period according to temperature
        b_for_sporulation: float
            Shape parameter for the calculation of the progress 
            of the sporulation period according to temperature
        temp_min_for_sporulation: float
            Minimal temperature for sporulation (°C)
        temp_max_for_sporulation: float
            Maximal temperature for sporulation (°C)
        m_for_sporulation: float
            Shape parameter for the calculation of the normalized temperature for sporulation
        n_for_sporulation: float
            Shape parameter for the calculation of the normalized temperature for sporulation
        nb_du_max: int
            Number of dispersal unit produced by a sporulating
            surface unit by time step (dispersal unit / cm2 / hour)
        treshold_nb_du_to_empty: int
            Treshold of dispersal unit on a surface to consider it empty
        """
        self.name = "PowderyMildew"
        self.LATENT = LATENT
        self.SPORULATING = SPORULATING
        self.EMPTY = EMPTY
        self.DEAD = DEAD
        self.temp_min_for_infection = temp_min_for_infection
        self.temp_max_for_infection = temp_max_for_infection
        self.m_for_infection = m_for_infection
        self.n_for_infection = n_for_infection
        self.max_infection_rate = max_infection_rate
        self.decay_rate = decay_rate
        self.a_RH_effect = a_RH_effect
        self.b_RH_effect = b_RH_effect
        self.RH_opt_for_infection = RH_opt_for_infection
        self.c_wetness_effect = c_wetness_effect
        self.d_wetness_effect = d_wetness_effect
        
        self.diameter_max = diameter_max
        self.diameter_min = diameter_min
        self.leaf_age_effect = leaf_age_effect
        self.temp_min_for_growth = temp_min_for_growth
        self.temp_max_for_growth = temp_max_for_growth
        self.m_for_growth = m_for_growth
        self.n_for_growth = n_for_growth
        self.growth_rate = growth_rate
        self.half_growth_time = half_growth_time
        
        self.temp_min_for_latency = temp_min_for_latency
        self.temp_max_for_latency = temp_max_for_latency
        self.m_for_latency = m_for_latency
        self.n_for_latency = n_for_latency
        self.min_latency_duration = min_latency_duration
        
        self.a_for_sporulation = a_for_sporulation
        self.b_for_sporulation = b_for_sporulation
        self.temp_min_for_sporulation = temp_min_for_sporulation
        self.temp_max_for_sporulation = temp_max_for_sporulation
        self.m_for_sporulation = m_for_sporulation
        self.n_for_sporulation = n_for_sporulation
        self.nb_du_max = nb_du_max
        self.treshold_nb_du_to_empty = treshold_nb_du_to_empty

        # TODO
        self.dt = 10
    def __call__(self):
        if PowderyMildew.fungus is None:
            PowderyMildew.fungus = self
        if PowderyMildewDU.fungus is None:
            PowderyMildewDU.fungus = self
        
        return PowderyMildew()

def powdery_mildew(**kwds):
    return PowderyMildewParameters(**kwds)
