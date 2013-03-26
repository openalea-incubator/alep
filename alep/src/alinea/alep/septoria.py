""" Classes of dispersal unit, lesion and ring specific of wheat septoria.

"""
# Imports #########################################################################
from alinea.alep.cycle2 import *
from random import random, randint
from math import floor, ceil

# Dispersal unit ##################################################################
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

# Lesion ##########################################################################
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
            
            remaining_age = self.rings[-1].age_dday - Dt
            self.rings[-1].status += 1
            new_ring = SeptoriaRing(lesion = self, status = self.fungus.IN_FORMATION, dt=1.)
            self.rings.append(new_ring)
            self.rings[-1].age_dday += remaining_age
            
        # Compute emissions
        
        if leaf.rain_intensity:
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
           if ring.is_stock_available(leaf, lesion = self):
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
            self.rings[-1].growth_control(lesion=self, growth_offer=growth_offer)
    
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

# Rings ###########################################################################
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
        self.first_rain_event = False
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
        
        if leaf.rain_intensity > 0.:
            self.first_rain_event = True if not self.first_rain_event else False
        else:
            self.first_rain_event = False
        
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
           
    def growth_control(self, lesion=None, growth_offer = 0.):
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
        # leaf = lesion.position
        # if leaf != None:
            # hs = leaf.healthy_surface
            # leaf.healthy_surface = hs-growth_offer if hs > growth_offer else 0. 
        # FIXME : Should not need the condition above...
        
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
        fungus = lesion.fungus
        emissions = []
            
        nb_du = len(self.stock_du)
        nb_emitted = min(nb_du, int(floor((leaf.rain_intensity / 10) * nb_du)))
        for i in range(nb_emitted):
            idx = randint(0, nb_du-1)
            emissions.append(self.stock_du.pop(idx))
            nb_du -= 1
            if nb_du == 0 :
                self.status == fungus.EMPTY
        
        return emissions
    
    def is_stock_available(self, leaf, lesion=None):
        """ Check if the stock of DU can be emitted.
        
        DU is free for liberation if :
            - there is DU in stock_du
            - ring is sporulating
            - relative humidity is above a treshold
            - it is the first hour of rain
            
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
        True or False :
            Availability of the stock
        
        """
        fungus = lesion.fungus
        
        if (self.stock_du and self.status == fungus.SPORULATING and 
            leaf.relative_humidity >= fungus.rh_min and
            self.first_rain_event):
            return True
        else:
            return False
        
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

# Fungus parameters (e.g. .ini): config of the fungus #############################
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
            Basis temperature for the accumulation of degree days (degrees celsius)
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

    def __call__(self, nb_spores = None, position = None):
        if Septoria.fungus is None:
            Septoria.fungus = self
        if SeptoriaDU.fungus is None:
            SeptoriaDU.fungus = self
        return Septoria(nb_spores=nb_spores, position=position)

def septoria(**kwds):
    return SeptoriaParameters(**kwds)

# Useful functions ################################################################
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