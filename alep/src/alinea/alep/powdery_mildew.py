""" Classes of dispersal unit, lesion and ring specific of grapevine powdery mildew.

"""
# Imports #########################################################################
from alinea.alep.cycle2 import *
from random import random, randint
from math import (exp, pi, floor, ceil)

# Dispersal unit ##################################################################
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
        super(PowderyMildewDU, self).__init__( position=position, nb_spores=nb_spores, status=status)
        
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

# Lesion ##########################################################################
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

# Ring ############################################################################
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

# Fungus parameters (e.g. .ini): config of the fungus #############################
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
            Minimal temperature for infection (degrees celsius)
        temp_max_for_infection: float
            Maximal temperature for infection (degrees celsius)
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
            Minimal temperature for growth (degrees celsius)
        temp_max_for_growth: float
            Maximal temperature for growth (degrees celsius)
        m_for_growth: float
            Shape parameter for the calculation of the normalized temperature for growth
        n_for_growth: float
            Shape parameter for the calculation of the normalized temperature for growth
        half_growth_time: float
            Time before 50% of colony growth (hours)
        temp_min_for_latency: float
            Minimal temperature for latency (degrees celsius)
        temp_max_for_latency: float
            Maximal temperature for latency (degrees celsius)
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
            Minimal temperature for sporulation (degrees celsius)
        temp_max_for_sporulation: float
            Maximal temperature for sporulation (degrees celsius)
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

def temp_norm_function(temp_mean, temp_min, temp_max, m, n):
    """ Compute the normalized temperature function.
    
    Compute the normalized temperature function as in Calonnec et al., 2008
    
    Parameters
    ----------
    temp_mean: float
        Mean temperature during the time step (degrees celsius)
    temp_min: float
        Minimum temperature during the time step (degrees celsius)
    temp_max: float
        Maximum temperature during the time step (degrees celsius)
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
    temp_norm =  min(max(0.,(temp_mean-temp_min)/(temp_max-temp_min)), 1.)
    
    # Calculation of the normalized temperature function
    temp_norm_function = (temp_norm**n)*((1-temp_norm)**m)*((m+n)**(m+n))/(n**n*m**m)
    return temp_norm_function