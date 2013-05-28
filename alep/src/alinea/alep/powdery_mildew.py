""" Classes of dispersal unit, lesion of grapevine powdery mildew.

"""
# Imports #########################################################################
from alinea.alep.fungal_objects import *
from random import random, randint
from math import (exp, pi, floor, ceil, sqrt)

# Dispersal unit ##################################################################
class PowderyMildewDU(DispersalUnit):
    """ Define a dispersal unit specific of powdery mildew.
    
    """
    fungus = None
    def __init__(self, nb_spores=None, position=None, status=None):
        """ Initialize the dispersal unit of powdery mildew.
        
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
        super(PowderyMildewDU, self).__init__(nb_spores=nb_spores, position=position, status=status)
        # Age of the dispersal unit in degree days
        self.age_dday = 0.
        # Viability of the dispersal unit:
        # Scale from 1:viable to 0:not viable
        # self.viability = 1.
        
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
        
        # Parameters for calculation
        f = self.fungus
        t_min = f.temp_min_for_infection
        t_max = f.temp_max_for_infection
        m = f.m_for_infection
        n = f.n_for_infection
        max_rate = f.max_infection_rate
        decay_rate = f.decay_rate
        a_RH_effect = f.a_RH_effect
        b_RH_effect = f.b_RH_effect
        RH_opt = f.RH_opt_for_infection
        c_wetness_effect = f.c_wetness_effect
        d_wetness_effect = f.d_wetness_effect
        nb_spores = self.nb_spores
       
        if nb_spores == 0.:
            self.disable()
            return
        
        # Compute dispersal unit age (and viability?)
        self.update_age_dday(dt, leaf)
        
        if self.is_ready_to_infect():
            # Temperature factor
            t_norm_function = temp_norm_function(temp, t_min, t_max, m, n)
            temp_factor = (max_rate * t_norm_function * 
                            exp(-decay_rate * leaf.age))
            
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
                # self.update_viability(dt, leaf)
                # if not self.is_viable():
                    # self.disable()
                
    def update_age_dday(self, dt=1., leaf=None):
        """ Update the age of the lesion.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        self.age_dday += self.compute_delta_ddays(dt, leaf)
        # self.age_dday += self.compute_delta_ddays_from_weather(leaf)
        
    def compute_delta_ddays(self, dt=1., leaf=None):
        """ Compute delta degree days in dt.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        
        Returns
        -------
        ddday: float
            Delta degree days in time step
        """
        f = self.fungus
        # Calculation
        if dt != 0.:
            ddday = max(0,(leaf.temp - f.basis_for_dday)/(24./dt))
        else:
            ddday = 0.
        return ddday 
    
    def compute_delta_ddays_from_weather(self, leaf=None):
        """ Compute delta degree days from weather data since last call.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
            
        Returns
        -------
        ddday: float
            Delta degree days in time step
        """
        f = self.fungus
        temp_list = np.array(leaf.temp_list)
        dt = len(temp_list)
        # Calculation
        if dt != 0.:
            ddday = max(0, sum((temp_list - f.basis_for_dday))/24.)
        else:
            ddday = 0.
        return ddday
    
    def update_viability(self, dt=1., leaf=None):
        """ Update the viability of the dispersal unit.
        
        The hypothesis is made that a dispersal unit is viable 5 days (120h).
        The decrease of viability by time step is doubled if the leaf is wet.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        f = self.fungus
        if leaf.wetness:
            dt *= 2.
        viability_loss = dt/f.viability_length
        self.viability -= viability_loss if self.viability > viability_loss else 0.
    
    def is_viable(self):
        """ Check if the dispersal unit is still viable.
        
        Parameters
        ----------
            None
        
        Returns
        -------
        True or False:
            True if the dispersal unit is viable, False otherwise
        """
        return self.viability > 0.
                    
    def is_ready_to_infect(self):
        """ Check if the dispersal unit is ready to infect according to its age.
        
        Parameters
        ----------
            None
        
        Returns
        -------
        True or False:
            True if the dispersal unit is ready for infection, False otherwise
        """
        f = self.fungus
        return self.age_dday >= f.degree_days_to_infect

# Lesion ##########################################################################
class PowderyMildew(Lesion):
    """ Powdery mildew lesion implemented as in Calonnec et al., 2008 for the most part. """
    
    def __init__(self, nb_spores=None, position=None):
        """ Initialize the lesion of powdery mildew.
        
        Parameters
        ----------
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        position: non defined
            Position of the dispersal unit on the phyto-element
        """
        super(PowderyMildew, self).__init__(nb_spores=nb_spores, position=position)
        f = self.fungus
        # Duration of the time step
        self.dt = 0.
        # Age of the lesion (in hours)
        self.age = 0.
        # Surface of the lesion
        self.surface = 0.
        # Diameter of the lesion
        self.diameter = 0.
        # Growth demand of the lesion
        self.growth_demand = None
        # Temperature around the lesion
        self.temp = None
        # Maximum diameter of the lesion according to its age
        self.diameter_max = f.diameter_max        
        # Status of the lesion
        self.status = f.LATENT
        # Cumulation of wetness periods
        self.cumul_wetness = 0.
        # Progression in latency period
        self.latency_progress = 0.
        # Progression in sporulation period
        self.sporulation_progress = 0.
        # Spore production activity
        self.production_is_active = False
        # Stock of spores
        self.stock_spores = 0.
    
    def update(self, dt=1., leaf=None):
        """ Update the status of the lesion and create a new growth ring if needed.
 
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        """       
        # Update lesion age
        self.dt = dt
        self.age += dt

        # Update temperature around the lesion
        self.temp = leaf.temp

        # Update lesion surface
        if self.growth_is_active:
            self.update_diameter_max(leaf)
            self.update_growth_demand()

        # Update lesion status
        if self.is_latent():
            self.progress_in_latency()
        elif self.is_sporulating():
            self.progress_in_sporulation()
    
    def update_diameter_max(self, leaf=None):
        """ Compute maximum diameter of the lesion according to its age.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        # Parameters for the calculation
        f = self.fungus
        diameter_max = f.diameter_max
        diameter_min = f.diameter_min
        leaf_age_effect = f.leaf_age_effect
        
        # Maximum diameter of the lesion according to leaf age
        kmax = diameter_min + (diameter_max - diameter_min) * exp(-leaf_age_effect * leaf.age)
        self.diameter_max = kmax
        
    def update_growth_demand(self):
        """ Compute lesion growth demand for the time step.
        
        Parameters
        ----------
            None
        """
        # Parameters for the calculation
        f = self.fungus
        temp = self.temp
        dt = self.dt
        diameter = self.diameter
        t_min = f.temp_min_for_growth
        t_max = f.temp_max_for_growth
        m = f.m_for_growth
        n = f.n_for_growth
        r = f.growth_rate
        half_time = f.half_growth_time
        
        # Growth demand in diameter
        kmax = self.diameter_max
        age = self.age
        t_norm_function = temp_norm_function(temp, t_min, t_max, m, n)
        age_factor = r*exp(r*(half_time-age))/(1+exp(r*(half_time-age)))**2
        diameter_demand = diameter + kmax * t_norm_function * age_factor
       
        # Growth demand in surface
        gd = (pi * diameter_demand**2 /4) - self.surface
        self.growth_demand = gd
    
    def progress_in_latency(self):
        """ Compute the ageing of the lesion in latency period and 
            set the status of the lesion to SPORULATING when needed.
               
        Parameters
        ----------
            None
        """
        # Parameters for the calculation
        f = self.fungus
        assert(self.status == f.LATENT)
        temp = self.temp
        t_min = f.temp_min_for_latency
        t_max = f.temp_max_for_latency
        m = f.m_for_latency
        n = f.n_for_latency
        min_latency_duration = f.min_latency_duration
        
        # Latency progress
        t_norm_function = temp_norm_function(temp, t_min, t_max, m, n)
        self.latency_progress += t_norm_function / min_latency_duration
        
        # Status update
        if self.latency_progress >= 1.:
            self.status = f.SPORULATING
            self.activate_production()
    
    def progress_in_sporulation(self):
        """ Compute the ageing of the lesion in latency period and 
            set the status of the lesion to SPORULATING when needed.
               
        Parameters
        ----------
            None
        """
        # Parameters for the calculation
        f = self.fungus
        assert(self.status == f.SPORULATING)
        temp = self.temp
        t_min = f.temp_min_for_sporulation
        t_max = f.temp_max_for_sporulation
        a = f.a_for_sporulation
        b = f.b_for_sporulation
        
        # Progress of the sporulating period
        if t_max > temp > t_min:
            self.sporulation_progress += a * exp(b * temp)
        
        # Production activity update
        if self.sporulation_progress >= 1.:
            self.disable_production()
    
    def control_growth(self, growth_offer=0.):
        """ update surface of the lesion according to its growth demand
            and available surface on the leaf.
        
        Parameters
        ----------
        growth_offer: float
            Minimum between 'growth_demand' and the surface available on
            the leaf for the lesion to grow (cm2)
        """
        # 1. Update surface
        if self.growth_is_active:
            # Growth offer is added to surface alive
            self.surface += growth_offer
            self.diameter = sqrt(self.surface*4/pi)
            # Maximum diameter according to lesion age and leaf age
            kmax = self.diameter_max
            # Check if any interruption of growth:
            if growth_offer < self.growth_demand: # or round(self.diameter,3) >= round(kmax,3):
                # /!\ TODO : Check condition above
                self.disable_growth()

        # 2. Update stock of spores
        if self.production_is_active:
            self.update_stock()
    
    def update_stock(self):
        """ Update the stock of spores available on the lesion
            according only to surface.
        
        Not used in Calonnec version.
        
        Parameters
        ----------
            None
        """
        # Parameters for the calculation
        f = self.fungus
        dt = self.dt
        surface = self.surface
        beta = f.beta_for_sporulation
        gamma = f.gamma_for_sporulation
        
        # Production of spores during dt
        production = beta*exp(gamma*surface)*dt/24
        
        self.stock_spores += production
    
    def update_stock_temperature_dependent(self):
        """ Update the stock of spores available on the lesion 
            according to surface and temperature.
        
        Not used in Calonnec version.
        
        Parameters
        ----------
            None
        """
        # Parameters for the calculation
        f = self.fungus
        dt = self.dt
        temp = self.temp
        surface = self.surface
        beta = self.beta_for_sporulation
        gamma = gamma_for_sporulation
        t_min = f.temp_min_for_sporulation
        t_max = f.temp_max_for_sporulation
        m = f.m_for_sporulation
        n = f.n_for_sporulation
                
        # Production of spores during dt
        max_rate = beta*exp(gamma*surface)*dt
        t_norm_function = temp_norm_function(temp, t_min, t_max, m, n)
        production = max_rate * t_norm_function
        
        self.stock_spores += production
    
    def emission(self, leaf=None):
        """ Create a list of dispersal units emitted by the ring.
        
        Computed with the formalism of Willocquet and Clerjeau, 1998.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        # Parameters for the calculation
        wind_speed = leaf.wind_speed
        f = self.fungus
        stock = self.stock_spores
        a = f.a_for_emission
        b = f.b_for_emission
        r = f.r_for_emission
        
        emissions = []
        if stock > 0.:
            # Dispersal rate
            dispersal_rate = min(1, exp(r*wind_speed+b) / (a*(1+exp(r*wind_speed+b))))
            # Number of spores emitted
            nb_spores = int(dispersal_rate * stock)
            
            if nb_spores > 0 :
                # One only spore by DU
                nb_spores_by_DU = 1
                nb_DU_emitted = nb_spores
                # Return emissions
                PowderyMildewDU.fungus = f
                emissions = [PowderyMildewDU(nb_spores=nb_spores_by_DU, status='emitted')
                                        for i in range(nb_DU_emitted)]
                # Update dtock of spores                       
                self.stock_spores -= nb_spores if self.stock_spores>nb_spores else 0.
        
        # Disable the lesion if not producing spore anymore and stock is empty
        if (self.is_sporulating() and 
            not self.production_is_active and
            self.stock_spores==0.):
            self.disable()

        return emissions
        
    def senescence_response(self):
        """ Disable the lesion if it is on senescent tissue.
        
        Parameters
        ----------
        """
        self.disable()
        
    def is_latent(self):
        """ Check if the status of the lesion is LATENT.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        return self.status==f.LATENT
    
    def is_sporulating(self):
        """ Check if the status of the lesion is SPORULATING.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        return self.status==f.SPORULATING

    def disable_production(self):
        """ Shut down lesion production activity (turn it to False)
        
        Parameters
        ----------
            None
        """
        self.production_is_active = False

    def activate_production(self):
        """ Turn lesion production activity on (turn it to True)
        
        Parameters
        ----------
            None
        """
        self.production_is_active = True
        
# Fungus parameters (e.g. .ini): config of the fungus #############################
class PowderyMildewParameters(Parameters):
    def __init__(self,
                 LATENT = 0,
                 SPORULATING = 1,
                 EMPTY = 2,
                 DEAD = 3,
                 degree_days_to_infect = 20.,
                 basis_for_dday = 0.,
                 viability_length = 120.,
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
                 diameter_max = 1.8,
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
                 beta_for_sporulation = 874.,
                 gamma_for_sporulation = 0.314,
                 
                 temp_min_for_sporulation = 15.,
                 temp_max_for_sporulation = 33.,
                 m_for_sporulation = 1.,
                 n_for_sporulation = 1.,
                                  
                 a_for_emission = 0.71, 
                 b_for_emission = -5.8, 
                 r_for_emission = 0.41, 
                 treshold_nb_du_to_empty = 1,
                 
                 *args, **kwds):
        """ Parameters for powdery mildew.
        
        Parameters
        ----------
        degree_days_to_infect: float
            Degree days needed to achieve infection
            (/!\ personal suggestion)
        basis_for_dday: float
            Base temperature for degree days cumulation
            (/!\ personal suggestion)
        viability_length: float
            Length of the period of viability of a dispersal unit deposited on a leaf
            (/!\ personal suggestion)
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
        beta_for_sporulation: float
            Shape parameter for the calculation of spore production
        gamma_for_sporulation: float
            Shape parameter for the calculation of spore production
        temp_min_for_sporulation: float
            Minimal temperature for sporulation (degrees celsius)
        temp_max_for_sporulation: float
            Maximal temperature for sporulation (degrees celsius)
        m_for_sporulation: float
            Shape parameter for the calculation of the normalized temperature for sporulation
        n_for_sporulation: float
            Shape parameter for the calculation of the normalized temperature for sporulation
        a_for_emission: float
            Shape parameter to compute emission of spores because of wind
        b_for_emission: float
            Shape parameter to compute emission of spores because of wind
        r_for_emission: float
            Shape parameter to compute emission of spores because of wind
        treshold_nb_du_to_empty: int
            Treshold of dispersal unit on a surface to consider it empty
        """
        self.name = "powdery_mildew"
        self.LATENT = LATENT
        self.SPORULATING = SPORULATING
        self.EMPTY = EMPTY
        self.DEAD = DEAD
        self.degree_days_to_infect = degree_days_to_infect
        self.basis_for_dday = basis_for_dday
        self.viability_length = viability_length
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
        self.beta_for_sporulation = beta_for_sporulation
        self.gamma_for_sporulation = gamma_for_sporulation
        
        self.temp_min_for_sporulation = temp_min_for_sporulation
        self.temp_max_for_sporulation = temp_max_for_sporulation
        self.m_for_sporulation = m_for_sporulation
        self.n_for_sporulation = n_for_sporulation
                
        self.a_for_emission = a_for_emission
        self.b_for_emission = b_for_emission
        self.r_for_emission = r_for_emission
        self.treshold_nb_du_to_empty = treshold_nb_du_to_empty

    def __call__(self, nb_spores=None, position=None):
        if PowderyMildew.fungus is None:
            PowderyMildew.fungus = self
        if PowderyMildewDU.fungus is None:
            PowderyMildewDU.fungus = self
        return PowderyMildew(nb_spores=nb_spores, position=position)

def powdery_mildew(**kwds):
    return PowderyMildewParameters(**kwds)

class Disease(object):
    name = 'powdery_mildew'

    @classmethod
    def parameters(cls, **kwds):
        return powdery_mildew(**kwds)
    
    @classmethod
    def dispersal_unit(cls, **kwds):
        PowderyMildewDU.fungus=cls.parameters(**kwds)
        return PowderyMildewDU
    
    @classmethod
    def lesion(cls, **kwds):
        PowderyMildew.fungus=cls.parameters(**kwds)
        return PowderyMildew
    
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