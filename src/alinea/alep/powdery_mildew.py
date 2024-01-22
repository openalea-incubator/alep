""" Classes of dispersal unit, lesion of grapevine powdery mildew.

..:todo: - adapt emission to new protocol
         - build PowderyMildewInspector

"""
# Imports #########################################################################
from alinea.alep.fungus import *
import numpy as np
from random import random, randint
from math import (exp, pi, floor, ceil, sqrt)

# Dispersal unit ##################################################################
class PowderyMildewDU(DispersalUnit):
    """ Define a dispersal unit specific of powdery mildew.
    
    """
    fungus = None
    def __init__(self, mutable=False):
        """ Initialize the dispersal unit (DU).
        
        :Parameters:
         - 'mutable' (bool) - True if each DU has its own parameters (intra-population variability),
         False if all DU of the same fungus share the same parameters.
        """
        super(PowderyMildewDU, self).__init__(mutable=mutable)
        # Age of the dispersal unit in degree days
        self.age_dday = 0.
        # Viability of the DU : Decreases with time down to 0
        self.viability = 1.
        
    def infect(self, dt=1, leaf=None):
        """ Compute infection by the dispersal unit of powdery mildew.
        
        :Parameters:
         - 'dt' (int): Time step of the simulation (in hours)
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
        """
        if self.is_active:
                
            # External variables
            if leaf.green_area== 0.:
                self.disable()
                return
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
            
            if self.nb_dispersal_units == 0.:
                self.disable()
                return
            
            # Compute dispersal unit age (and viability?)
            self.update_age_dday(dt, leaf)

            if self.is_ready_to_infect():
                # Temperature factor
                t_norm_function = temp_norm_function(temp, t_min, t_max, m, n)
                temp_factor = (max_rate * t_norm_function * exp(-decay_rate * leaf.age))
                
                # Relative humidity factor
                RH_factor = min(1., a_RH_effect * relative_humidity + b_RH_effect)
                
                # Wetness factor
                if relative_humidity >= RH_opt or leaf_wet:
                    wetness_factor = min(1., c_wetness_effect - d_wetness_effect * temp)
                else:
                    wetness_factor = 1.
                
                # Probability of infection
                proba_infection = temp_factor * RH_factor * wetness_factor * spores_factor
                proba_loss = 1 - self.viability
                if f.group_dus:
                    init_nb_dus = self.nb_dispersal_units
                    nb_les = np.random.binomial(init_nb_dus, proba_infection)
                    self.create_lesion(nb_les, leaf)
                    if init_nb_dus > nb_les:
                        nb_dead = np.random.binomial(init_nb_dus - nb_les, proba_loss)
                        self.nb_dispersal_units -= nb_dead
                else:
                    if np.random.random() < proba_infection:
                        self.create_lesion(1, leaf)
                    elif np.random.random() < proba_loss:
                        self.disable()
                
    def update_age_dday(self, dt=1., leaf=None):
        """ Update the age of the lesion.
        
        :Parameters:
         - 'dt' (int): Time step of the simulation (in hours)
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
        """
        self.age_dday += self.compute_delta_ddays(dt, leaf)
        
    def compute_delta_ddays(self, dt=1., leaf=None):
        """ Compute delta degree days in dt.
        
        :Parameters:
         - 'dt' (int): Time step of the simulation (in hours)
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
        
        :Returns:
         - 'ddday' (float): Delta degree days in time step
        """
        f = self.fungus
        if dt != 0.:
            ddday = max(0,(leaf.temp - f.basis_for_dday*dt)/(24./dt))
        else:
            ddday = 0.
        return ddday 
        
    def update_viability(self, dt=1., leaf=None):
        """ Update the viability of the dispersal unit.
        
        The hypothesis is made that a dispersal unit is viable 5 days (120h).
        The decrease of viability by time step is doubled if the leaf is wet.
        
        :Parameters:
         - 'dt' (int): Time step of the simulation (in hours)
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
        """
        f = self.fungus
        if leaf.wetness:
            dt *= 2.
        viability_loss = dt/f.viability_length
        self.viability -= viability_loss if self.viability > viability_loss else 0.
    
    @property
    def is_viable(self):
        """ Check if the dispersal unit is still viable.
        """
        return self.viability > 0.

    @property
    def is_ready_to_infect(self):
        """ Check if the dispersal unit is ready to infect according to its age.
        """
        f = self.fungus
        return self.age_dday >= f.degree_days_to_infect

# Lesion ##########################################################################
class PowderyMildewLesion(Lesion):
    """ Powdery mildew lesion implemented as in Calonnec et al., 2008 for the most part. """
    
    def __init__(self, nb_lesions=1,nb_spores=None, position=None, mutable=False):
        """ Initialize the lesion of powdery mildew.
        
        Parameters
        ----------
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        position: non defined
            Position of the dispersal unit on the phyto-element
        """
        super(PowderyMildewLesion, self).__init__(mutable=mutable)
        self.position = position
        self.nb_spores = nb_spores
        self.nb_lesions=nb_lesions
        f = self.fungus
        # Duration of the time step
        self.dt = 0.
        # Age of the lesion (in hours)
        self.age = 0.
        # Surface of the lesion
        self.surface = 0.
        # Diameter of the lesion
        self.diameter = 0.
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
        
        # Temporary variables for model inspection
        self.age_leaf_at_infection = None
        self.surface_evolution = []        
        self.temperatures = []        
    
    def update(self, dt=1., leaf=None):
        """ Update the status of the lesion and calculate growth demand.
 
        Parameters
        ----------
        :Parameters:
         - 'dt' (int): Time step of the simulation (in hours)
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
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
    
        # Temporary:
        if self.age_leaf_at_infection == None:
            self.age_leaf_at_infection = leaf.age
        self.temperatures.append(leaf.temp)
    
    def update_diameter_max(self, leaf=None):
        """ Compute maximum diameter of the lesion or each lesion in cohort according to its age.
        
        :Parameters:
         - 'dt' (int): Time step of the simulation (in hours)
         - 'leaf' (Leaf sector node of a MTG): A leaf sector with properties (e.g. area, 
            green area, senescence, wetness, temperature, other DUs and lesions, ...)
        """
        # Parameters for the calculation
        f = self.fungus
        diameter_max = f.diameter_max
        diameter_min = f.diameter_min
        leaf_age_effect = f.leaf_age_effect
        
        # Maximum diameter of a lesion in cohort according to leaf age
        kmax = diameter_min + (diameter_max - diameter_min) * exp(-leaf_age_effect * leaf.age)
        self.diameter_max = kmax
        
    def update_growth_demand(self):
        """ Compute lesion growth demand for the time step.
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
        age = self.age * dt/24.
        t_norm_function = temp_norm_function(temp, t_min, t_max, m, n)
        age_factor = (dt/24.)*r*exp(r*(half_time-age))/(1+exp(r*(half_time-age)))**2
        diameter_demand = diameter + kmax * t_norm_function * age_factor
        
        # Growth demand in surface
        self.growth_demand = self.nb_lesions * (pi * diameter_demand**2 /4) - self.surface
    
    def progress_in_latency(self):
        """ Compute the ageing of the lesion in latency period and 
            set the status of the lesion to SPORULATING when needed.
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
        
        :Parameters:
         - 'growth_offer' (float): Minimum between 'growth_demand' and the surface available on
            the leaf for the lesion to grow (cm2)
        """
        # 1. Update surface
        if self.growth_is_active:
            # Growth offer is added to surface alive
            self.surface += growth_offer
            self.diameter = sqrt(self.surface*4/pi)
            
            # Temporary
            self.surface_evolution.append(self.surface)
            
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
        
        Not used in Calonnec model.
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
        
        Not used in Calonnec model.
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
        max_rate = beta*exp(gamma*surface)*dt/24
        t_norm_function = temp_norm_function(temp, t_min, t_max, m, n)
        production = max_rate * t_norm_function
        
        self.stock_spores += production
    
    def emission(self, leaf=None):
        """ Create a list of dispersal units emitted by the ring.
        
        Computed with the formalism of Willocquet et al., 1998.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. area, green area, healthy area,
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
                emissions = [PowderyMildewDU(nb_spores=nb_spores_by_DU, status='emitted',
                                        position=self.position) for i in range(nb_DU_emitted)]
                # Update dtock of spores                       
                self.stock_spores -= nb_spores if self.stock_spores>nb_spores else 0.
        
        # Disable the lesion if not producing spore anymore and stock is empty
        if (self.is_sporulating() and 
            not self.production_is_active and
            self.stock_spores==0.):
            self.disable()

        return emissions
    
    def reduce_stock(self, nb_spores_emitted):
        """ Reduce the stock of spores after emission.
        
        Parameters
        ----------
        nb_spores_emitted: int
            Number of spores emitted
        """
        self.stock_spores -= nb_spores_emitted
        if self.stock_spores < self.fungus.treshold_spores:
            # print("coucou")
            self.stock_spores = 0.
    
        # Disable the lesion if not producing spore anymore and stock is empty
        if (self.is_sporulating() and 
            not self.production_is_active and
            self.stock_spores==0.):
            self.status = self.fungus.EMPTY
            self.disable()
    
    def become_senescent(self, **kwds):
        """ Set is_senescent to true.
        """
        self.is_senescent = True
    
    def senescence_response(self):
        """ Disable the lesion if it is on senescent tissue.
        """
        self.disable()
        
    def is_latent(self):
        """ Check if the status of the lesion is LATENT.
        """
        f = self.fungus
        return self.status==f.LATENT
    
    def is_sporulating(self):
        """ Check if the status of the lesion is SPORULATING.
        """
        f = self.fungus
        return self.status==f.SPORULATING

    def is_empty(self):
        """ Check if the status of the lesion is EMPTY.
        """
        f = self.fungus
        return self.status==f.EMPTY    
        
    def disable_production(self):
        """ Shut down lesion production activity (turn it to False)
        """
        self.production_is_active = False

    def activate_production(self):
        """ Turn lesion production activity on (turn it to True)
        """
        self.production_is_active = True
        
# Fungus parameters: config of the fungus ##########################################################
powdery_mildew_parameters = dict(name='powdery_mildew',
                                 LATENT = 0,
                                 SPORULATING = 1,
                                 EMPTY = 2,
                                 DEAD = 3,
                                 group_dus = True,
                                 degree_days_to_infect = 20.,
                                 basis_for_dday = 0.,
                                 viability_length = 120.,
                                 temp_min_for_infection = 5.,
                                 temp_max_for_infection = 33.,
                                 m_for_infection = 0.338,
                                 n_for_infection = 1.055,
                                 max_infection_rate = 0.53,
                                 decay_rate = 0.147,
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
                                 growth_rate = 0.2,
                                 half_growth_time = 13., 
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
                                 treshold_spores = 1)
                                 
# Powdery mildew fungus ############################################################################
class PowderyMildewFungus(Fungus):
    """ Define a fungus model with dispersal unit class, lesion class and set of parameters """
    def __init__(self,
                 Lesion = PowderyMildewLesion,
                 DispersalUnit = PowderyMildewDU,
                 parameters = powdery_mildew_parameters):
        super(PowderyMildewFungus, self).__init__(Lesion = Lesion,
                                              DispersalUnit = DispersalUnit,
                                              parameters = parameters)
    
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