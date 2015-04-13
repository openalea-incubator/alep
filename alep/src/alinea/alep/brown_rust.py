# -*- coding: latin1 -*- 
""" Classes of dispersal unit and lesion of wheat brown rust.
"""
from alinea.alep.fungal_objects import *
import numpy as np

# Dispersal unit ###################################################################################
class BrownRustDU(DispersalUnit):
    """ Define a dispersal unit (or cohort of DUs if group_dus == True) of brown rust """
    def __init__(self, mutable = False):
        """ Initialize the dispersal unit of brown rust """
        super(BrownRustDU, self).__init__(mutable = mutable)
        # Cumulation of temperature conditions
        self.temperature_sequence = []
        # Cumulation of wetness conditions
        self.wetness_sequence = []
    
    def infect(self, dt = 1, leaf = None):
        """ Compute infection success by the dispersal unit 
        
        Response to temperature:
        Pivonia S and Yang X. B. 2006. Relating Epidemic Progress from a General Disease Model
        to Seasonal Appearance Time of Rusts in the United States: Implications for Soybean Rust.
        Phytopathology 96:400-407.
        
        Response to wetness duration:
        Caubel J, Launay M, Lannou C, Brisson N. 2012. Generic response functions to simulate 
        climate-based processes in models for the development of airborne fungal crop pathogens. 
        Ecological Modelling 242: 92–104.
        """
        if self.position is not None and leaf.senesced_length is not None:
            self.disable_if_senescent(leaf = leaf)
        
        f = self.fungus
        # Accumulate climatic data on the leaf sector during the time step
        self.temperature_sequence += leaf.temperature_sequence.tolist()
        self.wetness_sequence += leaf.wetness_sequence.tolist()
        
        # Infection success
        temps = self.temperature_sequence
        wets = self.wetness_sequence
        if len(temps) >= f.infection_delay:
            # Response to temperature
            temp_mean = np.mean(temps)
            beta = (f.temp_max_inf - f.temp_opt_inf)/(f.temp_opt_inf - f.temp_min_inf)
            alpha = 1./((f.temp_opt_inf - f.temp_min_inf)*(f.temp_max_inf - f.temp_opt_inf)**beta)
            temp_factor = max(0., alpha*(temp_mean-f.temp_min_inf)*(f.temp_max_inf - temp_mean)**beta)

            # Response to wetness
            wet_duration = len([w for w in wets if w == True])
            wet_factor = 1 - np.exp(-(f.A_wet_infection*(wet_duration - f.wetness_min))**f.B_wet_infection)

            # Combination
            proba_infection = temp_factor * wet_factor
            dry_duration = len(wets) - wet_duration
            loss_rate = min(1., dry_duration / f.loss_delay if f.loss_delay > 0. else 0.)
            if f.group_dus:
                nb_les = 0
                nb_dead = 0
                for i_du in range(self.nb_dispersal_units):
                    if np.random.random() < proba_infection:
                        nb_les += 1
                    elif np.random.random() < loss_rate:
                        nb_dead += 1
                if nb_les > 0:
                    du = self.fungus.dispersal_unit(mutable = self.mutable)
                    du.set_position(self.position[:nb_les])
                    du.create_lesion(leaf)
                if nb_les + nb_dead > 0:
                    self.set_position(self.position[nb_les+nb_dead:])
                if self.nb_dispersal_units == 0.:
                    self.disable()
                    return
            else:
                if np.random.random() < proba_infection:
                    self.create_lesion(leaf)
                elif np.random.random() < loss_rate:
                    self.disable()
                    return
                
                    dry_duration = len(wets) - wet_duration
                    loss_rate = min(1., dry_duration / f.loss_delay if f.loss_delay > 0. else 0.)
                    if self.fungus.group_dus:
                        nb_disabled = sum([np.random.random() < loss_rate 
                                           for i in range(self.nb_dispersal_units)])
                        self.position = self.position[nb_disabled:]
                        if self.nb_dispersal_units == 0.:
                            self.disable()
                            return
                    else:
                        if np.random.random() < loss_rate:
                            self.disable()
                            return

    def infect_single_du(self, proba_infection = 0.):
        """ Calculate if """
                        
    def disable_if_senescent(self, leaf = None):
        """ Compare position of DU to position of senescent 
            because the DU can only infect green tissues 
        """
        if self.fungus.group_dus == True:
            self.position = filter(lambda x: x[0]>leaf.senesced_length, self.position)
            if self.nb_dispersal_units == 0.:
                self.disable()
                return
        elif self.position[0][0] <= leaf.senesced_length:
            self.disable()
            return
    
    def set_position(self, position=None):
        """ Set the position of the DU to position given in argument
            (force iterable to manage cohorts)
        """
        if position is not None and len(position) > 0 and not is_iterable(position[0]):
            self.position = [position]
        else:
            self.position = position
    
    @property
    def nb_dispersal_units(self):
        """ Get number of dispersal units in cohort if group_dus == True.
            Each DU in cohort has its own position. """
        if self.fungus.group_dus == True:
            return len(self.position) if self.position is not None else 1.
        else:
            return 1.
        
# Lesion ###########################################################################################
class BrownRustLesion(Lesion):
    """ Define a lesion of brown rust """
    def __init__(self, mutable = False):
        """ Initialize the lesion of brown rust """
        super(BrownRustLesion, self).__init__(mutable = mutable)
        # Status of the lesion
        self.status = 0.
        # Age in linear degree days (basis 0°C)
        self.age_tt = 0.
        # Physiological age in latency
        self.age_physio_lat = 0.
        # Surfaces in each stage
        self.surface_lat = 0.
        self.surface_spo = 0.
        self.surface_dead = 0.
        # In case of group_dus == True, check if call of senescence_response is still needed
        self.senescence_response_completed = False
    
    def update(self, dt = 1, leaf = None):
        """ Update the growth demand and the status of the lesion """
        if len(leaf.temperature_sequence) >= self.fungus.latency_min:
            raise ValueError, "Time step is too long, it must be inferior to latency min"
        
        # Manage senescence              
        if (leaf.senesced_length is not None and self.position is not None and
            any([x[0]<=leaf.senesced_length for x in self.position])):
            self.senescence_response(leaf.senesced_length)
            
        if self.is_active:
            # Update status
            if self.is_sporulating:
                if self.surface_lat > 0.:
                    self.latency(leaf_temperature = leaf.temperature_sequence)
                self.sporulation(leaf_temperature = leaf.temperature_sequence)
            else:
                self.latency(leaf_temperature = leaf.temperature_sequence)
                
    def get_effective_temp(self, T = 0.):
        """ Calculate effective temperature on a piecewise linear function of temperature.
            Limit growth and latency for temperatures above optimal.
            
            Effective temperature = piecewise linear function of mean temperature in time step:
            Pivonia S and Yang X. B. 2006. Relating Epidemic Progress from a General Disease Model
            to Seasonal Appearance Time of Rusts in the United States: Implications for Soybean Rust.
            Phytopathology 96:400-407. 
        """
        f = self.fungus
        temp_opt = f.temp_opt_lat
        temp_min = f.temp_min_lat
        temp_max = f.temp_max_lat
        if temp_min < T < temp_opt:
            return temp_opt/((temp_opt - temp_min)/(T - temp_min))
        elif temp_opt <= T < temp_max:
            return temp_opt/((temp_max - temp_opt)/(temp_max - T))
        else:
            return 0.
    
    def delta_thermal_time_growth(self, leaf_temperature = [0.]):
        """ Calculate progress in effective temperature for growth process """
        if len(leaf_temperature) > 0.:
            return sum([max(0, self.get_effective_temp(temp))*1/24. for temp in leaf_temperature])
        else:
            return 0.
    
    def logistic(self, x0, x):
        """ Calculate y for x with logistic curve """
        f = self.fungus
        return f.Smax / (1. + np.exp( -f.k * (x - x0)))
    
    def latency(self, leaf_temperature = [0.]):
        """ Calculate growth demand and surface exchanged to sporulating stage 
        
        Usage of effective temperature to cumulate thermal time (see 'self.get_effective_temp')
        
        Growth demand = logistic function of cumulation of thermal time:
        Audsley, E., A. Milne, and N. Paveley. A Foliar Disease Model for Use in Wheat Disease 
        Management Decision Support Systems. Annals of Applied Biology 147, 
        no. 2 (October 2005): 161–72. 
        
        Surface entering sporulation = f(cumul effective temperatures) parallel to synthesis of 
        tissues entering latency.        
        """
        f = self.fungus
        dtt = self.delta_thermal_time_growth(leaf_temperature = leaf_temperature)
        x1 = self.age_tt
        x2 = x1 + dtt
        
        # Calculate growth demand
        if self.growth_is_active:
            self.growth_demand = self.nb_lesions * (self.logistic(f.x0, x2) - self.logistic(f.x0, x1))
        
        # Calculate progress to sporulation
        to_spo = f.sporulating_capacity * round(self.nb_lesions * 
                (self.logistic(f.x0+f.latency_min, x2) - self.logistic(f.x0+f.latency_min, x1)), 6)
        if not self.growth_is_active:
            to_spo *= self.surface/self._surface_max
        to_spo = min(self.surface_lat, to_spo)
        self.age_tt = x2
        if to_spo > 0.:
            if self.is_latent:
                self.status += 1
            self.surface_lat -= to_spo
            self.surface_spo += to_spo

    def sporulation(self, leaf_temperature = [0.]):
        """ Calculate balance of surface entering sporulation and surface empty,
            and stock of spores 
        """
        pass

    def control_growth(self, growth_offer = 0.):
        """ Limit lesion growth to the surface available on the leaf ('growth_offer')"""
        if growth_offer == 0.:
            self.disable_growth()
            self.growth_demand = 0.
            return
        
        if self.growth_is_active:
            f = self.fungus
            # Assign growth offer
            self.surface_lat += growth_offer
            
            # If lesion has reached max size, disable growth
            if round(self.surface, 4) >= round(self._surface_max, 4):
                self.disable_growth()
            
            self.growth_demand = 0.
    
    def emission(self, emission_rate = 1e4):
        """ Generate a list of dispersal unit emitted by the lesion """
        return []
    
    def senescence_response(self, senesced_length = 0.):
        """ Calculate surface alive and dead after senescence """
        if not self.is_senescent:
            self.become_senescent()
        if not self.senescence_response_completed:
            # Get ratio of lesions senesced in cohort, if individual lesion ratio_sen = 1
            nb_sen = len(filter(lambda x: x[0]<=senesced_length, self.position))
            ratio_sen = float(nb_sen)/(self.nb_lesions)
            # Exchange surfaces between stages
            lat_to_dead = self.surface_lat * ratio_sen
            self.surface_lat -= lat_to_dead
            spo_to_dead = self.surface_spo * ratio_sen
            self.surface_spo -= spo_to_dead
            self.surface_dead += lat_to_dead + spo_to_dead
            # Update number of lesions alive in cohort
            self.position = filter(lambda x: x[0]>senesced_length, self.position)
            if self.nb_lesions == 0.:
                assert self.surface_alive == 0.
                self.senescence_response_completed = True
                self.disable()
    
    def set_position(self, position=None):
        """ Set the position of the lesion to position given in argument
            (force iterable to manage cohorts)
        """
        if position is not None and not is_iterable(position[0]):
            self.position = [position]
        else:
            self.position = position
    
    @property
    def nb_lesions(self):
        """ Get number of lesions in cohort if group_dus == True.
            Each lesion in cohort has its own position. """
        if self.fungus.group_dus == True:
            return len(self.position) if self.position is not None else 1.
        else:
            return 1.
            
    @property
    def is_latent(self):
        """ Check if lesion is latent """
        return self.status == self.fungus.LATENT
    
    @property
    def is_sporulating(self):
        """ Check if lesion is sporulating """
        return self.status == self.fungus.SPORULATING
        
    @property
    def is_empty(self):
        """ Check if lesion is empty """
        return self.status == self.fungus.EMPTY
    
    @property
    def surface_alive(self):
        """ Get surface alive (non senescent) on the lesion """
        return self.surface_lat + self.surface_spo
    
    @property
    def surface(self):
        """ Get total surface of the lesion : non senescent and senescent """
        return self.surface_alive + self.surface_dead
        
    @property
    def _surface_max(self):
        """ Calculate the surface max for a cohort of lesion. """
        return self.fungus.Smax * self.nb_lesions + self.surface_dead
        
# Fungus parameters: config of the fungus ##########################################################
brown_rust_parameters = dict(LATENT = 0,
                             SPORULATING = 1,
                             EMPTY = 2,
                             group_dus = True,
                             temp_opt_inf = 15.,
                             temp_max_inf = 30.,
                             temp_min_inf = 2.,
                             wetness_min = 0.,
                             A_wet_infection = 0.11,
                             B_wet_infection = 3.152,
                             infection_delay = 8.,
                             loss_delay = 48.,
                             Smax = 0.22,
                             k = 0.032,
                             x0 = 254.,
                             temp_opt_lat = 27.,
                             temp_max_lat = 40.,
                             temp_min_lat = 0.,
                             latency_min = 175.,
                             sporulating_capacity = 0.2)
                             
# Brown rust fungus ################################################################################
class BrownRustFungus(Fungus):
    """ Define a fungus model with dispersal unit class, lesion class and set of parameters """
    def __init__(self, name = 'brown_rust', 
                 Lesion = BrownRustLesion,
                 DispersalUnit = BrownRustDU, 
                 parameters = brown_rust_parameters):
        super(BrownRustFungus, self).__init__(name = name, 
                                              Lesion = Lesion, 
                                              DispersalUnit = DispersalUnit,
                                              parameters = parameters)

import collections
def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)
                                              
def get_proba_inf(T):
    # Pivonia and Yang
    temp_opt = 15.
    temp_max = 30.
    temp_min = 2.
    T = min(temp_max, T)
    beta = (temp_max - temp_opt)/(temp_opt - temp_min)
    alpha = 1/((temp_opt - temp_min)*(temp_max - temp_opt)**beta)
    return max(0, alpha*(T-temp_min)*(temp_max - T)**beta)
    
def get_progress(T):
    latency_min = 144.
    return get_efficency(T)/latency_min

def get_effective_T(T):
    # Pivonia and Yang
    temp_opt = 27.
    temp_max = 40.
    temp_min = 0.
    if temp_min < T < temp_opt:
        return 1./((temp_opt - temp_min)/(T - temp_min))
    elif temp_opt <= T < temp_max:
        return 1./((temp_max - temp_opt)/(temp_max - T))
    else:
        return 0.
        
def plot_progress():
    import matplotlib.pyplot as plt
    temp = range(40)
    plt.plot(temp, [get_progress(T) for T in temp])
        
def logistic_growth(x, x0 = 254.):
    # Inspire de Audsley
    Kmax = 1.
    k = 0.032
    return Kmax/(1+np.exp(-k*(x - x0)))
    
def plot_logistic_growth():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(range(450), [logistic_growth(x) for x in range(450)])
    # Points de Audsley
    ax.plot([162], [0.05], 'ro')
    ax.plot([324], [0.9], 'ro')
    
def compare_logistic_growth_and_spo(nb_steps = 1000):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(range(nb_steps), [logistic_growth(x) for x in range(nb_steps)])
    ax.plot(range(nb_steps), [logistic_growth(x, 254.+175.) for x in range(nb_steps)])
    # Points de Audsley
    ax.plot([162], [0.05], 'ro')
    ax.plot([324], [0.9], 'ro')