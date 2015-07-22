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
        # Number of dispersal_units
        self.nb_dispersal_units = 1

    def infect(self, dt = 1, leaf = None):
        """ Compute infection success by the dispersal unit

        Response to temperature:
        Pivonia S and Yang X. B. 2006. Relating Epidemic Progress from a General Disease Model
        to Seasonal Appearance Time of Rusts in the United States: Implications for Soybean Rust.
        Phytopathology 96:400-407.

        Response to wetness duration:
        Caubel J, Launay M, Lannou C, Brisson N. 2012. Generic response functions to simulate
        climate-based processes in models for the development of airborne fungal crop pathogens.
        Ecological Modelling 242: 92â??104.
        """
        if self.is_active:
            f = self.fungus
            if leaf.green_area== 0.:
                self.disable()
                return
                
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
                wet_factor = np.exp(-f.B_wet_infection * np.exp(-f.k_wet_infection*wet_duration))
                dry_duration = len(wets) - wet_duration
                loss_rate = min(1., dry_duration / f.loss_delay if f.loss_delay > 0. else 0.)
    
                # Combination
                proba_infection = temp_factor * wet_factor * f.proba_inf
                if f.group_dus:
                    init_nb_dus = self.nb_dispersal_units
                    nb_les = np.random.binomial(init_nb_dus, proba_infection)
                    self.create_lesion(nb_les, leaf)
                    if init_nb_dus > nb_les:
                        nb_dead = np.random.binomial(init_nb_dus - nb_les, loss_rate)
                        self.nb_dispersal_units -= nb_dead                
                        if self.nb_dispersal_units == 0.:
                            self.disable()
                            return
                else:
                    if np.random.random() < proba_infection:
                        self.create_lesion(1, leaf)
                    elif np.random.random() < loss_rate:
                        self.disable()
                        return

    def create_lesion(self, nb_lesions = 1, leaf = None, **kwds):
        green_length = leaf.green_length
        length = leaf.length
        if green_length>0 and nb_lesions>0:        
            les = self.fungus.lesion(mutable = self.mutable)
            les.__dict__.update(kwds)
            les.set_position([[length - np.random.random()*green_length, 0] 
                                for i in range(nb_lesions)])
            self.nb_dispersal_units -= nb_lesions
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
        else:
            self.disable()
    
    def set_status(self, status = 'deposited'):
        self.status = status
        
    def set_nb_dispersal_units(self, nb_dispersal_units = 1):
        self.nb_dispersal_units = nb_dispersal_units

# Lesion ###########################################################################################
class BrownRustLesion(Lesion):
    """ Define a lesion of brown rust """
    def __init__(self, mutable = False):
        """ Initialize the lesion of brown rust """
        super(BrownRustLesion, self).__init__(mutable = mutable)
        
        # TEMP
        if self.nb_lesions >1:
            import pdb
            pdb.set_trace()
        
        # Status of the lesion
        self.status = 0.
        # Age in thermal time
        self.age_tt = 0.
        # Age in sporulation
        self.age_sporulation = 0.
        # Surfaces in each stage
        self.surface_chlo = 0.
        self.surface_spo = 0.
        self.surface_sink = 0.
        self.surface_empty = 0.
        self.surface_dead = 0.
        # In case of group_dus == True, check if call of senescence_response is still needed
        self.senescence_response_completed = False
        # Number of lesions in senescence
        self.nb_lesions_sen = 0
        # Stock of spores on lesion
        self.stock_spores = 0.
        # Potential surface if no compeptition
        self.potential_surface = 0.
        
    def update(self, dt = 1, leaf = None):
        """ Update the growth demand and the status of the lesion """
        # Eliminate extra lesions
        if round(self.surface, 10) == 0. and not self.growth_is_active:
            self.disable()
        
        # Manage senescence
        if (leaf.senesced_length is not None and self.position is not None and
            any([x[0]<=leaf.senesced_length for x in self.position])):
            self.senescence_response(leaf.senesced_length)

        if self.is_active:
            # Calculate progress in thermal time
            self.update_age_and_status(leaf_temperature = leaf.temperature_sequence)
            
            # Calculate growth demand
            self.update_growth()         
            
            # Calulate production of spores and necrosis
            if self.is_sporulating:
                self.update_sporulation(leaf)
                
            # Update potential surface
            self.potential_surface += self.growth_demand


    def get_effective_temp(self, T = 0.):
        """ Calculate effective temperature on a piecewise linear function of temperature.
            Limit growth and latency for temperatures above optimal.

            Effective temperature = piecewise linear function of mean temperature in time step:
            Pivonia S and Yang X. B. 2006. Relating Epidemic Progress from a General Disease Model
            to Seasonal Appearance Time of Rusts in the United States: Implications for Soybean Rust.
            Phytopathology 96:400-407.
        """
        f = self.fungus
        temp_opt = f.temp_opt_chlo
        temp_min = f.temp_min_chlo
        temp_max = f.temp_max_chlo
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

    def update_age_and_status(self, leaf_temperature = [0.]):
        self.dtt = self.delta_thermal_time_growth(leaf_temperature = leaf_temperature)
        self.age_tt += self.dtt
        
        if (self.is_chlorotic and self.age_tt >= self.fungus.latency):
            self.status += 1
            
    def logistic(self, x0, x):
        """ Calculate y for x with logistic curve """
        f = self.fungus
        return f.Smax / (1. + np.exp( -f.k * (x - x0)))

    def update_growth(self):
        f = self.fungus
        x1 = self.age_tt - self.dtt
        x2 = self.age_tt
        self.growth_demand = self.nb_lesions_non_sen * (self.logistic(f.x0, x2) - \
                        self.logistic(f.x0, x1))
        
    def update_sporulation(self, leaf):
        f = self.fungus
        self.age_sporulation += self.dtt        
        self.update_necrosis(leaf)
        self.stock_spores += self.surface_spo * f.production_max_by_tt * \
                                 self.dtt * f.conversion_mg_to_nb_spo
        
    def update_necrosis(self, leaf=None):
           
        def gompertz_necro(date, dens):
            a = 3.08e-5
            b = 0.00305171209440162
            c = 0.0089309782201911596
            d = 10.448840245012178
            A = a*dens + b
            B = c*dens + d
            return np.exp(-B * np.exp(-A*date)) 
        
        if self.surface_sink > 0.:
            f = self.fungus
            a_spo = self.age_sporulation
            nb_les = sum([l.nb_lesions_non_sen for l in leaf.lesions])
            dens = nb_les/leaf.green_area if leaf.green_area > 0 else 0.
            ratio_empty = gompertz_necro(a_spo, dens)-gompertz_necro(a_spo-self.dtt, dens)
            ratio_les = self.nb_lesions_non_sen/float(nb_les)
            total_smax = self.nb_lesions_non_sen*(1-np.exp(-dens*f.Smax))/dens
            smax = ratio_les*total_smax
            empty_sink = min(self.surface_sink, smax * f.ratio_sink * ratio_empty)
            empty_chlo =  min(self.surface_chlo, smax * f.ratio_chlo * ratio_empty)
            empty_spo =  min(self.surface_spo, smax * f.ratio_spo * ratio_empty)

            self.surface_sink -= empty_sink
            self.surface_chlo -= empty_chlo
            self.surface_spo -= empty_spo
            self.surface_empty += empty_sink + empty_chlo + empty_spo


    def control_growth(self, growth_offer = 0.):
        """ Limit lesion growth to the surface available on the leaf ('growth_offer')"""
        if self.growth_is_active:
            # Assign growth offer
            f = self.fungus
            if growth_offer >= 0:
                self.surface_sink += (1 - f.ratio_chlo - f.ratio_spo) * growth_offer
                self.surface_chlo += f.ratio_chlo * growth_offer
                self.surface_spo += f.ratio_spo * growth_offer
            else:
                if self.surface_alive>0:
                    self.surface_sink += growth_offer*self.surface_sink/self.surface_alive
                    self.surface_chlo += growth_offer*self.surface_chlo/self.surface_alive                 
                    self.surface_spo += growth_offer*self.surface_spo/self.surface_alive  
                    self.surface_empty += growth_offer*self.surface_empty/self.surface_alive  
                    self.surface_dead -= growth_offer

            # If lesion has reached max size, disable growth
            if round(self.surface, 4) >= round(self._surface_max, 4):
                self.disable_growth()

            self.growth_demand = 0.

    def emission(self, emission_rate = 1e4):
        """ Generate a list of dispersal unit emitted by the lesion """
        nb_dus =  int(self.stock_spores)
        self.stock_spores = 0.
        return nb_dus

    def senescence_response(self, senesced_length = 0.):
        """ Calculate surface alive and dead after senescence """
        if not self.is_senescent:
            self.become_senescent()
        if not self.senescence_response_completed:
            # Get ratio of lesions senesced in cohort, if individual lesion ratio_sen = 1
            nb_sen = len(filter(lambda x: x[0]<=senesced_length, self.position))
            nb_new_sen = nb_sen - self.nb_lesions_sen
            ratio_sen = float(nb_new_sen)/(self.nb_lesions_non_sen)

            # Exchange surfaces between stages
            sink_to_dead = self.surface_sink * ratio_sen
            self.surface_sink -= sink_to_dead
            chlo_to_dead = self.surface_chlo * ratio_sen
            self.surface_chlo -= chlo_to_dead
            spo_to_dead = self.surface_spo * ratio_sen
            self.surface_spo -= spo_to_dead
            empty_to_dead = self.surface_empty * ratio_sen
            self.surface_empty -= empty_to_dead
            self.surface_dead += sink_to_dead + chlo_to_dead + spo_to_dead + empty_to_dead
            
            # Update potential surface
            self.potential_surface *= (1-ratio_sen)
            
            # Update number of lesions in senescence
            self.nb_lesions_sen += nb_new_sen  
            
            # Stop developemnt when last lesion of cohort is reached
            if self.nb_lesions_non_sen == 0.:
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
            
    def disable_growth(self):
        """ Disable growth of the lesion.
        """
        self.growth_is_active = False
        self.growth_demand = 0.
        if self.surface < 1e-5:
            self.surface_chlo = 0.
            self.surface_spo = 0.
            self.surface_sink = 0.
            self.surface_empty = 0.
            self.surface_dead = 0.
            self.disable()

    def disappear(self):
        """ Kill the lesion and pass all surfaces to 0. 
        """
        self.growth_is_active = False
        self.growth_demand = 0.
        self.surface_chlo = 0.
        self.surface_spo = 0.
        self.surface_sink = 0.
        self.surface_empty = 0.
        self.surface_dead = 0.
        self.disable()
            
    @property
    def nb_lesions(self):
        """ Get number of lesions in cohort if group_dus == True.
            Each lesion in cohort has its own position. """
        if self.fungus.group_dus == True:
            return len(self.position) if self.position is not None else 1.
        else:
            return 1.
            
    @property
    def nb_lesions_non_sen(self):
        """ Get number of non senescent lesions in cohort. """
        if self.position is None:
            return None
        else:
            return self.nb_lesions - self.nb_lesions_sen

    @property
    def is_chlorotic(self):
        """ Check if lesion is chlorotic """
        return self.status == self.fungus.CHLOROTIC

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
        return self.surface_sink + self.surface_chlo + self.surface_spo + self.surface_empty

    @property
    def surface(self):
        """ Get total surface of the lesion : non senescent and senescent """
        return self.surface_alive + self.surface_dead

    @property
    def _surface_max(self):
        """ Calculate the surface max for a cohort of lesion. """
        return self.fungus.Smax * self.nb_lesions_non_sen + self.surface_dead
        
    @property
    def surface_non_senescent(self):
        """ calculate the surface of the lesion non affected by senescence. """
        return self.surface_alive
        
# Fungus parameters: config of the fungus ##########################################################
brown_rust_parameters = dict(name = 'brown_rust',
                             CHLOROTIC = 0,
                             SPORULATING = 1,
                             EMPTY = 2,
                             group_dus = True,
                             proba_inf = 0.1,
                             temp_opt_inf = 15.,
                             temp_max_inf = 30.,
                             temp_min_inf = 2.,
                             wetness_min = 0.,
                             k_wet_infection = 0.47,
                             B_wet_infection = 30.,
                             infection_delay = 3.,
                             loss_delay = 48.,
                             Smax = 0.09,
                             k = 0.015,
                             x0 = 350.,
                             temp_opt_chlo = 27.,
                             temp_max_chlo = 40.,
                             temp_min_chlo = 0.,
                             ratio_sink = 0.33, 
                             ratio_chlo = 0.4,
                             ratio_spo = 0.27,
                             latency = 170.,
                             sporulating_capacity = 0.7,
                             production_max_by_tt = 0.05,
                             conversion_mg_to_nb_spo = 2e5,
                             delay_high_spo = 720.,
                             delay_total_spo = 1050.)

# Brown rust fungus ################################################################################
class BrownRustFungus(Fungus):
    """ Define a fungus model with dispersal unit class, lesion class and set of parameters """
    def __init__(self,
                 Lesion = BrownRustLesion,
                 DispersalUnit = BrownRustDU,
                 parameters = brown_rust_parameters):
        super(BrownRustFungus, self).__init__(Lesion = Lesion,
                                              DispersalUnit = DispersalUnit,
                                              parameters = parameters)

import collections
def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)
    
def group_duplicates_in_cohort(g):
    def _get_index_duplicates(seq):
        from collections import defaultdict
        dd = defaultdict(list)
        for i,item in enumerate(seq):
            dd[item].append(i)
        return [idx for key,idx in dd.iteritems() if len(idx)>1 and key==0.][-1]
    
    def group_lesions(les):
        new_l = les[0].fungus.lesion()
        new_l.position = sum([l.position for l in les], [])
        return new_l
    
    lesions = g.property('lesions')
    for vid, les in lesions.iteritems():
        ages = [l.age_tt for l in les]
        if len(les)!=len(set(ages)):
            idxs = _get_index_duplicates(ages)
            new_les = group_lesions([les[i] for i in idxs])
            les = les[:idxs[0]]
            les.append(new_les)
            lesions[vid] = les
    
def get_proba_inf_T(T):
    # Pivonia and Yang
    temp_opt = 15.
    temp_max = 30.
    temp_min = 2.
    T = min(temp_max, T)
    beta = (temp_max - temp_opt)/(temp_opt - temp_min)
    alpha = 1/((temp_opt - temp_min)*(temp_max - temp_opt)**beta)
    return max(0, alpha*(T-temp_min)*(temp_max - T)**beta)

def plot_proba_inf_T():
    import matplotlib.pyplot as plt
    temp = range(36)
    fig, ax = plt.subplots()
    ax.plot(temp, [get_proba_inf_T(T) for T in temp])
    ax.set_ylabel("Probability of infection", fontsize = 16)
    ax.set_xlabel("Temperature (degrees Celsius)", fontsize = 16)

def weibull(x, A = 0.11, B = 3.152, x_min = 0.):
    return 1 - np.exp(-(A*x - x_min)**B)
    
def gompertz(x, k = 0.47, B = 30.):
    return np.exp(-B * np.exp(-k*x))
    
def get_proba_inf_WD(WD, A = 0.11, B = 3.152):
    # Caubel
    wetness_min = 8
    if WD >= wetness_min:
        return weibull(WD, A = A, B = B, x_min = 0.)
    else:
        return 0.

def compare_proba_inf_WD_weibull(A = 0.11, B = 3.152):
    import matplotlib.pyplot as plt
    wetness = np.arange(0,30,0.1)
    y1 = [get_proba_inf_WD(WD) for WD in wetness]
    y2 = [weibull(WD, A = A, B = B, x_min = 0.) for WD in wetness]
    fig, ax = plt.subplots()
    ax.plot(wetness, y1, 'b')
    ax.plot(wetness, y2, 'r')
    
def compare_proba_inf_WD_gomp(k = 0.47, B = 30.):
    import matplotlib.pyplot as plt
    wetness = np.arange(0,30,0.1)
    y1 = [get_proba_inf_WD(WD) for WD in wetness]
    y2 = [gompertz(WD, k=k, B=B) for WD in wetness]
    fig, ax = plt.subplots()
    ax.plot(wetness, y1, 'b')
    ax.plot(wetness, y2, 'r')

def plot_proba_inf_WD():
    import matplotlib.pyplot as plt
    wetness = np.arange(0,30,0.1)
    fig, ax = plt.subplots()
    ax.plot(wetness, [gompertz(WD) for WD in wetness])
    ax.set_ylabel("Probability of infection", fontsize = 16)
    ax.set_xlabel("Wetness duration (hours)", fontsize = 16)

def get_progress(T):
    delay_to_chlo = 144.
    return get_effective_T(T)/delay_to_chlo

def get_effective_T(T):
    # Pivonia and Yang
    temp_opt = 27.
    temp_max = 40.
    temp_min = 0.
    if temp_min < T < temp_opt:
        return temp_opt/((temp_opt - temp_min)/(T - temp_min))
    elif temp_opt <= T < temp_max:
        return temp_opt/((temp_max - temp_opt)/(temp_max - T))
    else:
        return 0.

def plot_effective_T():
    import matplotlib.pyplot as plt
    temp = range(41)
    fig, ax = plt.subplots()
    ax.plot(temp, [get_effective_T(T) for T in temp])
    ax.set_ylabel("Teff (degrees Celsius)", fontsize = 16)
    ax.set_xlabel("Temperature (degrees Celsius)", fontsize = 16)

def plot_progress():
    import matplotlib.pyplot as plt
    temp = range(40)
    fig, ax = plt.subplots()
    ax.plot(temp, [get_progress(T) for T in temp])

def logistic_growth(x, x0 = 254.):
    # Inspire de Audsley
    #(parametre Kmax deduit de donnees de Robert 2004 apres test sur competition)
    Kmax = 0.22
    k = 0.032
    return Kmax/(1+np.exp(-k*(x - x0)))

def plot_logistic_growth():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(range(450), [logistic_growth(x) for x in range(450)])
    # Points de Audsley
    # ax.plot([162], [0.05*0.22], 'ro')
    # ax.plot([324], [0.9*0.22], 'ro')
    ax.set_ylabel("Surface of 1 lesion (cm2)", fontsize = 18)
    ax.set_xlabel("Thermal time (Teff)", fontsize = 18)

def compare_logistic_growth_and_spo(nb_steps = 1000):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(range(nb_steps), [logistic_growth(x) for x in range(nb_steps)])
    ax.plot(range(nb_steps), [logistic_growth(x, 254.+175.) for x in range(nb_steps)])
    # Points de Audsley
    ax.plot([162], [0.05*0.22], 'ro')
    ax.plot([324], [0.9*0.22], 'ro')

