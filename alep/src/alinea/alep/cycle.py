# -*- coding: latin1 -*- 

##
##

""" A Lesion is composed by successive growth Rings.

A ring is defined as a multi-state automaton.
The different states are:
    * DEPOSIT
    * EMERGENT
    * INCUBATING
    * CHLOROTIC
    * SPORULATING
    * EMPTY
    * DEAD    
    
The Lesion class can call 2 types of Ring class: Septoria and PowderyMildew.

"""
from random  import random
from math import exp
from math import pi

# Development stages of the lesion
DEPOSIT = 0
EMERGENT = 1

# TODO G.Garin, 24/09/2012:
# Ne pas laisser les variables globales d'etat des rings a cet endroit.
# Doit etre defini dans chaque maladie. La lesion ne doit pas les connaitre.
# OK pour etats apres emergents. Pb pour etats avant : Cf initialisation d'un 
# nouvel anneau.

class Lesion(object):
    """ 
    """
    
    def __init__(self, fungus):
        """ Initialize the lesion. 
        
        :Parameters:
          - `fungus` (function): returns a class of specific parameters for 
          the chosen fungus (e.g. 'septoria()' or 'powderyMildew()').

        """
        self.rings = []
        self.fungus = fungus
    
    def create(self):
        """ Create a new lesion.
        
        :Parameters:
          - `fungus` (function): returns a class of specific parameters for 
          the chosen fungus (e.g. 'septoria()' or 'powderyMildew()').

        """
        new_ring = self.fungus_factory()
        new_ring.status = EMERGENT
        self.rings.append(new_ring)
        # TODO G.Garin, 23/11/2012:
        # Changer car ne doit plus fonctionner avec des anneaux. 
        # Doit aussi pouvoir etre plus parametrable (ie construire des lesions
        # d'un etat donne, d'un age donne et d'une surface donnee).
        # En outre, si c'est la lesion qui gere la reussite de l'infection,
        # alors elle devra appeler la fonction 'create'.
    
    def update(self, dt, leaf, environment, **kwds):
        """ Update the status of the lesion and create a new growth ring if needed.
        
        If it has enough space, a lesion grows in a concentric way until it reaches
        its maximum size. At each time step, it forms a new ring of mycelium. 
        This ring passes through different stages of an automaton (see 'Ring' class).
        The status of the lesion is updated according the status of its first ring.
         
        :Parameters:
          - `dt` (float): delta time.
          - `leaf` (class): a leaf sector with properties (e.g. healthy surface,
          senescence, rain intensity, wetness, temperature, lesions).

        """        
        environment['lesion'] = self
        
        # Update the status of each ring
        for ring in self.rings:
            ring.update(dt, leaf, environment=environment, **kwds)
            
        # Create a new ring 
        if self.can_form_new_ring(leaf, environment):
            new_ring = self.fungus_factory()
            if self.status >= EMERGENT:
                new_ring.status = EMERGENT

            self.rings.append(new_ring)

        # A etoffer (arret croissance a cause de la senescence, ...?)
        
        # TODO G.Garin, 24/09/2012:
        # Probleme de la variable globale 'EMERGENT'. 
        
    def is_dead(self):
        """ Update the status of all the rings to 'DEAD' if the lesion is dead. """
        return all(ring.is_dead() for ring in self.rings)
         
    def can_form_new_ring(self, leaf, environment):
        """ Check if the lesion can form a new ring.
        
        A new ring is not created :
            * if the first du (which is modelled as a ring) is dead
            * if the infection is not achieved
            * during the time step when the infection occurs
            * if the surface max of the lesion has been reached
            * if there is no green surface on the leaf
            
        """
        
        if self.fungus.name == "Septoria":
            return Septoria.can_form_new_ring(Septoria(), leaf, environment)
        
        elif self.fungus.name == "PowderyMildew":
            return PowderyMildew.can_form_new_ring(PowderyMildew(), leaf, environment)
            
        # TODO G.Garin, 24/09/2012:
        # Supprimer la fonction 'can_form_new_ring' de la classe lesion.
        # La mettre dans chaque classe de maladie. Il ne doit pas y avoir de 
        # 'if nom_du_champignon = "champignon"' dans la classe lesion.

    def fungus_factory(self):
        """ This factory will search for entry points.
        
        Thus, we will not have to change this file to add a new fungus.
        """
        if self.fungus.name == "Septoria":
            return Septoria()
        elif self.fungus.name == "PowderyMildew":
            return PowderyMildew()
            
    @property
    def surface(self):
        """ Compute the surface of the lesion.
        """
        surf = sum(ring.surface for ring in self.rings)
        return surf

    @property
    def status(self):
        """ Compute the status of the lesion.
        """
        if self.rings:
            return self.rings[0].status
    
    @property
    def age(self):
        """ Compute the age of the lesion.
        """
        if self.rings:
            return self.rings[0].age
            
    @property
    def age_dday(self):
        """ Compute the thermal age of the lesion.
        """
        if self.rings:
            return self.rings[0].age_dday

    @status.setter
    def status(self, value):
        """ Set the status of the lesion to the chosen value.
        """
        if self.rings:
            self.rings[0].status = value

 ####################################################################################################
def proba(p):
    """ p in 0,1 """
    return random() < p

class Ring(object):
    """ Ring of Lesion at a given age.
    """
    def is_dead(self):
        return self.status == DEAD

class Septoria(Ring):
    """ Ring of Lesion of Septoria at a given age.
    """
    DEPOSIT = 0
    EMERGENT = 1
    INCUBATING = 2
    CHLOROTIC = 3
    SPORULATING = 4
    EMPTY = 5
    DEAD = 6

    def __init__(self, status = DEPOSIT):
        """ Initialize the lesion. 
        
        :Parameters:
          - `status` (int): code for the initial status of the new ring.

        """
        super(Septoria, self).__init__()
        self.status = status
        self.surface = 0.
        self.age = 0.
        self.age_dday = 0.
        self.age_dday_chlorotic = 0.
        self.cumul_wetness = 0.
        self.cumul_rain_event = 0.
        self.rain_before = False
    
    def can_form_new_ring(self, leaf, environment=None):
        """ Check if the lesion can form a new ring.
        
        A new ring is not created :
            * if the first du (which is modelled as a ring) is dead
            * if the infection is not achieved
            * during the time step when the infection occurs
            * if the surface max of the lesion has been reached
            * if there is no green surface on the leaf
            
        """
        lesions = leaf.lesions
        healthy_surface = leaf.healthy_surface
        incubating_surface = sum([les.surface for les in lesions if les.status == self.INCUBATING])
        green_surface = healthy_surface + incubating_surface
        
        lesion = environment['lesion']
        
        if( leaf.temp < lesion.fungus.basis_for_dday or
            lesion.status == self.DEAD or
            lesion.status == self.DEPOSIT or
            (lesion.status == self.EMERGENT and lesion.age_dday == 0.) or
            lesion.surface == lesion.fungus.Smax or
            green_surface == 0. ):
            return False
        else:
            return True
                
    def update(self, dt, leaf, environment=None):
        """ Update the status of the ring.
        
        * Cumulate the age of the ring.
        * Assign leaf data to the ring in order to access it in the methods.
        * Call the property 'stage' depending on the current status of the ring.
        
        :Parameters:
          - `dt` (float): delta time.
          - `ddday` (float): delta degree day (base temperature = -2°C).
          - `leaf` (class): a leaf sector with properties (e.g. healthy surface,
          senescence, rain intensity, wetness, temperature, lesions).
          
        """
        lesion = environment['lesion']
        
        if self.DEAD > self.status > self.DEPOSIT:
            self.age += dt
            ddday = (leaf.temp - lesion.fungus.basis_for_dday)/24
            self.age_dday += ddday
            # NOTE G.Garin 24/09/2012 :
            # Dividing by 24 because of hourly time step.
        if self.status == self.CHLOROTIC:
            self.age_dday_chlorotic += ddday
        self.leaf = leaf
        self.stage(dt=dt, ddday=(leaf.temp - lesion.fungus.basis_for_dday)/24, environment=environment)
        
    def du_on_leaf(self, environment=None, **kwds):
        """ Compute the success of infection by the deposited dispersal unit,
        and the probability of washing by the rain.
        
        The dispersal unit is considered as the potential future first 'ring' by the model.
        * Washing by the rain : 
          Function on the rain intensity and the healthy surface on the leaf sector.
          (See equation (2) in Robert et al., 2008)
          
        * Success of infection : 
          - Cumulate the CONTINUOUS wetness duration of the ring.
          - On the time step, the infection occurs only between terminals of 
          temperature, and if the wetness duration of the ring is sufficient
          (see fungus parameters : 'temp_min', 'temp_max', 'wd_min').
          - Otherwise, if the ring is not wet:
              ¤ if the temperature is not between the terminals, nothing happens.
              ¤ if the ring has never been wet, it can die under a certain 
              probability (see fungus parameters : 'loss_rate').
              ¤ if the ring has already been wet, it dies automatically.

        :Parameters:
          - `environment` (dict): data from the environment in which we use here :
              ¤ `lesion` (class): properties of the lesion (e.g. parameters for the septoria).
        
        :Returns:
          - `self.wetness` (float) : continuous wetness duration of the ring (in hours).
          - `self.status` (int) : code for the status of the ring.
          (here, can stay in DEPOSIT, can become INCUBATING or DEAD)
        """
       
        assert(self.status == self.DEPOSIT)
        
        leaf = self.leaf
        leaf_wet = leaf.wetness # (boolean): True if the leaf sector is wet during this time step.
        temp = leaf.temp # (float) : mean temperature on the leaf sector during the time step (in °C).
        rain_intensity = leaf.rain_intensity # (float) : rain intensity on the leaf sector during the time step (in mm/h).
        healthy_surface = leaf.healthy_surface # (float) : healthy surface (=with no lesion) on the leaf sector during the time step (in cm^2).
        
        lesion = environment['lesion']
        
        if leaf_wet:
            self.cumul_wetness += 1
        elif self.cumul_wetness > 0: # self.wetness > 0 ou un peu plus ?
            assert not leaf_wet
            self.cumul_wetness = 0
            self.status = self.DEAD # The dispersal unit dies
        else:
            assert not leaf_wet
            assert self.cumul_wetness == 0
        
        # TODO: design a new equation : coherente (adimesionnelle)
        # et valable pour differents pas de temps...
        # A mettre ailleurs !
        washing_rate = rain_intensity / (healthy_surface + rain_intensity)
        
        if proba(washing_rate):
            self.status = self.DEAD # The dispersal unit dies.
        
        # Notion de temps hydro-thermique, cumul de temperature.
        if (lesion.fungus.temp_min <= temp <= lesion.fungus.temp_max) and self.cumul_wetness >= lesion.fungus.wd_min :
            self.status = self.EMERGENT # The ring emerges.
            # Ajouter proba inf.
        elif self.cumul_wetness == 0 :
            if proba(lesion.fungus.loss_rate): # TODO : Proba conditionnelle doit se cumuler.
                self.status = self.DEAD # The dispersal unit dies. 
    
    def emergent(self, ddday, environment=None, **kwds):
        """ Compute the surface of the ring according to the status of the lesion,
        and assign the adequate status to the ring.
        
        Once it has emerged, the ring surface remains constant all its lifespan. Its status evolves.
        * For the first ring : its surface is always 'epsilon' (see fungus parameters).
        * For the other rings : their surface is a function of the thermal time, and 
        of the presence of other lesions on the leaf sector.
        
        At the time of their emergence, the first rings are INCUBATING. Once the incubation
        of the lesion is over, the new rings emerge as CHLOROTIC.
          
         :Parameters:
          - `ddday` (float): delta degree day (base temperature = -2°C).
          - `environment` (dict): data from the environment in which we use here :
              ¤ `lesion` (class): properties of the lesion (e.g. status, parameters for the septoria).
          
        :Returns:
          - `self.surface`(float): surface of the freshly emerged ring (in cm^2).
          - `self.status` (int): code for the status of the ring
          (here can change to INCUBATING, or CHLOROTIC).       
        """     
        assert(self.status == self.EMERGENT)
        
        leaf = self.leaf
        lesions = leaf.lesions # (class) : Lesions on the leaf sector with their properties.
        nb_lesions = len(lesions) # (int) : number of lesions on the leaf sector.
        nb_incubating_lesions = len([les for les in lesions if les.status == self.INCUBATING])
        # (int) : number of lesions in incubation on the leaf sector.
        incubating_surface = sum([les.surface for les in lesions if les.status == self.INCUBATING])
        # (int) : surface of the lesions in incubation on the leaf sector (in cm^2).
        healthy_surface = leaf.healthy_surface
        # (int) : surface with no lesion on the leaf sector (in cm^2).
        
        lesion = environment['lesion']
        
        if healthy_surface > 0:
            if lesion.status == self.EMERGENT :
                self.surface = lesion.fungus.epsilon
                self.status = self.INCUBATING # The ring begins the incubation.
            elif lesion.status == self.INCUBATING :
                free_space = healthy_surface / nb_incubating_lesions # An incubating ring can emerge only on an healthy surface.
                self.surface = min(free_space, lesion.fungus.Smin * ddday / lesion.fungus.degree_days_to_chlorosis)
                self.status = self.INCUBATING # The ring appears in incubation.
            else:
                free_space = (healthy_surface + incubating_surface)/ (nb_lesions - nb_incubating_lesions)
                # A chlorotic ring can emerge only on green surface (= healthy surface + incubating surface)
                size_before_Smax = lesion.fungus.Smax - lesion.surface
                self.surface = min(free_space, size_before_Smax, lesion.fungus.growth_rate * ddday)
                self.status = self.CHLOROTIC # The ring is directly chlorotic.
            
            # TODO : The chlorotic lesions must grow before the incubating ones.
        
    def incubating(self, environment=None, **kwds):
        """ Set the status of the ring to CHLOROTIC when needed.
        
        An incubating ring becomes automatically chlorotic when the status
        of the very first ring (= status of the lesion) becomes chlorotic.
        It is a function of the thermal age.
        
        :Parameters:
          - `environment` (dict): data from the environment in which we use here :
              ¤ `lesion` (class): properties of the lesion (e.g. status, parameters for the septoria).
        
        :Returns:
          - `self.status` (int): code for the status of the ring
            (can stay INCUBATING or can change to CHLOROTIC)
        """
        assert(self.status == self.INCUBATING)
        
        lesion = environment['lesion']
        
        if self == lesion.rings[0]: # Particular case of the first ring.
            if self.age_dday >= lesion.fungus.degree_days_to_chlorosis:
                self.status = self.CHLOROTIC
                # TODO maybe ou TOTHINK : calcul du growth rate exact au prorata
        else : # Case of the following rings.
            if lesion.status == self.CHLOROTIC:
                self.status = self.CHLOROTIC # The ring becomes chlorotic.
     
    def chlorotic(self, environment=None, **kwds):
        """ Set the status of the ring to SPORULATING when needed.
        
        Each ring entering in the CHLOROTIC stage must wait 110 DD 
        to be SPORULATING.
        
        :Parameters:
          - `environment` (dict): data from the environment in which we use here :
              ¤ `lesion` (class): properties of the lesion (e.g. parameters for the septoria).
        
        :Returns:
          - `self.status` (int): code for the status of the ring
            (can stay CHLOROTIC or can change to SPORULATING)
        """
        assert(self.status == self.CHLOROTIC)
        
        lesion = environment['lesion']
        
        if self.age_dday_chlorotic >= lesion.fungus.degree_days_to_sporulation:
            self.status = self.SPORULATING # The ring begins the sporulation.

    def sporulating(self, environment=None, **kwds):
        """ Compute the number of rain events on the ring, 
        and update the status of the ring when needed.
        
        A sporulating ring bears fructifications containing dispersal units.
        These dispersal units are spread by the rain if the relative humidity
        is greater than or equal to 85%. It is assumed that the ring is EMPTY
        after 3 separate rain events.
        
        :Parameters:
          - `environment` (dict): data from the environment in which we use here :
              ¤ `lesion` (class): properties of the lesion (e.g. status, parameters for the septoria).
        
        :Returns:
          - self.status : developement stage of the ring 
            (here stays in sporulation or dies)
          - newdu : number
        """
        assert(self.status == self.SPORULATING)
        
        leaf = self.leaf
        rain_intensity = leaf.rain_intensity # (float) : rain intensity on the leaf sector during the time step (in mm/h).
        relative_humidity = leaf.relative_humidity # (float) : relative humidity on the leaf sector during the time step (in %).
        lesion = environment['lesion']

        if rain_intensity > 0 and relative_humidity >= lesion.fungus.rh_min and not self.rain_before:
            self.cumul_rain_event += 1
        
        if self.cumul_rain_event >= lesion.fungus.rain_events_to_empty:
            self.status = self.EMPTY # The ring is empty.
            
        if rain_intensity == 0:
            self.rain_before = False
        else:
            self.rain_before = True
        
    def empty(self, *args, **kwds):
        """ Assert if the status of the ring is 'EMPTY'.
        """
        assert(self.status == self.EMPTY)
    
    def dead(self, *args, **kwds):
        """ Assert if the status of the ring is 'DEAD'.

        """
        assert(self.status == self.DEAD)
        
    @property
    def stage(self):
        if self.status == self.DEPOSIT:
            return self.du_on_leaf
        elif self.status == self.EMERGENT:
            return self.emergent
        elif self.status == self.INCUBATING:
            return self.incubating
        elif self.status == self.CHLOROTIC:
            return self.chlorotic
        elif self.status == self.SPORULATING:
            return self.sporulating
        elif self.status == self.EMPTY:
            return self.empty
        elif self.status == self.DEAD:
            return self.dead
        else:
            return

        
class PowderyMildew(Ring):
    """ Ring of Lesion of PowderyMildew at a given age.
    """
    DEPOSIT = 0
    EMERGENT = 1
    LATENT = 2
    SPORULATING = 3
    EMPTY = 4
    DEAD = 5

    def __init__(self, status = DEPOSIT):
        """ Initialize the lesion. 
        
        :Parameters:
          - `status` (int): code for the initial status of the new ring.

        """
        super(PowderyMildew, self).__init__()
        self.status = status
        self.surface = 0.
        self.age = 0.
        self.latency_progress = 0.
        self.sporulation_progress = 0.
        self.nb_du = 0.
        
    def can_form_new_ring(self, leaf, environment=None):
        """ Check if the lesion can form a new ring.
        
        A new ring is not created :
            * if the first du (which is modelled as a ring) is dead
            * if the infection is not achieved
            * during the time step when the infection occurs
            * if the surface max of the lesion has been reached
            * if there is no green surface on the leaf
            
        """
        lesions = leaf.lesions
        healthy_surface = leaf.healthy_surface
        
        lesion = environment['lesion']
        
        if (lesion.status == self.DEAD or
            lesion.status == self.DEPOSIT or 
            (lesion.status == self.EMERGENT and lesion.age_dday == 0.) or
            healthy_surface == 0. ):
            return False
        else:
            return True
        
        # TODO : 
        #   - Pas de Smax pour oidium --> A ameliorer
        #   - Mieux si cette fonction depend pas du type de champignon.
        
    def update(self, dt, ddday, leaf, environment=None):
        """ Update the status of the ring.
        
        * Cumulate the age of the ring.
        * Assign leaf data to the ring in order to access it in the methods.
        * Call the property 'stage' depending on the current status of the ring.
        
        :Parameters:
          - `dt` (float): delta time.
          - `leaf` (class): a leaf sector with properties (e.g. healthy surface,
          senescence, rain intensity, wetness, temperature, lesions).
          
        """
        if self.DEAD > self.status > self.DEPOSIT:
            self.age += dt
        self.leaf = leaf
        self.stage(dt=dt, environment=environment)
        
    def du_on_leaf(self, environment=None, **kwds):
        """ Compute the success of infection by the deposited dispersal unit.
        
        The success of infection is a function of the temperature, the relative
        humidity and the wetness of the leaf during the time step.
        
        :Parameters:
          - 
        
        :Returns:
          - 
        """
        assert(self.status == self.DEPOSIT)
        
        # External variables
        leaf = self.leaf
        leaf_wet = leaf.wetness
        temp = leaf.temp
        relative_humidity = leaf.relative_humidity
        lesion = environment['lesion']
        
        # Raw parameters for the calculation
        temp_min = lesion.fungus.temp_min_for_infection
        temp_max = lesion.fungus.temp_max_for_infection
        m = lesion.fungus.m_for_infection
        n = lesion.fungus.n_for_infection
        max_infection_rate = lesion.fungus.max_infection_rate
        decay_rate = lesion.fungus.decay_rate
        a_RH_effect = lesion.fungus.a_RH_effect
        b_RH_effect = lesion.fungus.b_RH_effect
        RH_opt = lesion.fungus.RH_opt_for_infection
        c_wetness_effect = lesion.fungus.c_wetness_effect
        d_wetness_effect = lesion.fungus.d_wetness_effect
       
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
        
        # Infection rate 
        infection_rate = temp_factor * RH_factor * wetness_factor
        if proba(infection_rate):
            self.status = self.EMERGENT
        else:
            self.status = self.DEAD
        
        # TODO : A ameliorer :
        # - Gerer le delai avant l'infection --> Changement d'echelle temporelle
        # - Duree de vie d'une spore si elle n'a pas fait d'infection --> Introduction
        # d'un loss_rate.
           
    def emergent(self, environment=None, **kwds):
        """ Compute the surface of the ring and assign the adequate status to it.
        
                
         :Parameters:
          - 
          
        :Returns:
          - 
        """     
        assert(self.status == self.EMERGENT)
        
        # External variables
        leaf = self.leaf
        temp = leaf.temp
        healthy_surface = leaf.healthy_surface
        lesions = leaf.lesions
        nb_lesions = len(lesions)
        lesion = environment['lesion']
        
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
        self.surface = min(free_space, 
                           (pi/4) * (kmax * temp_norm_function_for_growth * 
                           growth_rate * exp(growth_rate * (half_growth_time - lesion.age)) /
                           ((1. + exp(growth_rate * (half_growth_time - lesion.age)))**2))**2)
        
        self.status = self.LATENT
        
    def latent(self, environment=None, **kwds):
        """ Set the status of the ring to SPORULATING when needed.
               
        :Parameters:
          - 
        
        :Returns:
          - 
        """
        assert(self.status == self.LATENT)
        
        # External variables
        leaf = self.leaf
        temp = leaf.temp
        lesion = environment['lesion']
       
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
            self.status = self.SPORULATING
     
    def sporulating(self, environment=None, **kwds):
        """ Compute the number of rain events on the ring, 
        and update the status of the ring when needed.
               
        :Parameters:
          -
        
        :Returns:
          - 
        """
        assert(self.status == self.SPORULATING)
        
        # External variables
        dispersal_rate = environment['dispersal_rate'] 
        leaf = self.leaf
        temp = leaf.temp
        lesion = environment['lesion']
        
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
            self.status = self.EMPTY
        
        # TODO : Calculs a revoir... + Integrer stochasticite.
        
    def empty(self, *args, **kwds):
        """ Assert if the status of the ring is 'EMPTY'.
        """
        assert(self.status == self.EMPTY)
    
    def dead(self, *args, **kwds):
        """ Assert if the status of the ring is 'DEAD'.

        """
        assert(self.status == self.DEAD)
        
    @property
    def stage(self):
        if self.status == self.DEPOSIT:
            return self.du_on_leaf
        elif self.status == self.EMERGENT:
            return self.emergent
        elif self.status == self.LATENT:
            return self.latent
        elif self.status == self.SPORULATING:
            return self.sporulating
        elif self.status == self.EMPTY:
            return self.empty
        elif self.status == self.DEAD:
            return self.dead
        else:
            return

def temp_norm_function(temp_mean, temp_min, temp_max, m, n):
    """ Compute the normalized temperature function.
    
    Compute the normalized temperature function as in Calonnec et al., 2008
    
    :Parameters:
      - `temp_mean` (float): mean temperature during the time step (°C).
      - `temp_min` (float): minimum temperature during the time step (°C).
      - `temp_max` (float): maximum temperature during the time step (°C).
      - `m` (float): shape parameter.
      - `n` (float): shape parameter.
    
    :Returns:
      - `temp_norm_function` : normalized temperature function.
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
                 basis_for_dday = -2,
                 temp_min = 10,
                 temp_max = 30,
                 wd_min = 10,
                 loss_rate = 1./120,
                 degree_days_to_chlorosis = 220,
                 degree_days_to_sporulation = 110,
                 epsilon = 0.001,
                 Smin = 0.03,
                 Smax = 0.3,
                 growth_rate = 0.0006,
                 rh_min = 85,
                 rain_events_to_empty = 3,            
                 *args, **kwds):
        """ Parameters for septoria.
        
        :Parameters:
            - `basis_for_dday` (float): basis temperature for the accumulation of degree days (°C)
            - `temp_min` (float): Minimal temperature for infection.
            - `temp_max` (float): Maximal temperature for infection.
            - `wd_min` (float): Minimal wetness duration for infection.
            - `loss_rate` (float): Loss rate of dispersal units in 1 hour.
            - `degree_days_to_chlorosis` (float): Thermal time between emergence 
            and chlorosis (i.e. incubation for the first rings).
            - `degree_days_to_sporulation` (float): Thermal time between chlorosis 
            and sporulation.
            - `epsilon` (float): Initial size of incubating lesion (cm2).
            - `Smin` (float): Initial size of chlorotic lesion (cm2).
            - `Smax` (float): Lesion maximum size (cm2).
            - `growth_rate` (float): Lesion growth rate (cm2.dday-1).
            - `rh_min` (float): Minimal relative humidity for sporulation.
            - `rain_events_to_empty` (int): number of rain events to empty a sporulating ring.
            
        """
        self.name = "Septoria"
        self.basis_for_dday = basis_for_dday
        self.temp_min = temp_min
        self.temp_max = temp_max
        self.wd_min = wd_min
        self.loss_rate = loss_rate
        self.degree_days_to_chlorosis = degree_days_to_chlorosis
        self.degree_days_to_sporulation = degree_days_to_sporulation
        self.epsilon = epsilon
        self.Smin = Smin
        self.Smax = Smax
        self.growth_rate = growth_rate
        self.rh_min = rh_min
        self.rain_events_to_empty = rain_events_to_empty

def septoria(**kwds):
    return SeptoriaParameters(**kwds)
    
class PowderyMildewParameters(Parameters):
    def __init__(self,
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
        
        :Parameters:
            - `temp_min_for_infection` (float): minimal temperature for infection (°C).
            - `temp_max_for_infection` (float): maximal temperature for infection (°C).
            - `m_for_infection` (float): shape parameter for the calculation of the
               normalized temperature for infection.
            - `n_for_infection` (float): shape parameter for the calculation of the
               normalized temperature for infection.
            - `max_infection_rate` (float): maximal infection rate ([0;1]).
            - `decay_rate` (float): leaf susceptibility decay with age ([0;1]).
            - `a_RH_effect` (float): shape parameter for the calculation of the effect
            of relative humidity on infection.
            - `b_RH_effect` (float): shape parameter for the calculation of the effect
            of relative humidity on infection.
            - `RH_opt_for_infection` (float): optimal relative humidity for infection (%).
            - `c_wetness_effect` (float): shape parameter for the calculation of the effect
            of leaf wetness on infection.
            - `d_wetness_effect` (float): shape parameter for the calculation of the effect
            of leaf wetness on infection.
            - `diameter_max` (float): maximum diameter of the lesion (cm).
            - `diameter_min` (float): minimum diameter of the lesion (cm).
            - `leaf_age_effect` (float): rate of colony growth with the leaves age ([0;1]).
            - `temp_min_for_growth` (float): minimal temperature for growth (°C).
            - `temp_max_for_growth` (float): maximal temperature for growth (°C).
            - `m_for_growth` (float): shape parameter for the calculation of the
               normalized temperature for growth.
            - `n_for_growth` (float): shape parameter for the calculation of the
               normalized temperature for growth.
            - `half_growth_time` (float): time before 50% of colony growth (hours).
            - `temp_min_for_latency` (float): minimal temperature for latency (°C).
            - `temp_max_for_latency` (float): maximal temperature for latency (°C).
            - `m_for_latency` (float): shape parameter for the calculation of the
               normalized temperature for latency.
            - `n_for_latency` (float): shape parameter for the calculation of the
               normalized temperature for latency.
            - `min_latency_duration` (float): minimum latency duration (hours).
            - `a_for_sporulation` (float) : shape parameter for the calculation of the
               progress of the sporulation period according to temperature.
            - `b_for_sporulation` (float) : shape parameter for the calculation of the
               progress of the sporulation period according to temperature.                 temp_min_for_sporulation
            - `temp_min_for_sporulation` (float) : minimal temperature for sporulation (°C).
            - `temp_max_for_sporulation` (float) : maximal temperature for sporulation (°C).
            - `m_for_sporulation` (float): shape parameter for the calculation of the
               normalized temperature for sporulation.
            - `n_for_sporulation` (float) : shape parameter for the calculation of the
               normalized temperature for sporulation.
            - `nb_du_max` (int) : number of dispersal unit produced by a sporulating
            surface unit by time step (dispersal unit / cm2 / hour).
            - `treshold_nb_du_to_empty` (int) : treshold of dispersal unit on a surface to 
            consider it empty.

        """
        self.name = "PowderyMildew"
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

def powdery_mildew(**kwds):
    return PowderyMildewParameters(**kwds)