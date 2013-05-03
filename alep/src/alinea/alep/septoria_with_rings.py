""" Class of lesion specific of wheat septoria, with rings whose surface grow then 
    stay the same when ageing.

"""
# Imports #########################################################################
from alinea.alep.fungal_objects import *
from alinea.alep.septoria import *
from random import randint, seed
import numpy as np
seed(1)

# Lesion ##########################################################################
class SeptoriaWithRings(Lesion):
    """ Define a lesion specific of septoria with growth rings.
    """
    def __init__(self, nb_spores=None, position=None):
        """ Initialize the lesion of septoria. 
        
        Parameters
        ----------
        position: non defined
            Position of the dispersal unit on the phyto-element
        nb_spores: int
            Number of spores aggregated in the dispersal unit

        """
        super(SeptoriaWithRings, self).__init__(nb_spores=nb_spores, position=position)
        # List of rings, add a ring at lesion formation
        self.rings = []
        ring = SeptoriaRing(lesion=self, status=self.fungus.INCUBATING)
        self.rings.append(ring)
        # Age of the center of the lesion (degree days)
        self.age_dday = 0.
        # Delta degree days during time step
        self.ddday = 0.
        # Position of senescence the time step before (Useful in case of senescence
        # to compute the time left for growth before senescence occur)
        self.old_position_senescence = None
        # Degree days before senescence if 'self.is_senescent'
        self.ddday_before_senescence = None       
        # Surface of disabled rings (empty or senescent)
        self.surface_dead = 0.
        # Surface of empty rings
        self.surface_empty = 0.
        # Growth activity of the lesion
        self.growth_is_active = True
        # Growth demand of the lesion (cm2)
        self._growth_demand = None
        # Is first hour of rain
        self.first_rain_hour = False
    
    def update(self, dt, leaf, **kwds):
        """ Update the status of the lesion and create a new growth ring if needed.
                
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)

        """
        assert self.is_active
        f = self.fungus
        
        # Compute delta degree days in dt
        # TODO : modify if list of temp coming from weather data
        # self.compute_delta_ddays(dt, leaf)
        self.compute_delta_ddays_from_weather(leaf)
        ddday = self.ddday
        
        # If senescence, compute length of growth period before senescence during time step
        if self.is_senescent:
            self.compute_time_before_senescence(ddday=ddday, leaf=leaf)
            ddday = self.ddday_before_senescence

        # Update the age of the lesion
        self.age_dday += ddday
            
        # Ageing of the rings / create new ones when needed
        nb_rings_initial = len(self.rings)
        # Note that new rings can be added in this time step, 
        # their update is managed by another module
        for i in range(nb_rings_initial):
        # for ring in self.rings:
            self.rings[i].update(ddday=ddday, lesion=self)
        
        # Update the perception of rain by the lesion
        if self.status == f.SPORULATING:
            if leaf.rain_intensity > 0. and leaf.relative_humidity >= f.rh_min:
                self.first_rain_hour = True if not self.first_rain_hour else False
            else:
                self.first_rain_hour = False
        
    def add_rings(self, list_of_rings=None):
        """ Add rings to list of rings if needed.
        
        Parameters
        ----------
        list_of_rings: list of ring instantiations
            Rings of lesion with properties (e.g. surface, status, age, etc.)       
        """
        self.rings += list_of_rings
    
    def compute_delta_ddays(self, dt=1., leaf=None):
        """ Compute delta degree days in dt.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions) 
        """
        f = self.fungus
        # Calculation
        if dt != 0.:
            ddday = max(0,(leaf.temp - f.basis_for_dday)/(24./dt))
        else:
            ddday = 0.
        # Save variable
        self.ddday = ddday 

    def compute_delta_ddays_from_weather(self, leaf=None):
        """ Compute delta degree days from weather data since last call.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions) 
        """
        f = self.fungus
        temp_list = np.array(leaf.temp_list)
        dt = len(temp_list)
        # Calculation
        if dt != 0.:
            ddday = max(0, sum((temp_list - f.basis_for_dday))/(24./dt))
        else:
            ddday = 0.
        # Save variable
        self.ddday = ddday
     
    def control_growth(self, growth_offer = 0.):
        """ Reduce surface of the last ring up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Surface available on the leaf for the ring to grow (cm2)
        """
        f = self.fungus
        total_growth_demand = self.growth_demand
        growth_offer_left = min(growth_offer, f.Smax - self.surface)
               
        rings_in_formation = [r for r in self.rings if r.is_in_formation(fungus=f)]
        if rings_in_formation:
            for ring in rings_in_formation:
                ring_growth_offer = min(ring.growth_demand, growth_offer_left)
                growth_offer_left -= ring.growth_demand if growth_offer_left>0. else 0.
                ring.control_growth(ring_growth_offer, lesion=self)
        
        # Update stock_spores
        new_rings_sporulating = [r for r in self.rings if (r.is_sporulating(fungus=f) and r.stock_spores==None)]
        # Note: at initiation, stock_spores=-1. It is only filled once.
        # The aim is to avoid refilling it when it falls back to 0.
        if new_rings_sporulating:
            for ring in new_rings_sporulating:
                ring.update_stock(lesion=self)
        
        # Disable growth when needed
        if growth_offer < total_growth_demand or self.surface == f.Smax:
            self.disable_growth()
            for r in self.rings:
                r.disable_growth()

        # Update surface dead
        self.update_surface_dead()
                
    def update_surface_dead(self):
        """ Update surface dead of the lesion and disable it if no ring active.
        
        Parameters
        ----------
            None
        """
        # Remove non active rings
        self.surface_dead += sum(ring.surface for ring in self.rings if not ring.is_active)
        self.surface_empty += sum(ring.surface for ring in self.rings if ring.is_empty(self.fungus))
        self.rings = [ring for ring in self.rings if ring.is_active]

        # Disable lesion when no ring left
        if not self.rings:
            self.disable()
        
    def compute_time_before_senescence(self, ddday=0., leaf=None):
        """ Compute length of growth period before senescence during time step.
        
        Parameters
        ----------
        ddday: float
            Delta degree days during 'dt'
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        old_pos = self.old_position_senescence
        new_pos = leaf.position_senescence
        speed = abs(new_pos-old_pos)/ddday if ddday > 0. else 0.
        
        self.ddday_before_senescence = abs(self.position[0]-old_pos)/speed if speed >0. else 0.
    
    def become_senescent(self, old_position_senescence=0.):
        """ The lesion will become senescent during this time step 
            and save the position of senescence the time step before.
            
        Parameters
        ----------
        old_position_senescence: float
            Position of senescence the time step before
        """
        # Turn on 'is_senescent'
        self.is_senescent = True
        # Save position of senescence the time step before.
        self.old_position_senescence = old_position_senescence
        
    def senescence_response(self, time_since_sen=0.):
        """ Kill rings affected by senescence and achieve ageing of the other rings
            up to the end of time step.
        
        During the time step, the growth of the lesion has been stopped when 
        senescence occured. Now, at the end of time step, we can suppress all
        the rings that do not survive senescence, and complete the ageing of the 
        others.
        
        Parameters
        ----------
        time_since_sen: float
            Time since growth interruption (degree days)
        """
        f = self.fungus
        
        dday_since_senescence = self.ddday - self.ddday_before_senescence 
        self.age_dday += dday_since_senescence
        
        # Stop growth
        self.disable_growth()
        
        for ring in self.rings:
            if ring.status < f.NECROTIC:
                # Disable all rings under NECROTIC state:
                ring.disable()
            else:
                # To be sure, reset all useless ring variables
                # Then update them to make them reach the end of time step
                ring.disable_growth()
                ring.update(ddday=dday_since_senescence, lesion=self)
        
        # Update surface dead
        self.update_surface_dead()
    
    def emission(self, leaf, **kwds):
        """ Create a list of dispersal units emitted by the entire lesion.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        
        Returns
        -------
        emissions: 
        """
        emissions = []
        for ring in self.rings:
            if ring.is_stock_available(lesion = self):
                ring_emissions = ring.emission(leaf, lesion=self)
                if ring_emissions:
                    emissions += ring_emissions
        return emissions    
        
    def is_dead(self):
        """ Update the status of all the rings to 'DEAD' if the lesion is dead.
        
        Parameters
        ----------
            None
        """
        return all(ring.is_dead() for ring in self.rings)
        
    def compute_all_surfaces(self):
        """ Compute lesion surfaces in different states.
        
        Parameters
        ----------
            None
            
        ..NOTE:: Temporary just for tests
        """
        f = self.fungus
        self.surface_inc = sum(r.surface for r in self.rings if r.is_incubating(f))
        self.surface_chlo = sum(r.surface for r in self.rings if r.is_chlorotic(f))
        self.surface_nec = sum(r.surface for r in self.rings if r.is_necrotic(f))
        self.surface_spo = sum(r.surface for r in self.rings if r.is_sporulating(f))+self.surface_empty
      
    @property
    def surface_alive(self):
        """ Compute the surface alive on the lesion.
        
        Parameters
        ----------
            None
            
        ..NOTE:: Temporary just for tests
        """
        return sum(ring.surface for ring in self.rings)
  
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
        """ Compute the growth_demand of the lesion.
        
        If a single ring is in formation, growth demand of the lesion equals growth 
        demand of this ring. If a new ring is emerging at this time step, two rings
        are in formation.        
        
        Parameters
        ----------
            None
            
        Returns
        -------
        surface: float
            Surface of the whole lesion (cm2)
        """
        fungus = self.fungus
        
        if self.rings:
            forming_rings = [r for r in self.rings if r.is_in_formation(fungus)]
            self._growth_demand = sum(r.growth_demand for r in forming_rings)
        else:
            self._growth_demand = 0.
        return self._growth_demand
    
    @property
    def stock_spores(self):
        """ Compute the stock of spores on the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        surface: float
            Surface of the whole lesion (cm2)
        """
        stock_spores = sum(ring.stock_spores for ring in self.rings if ring.stock_spores >0.)
        return stock_spores
    
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
            
    # @property
    # def age_dday(self):
        # """ Compute the thermal age of the lesion.
        
        # Parameters
        # ----------
            # None
            
        # Returns
        # -------
        # age_dday: float
            # Age of the lesion in degree days
        # """
        # if self.rings:
            # return self.rings[0].age_dday

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
            
    @growth_demand.setter
    def growth_demand(self, value):
        """ Set the growth demand of the lesion to the chosen value.
        
        Parameters
        ----------
        value : int
            Chosen value to set lesion growth demand
            
        Returns
        -------
            None
        """
        self._growth_demand = value
        
# Rings ###########################################################################
class SeptoriaRing(Ring):
    """ Ring of Lesion of Septoria at a given age.
    """
    def __init__(self, lesion, status):
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
        f = lesion.fungus        
        # Status of the ring
        self.status = status
        # Surface of the ring
        self.surface = 0.
        # Age of the ring
        self.age_dday = 0.
        # Delta degree days of growth during time step
        # self.ddday = 0.
        self.delta_growth = 0.
        # Delta age to complete the ring
        if len(lesion.rings)==0.:
            self.delta_age_ring = f.degree_days_to_chlorosis + f.delta_age_ring
        else:
            self.delta_age_ring = f.delta_age_ring
        # Activity of the ring
        self.is_active = True
        # Growth activity of the ring
        self.growth_is_active = True
        
        # See later for dispersion
        self.cumul_rain_event = 0.
        self.rain_before = False
        self.stock_spores = None

    def is_in_formation(self, fungus):
        """ Can keep growing!!! """
        ok = self.growth_is_active
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

    def update(self, ddday, lesion=None):
        """ Update the status of the ring.
        
        * Cumulate the age of the ring.
        * Assign leaf data to the ring in order to access it in the methods.
        * Call the property 'stage' depending on the current status of the ring.
        
        Parameters
        ----------
        ddday: float
            Delta degree days during 'dt'
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)            
        """
        f = lesion.fungus
        # Ageing of the ring
        self.age_dday += ddday
        if not self.is_in_formation(fungus=f):
            # Compute status of the ring
            self.stage(lesion=lesion)
        else:
            # Create new rings if needed
            if ddday > self.delta_age_ring:
                self.create_new_rings(ddday=(ddday - self.delta_age_ring), lesion=lesion)
                # Update delta age of growth
                self.delta_growth = self.delta_age_ring
                # Reset delta_age_ring
                self.delta_age_ring = None
            else:
                # Update delta age of growth
                self.delta_growth = ddday
                # Update delta age left until growth completion
                self.delta_age_ring -= ddday
            
            # Compute status of the ring
            self.stage(lesion=lesion)
                
    def create_new_rings(self, ddday=0., lesion=None):
        """ Add rings to lesion list of rings if needed.
        
        Parameters
        ----------
        ddday: float
            Delta degree days in 'dt' since apparition of new rings
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)   
        """
        f = lesion.fungus
        list_of_rings = []

        # Creation of the first new ring
        new_ring = SeptoriaRing(lesion=lesion, status=f.CHLOROTIC)
        list_of_rings.append(new_ring)
        list_of_rings[-1].age_dday = ddday 
        list_of_rings[-1].delta_growth = min(ddday, list_of_rings[-1].delta_age_ring)
        
        # Creation of the following rings with properties ('age_dday', 'delta_growth', 'delta_age_ring')
        while list_of_rings[-1].age_dday > list_of_rings[-1].delta_age_ring:
            ddday -= list_of_rings[-1].delta_age_ring
            list_of_rings[-1].delta_age_ring = None
            new_ring = SeptoriaRing(lesion=lesion, status=f.CHLOROTIC)
            list_of_rings.append(new_ring)
            list_of_rings[-1].age_dday = ddday
            list_of_rings[-1].delta_growth = min(ddday, list_of_rings[-1].delta_age_ring)
            
        list_of_rings[-1].delta_age_ring -= list_of_rings[-1].age_dday
                
        # Status of new rings
        for ring in list_of_rings:
            # ring.delta_age_ring -= ring.age_dday 
            ring.stage(lesion=lesion)
        
        # Attach each new ring to the lesion
        lesion.add_rings(list_of_rings)
            
    def in_formation(self, lesion=None, **kwds):
        """ Compute growth demand according to delta age of growth during time step.        
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
            
        Returns
        -------
            None
        """     
        f = lesion.fungus
        delta_growth = self.delta_growth
       
        if self.is_incubating(fungus=f):
            # First ring in incubation
            self.growth_demand = f.Smin * delta_growth / f.degree_days_to_chlorosis
        else:
            # All the other rings
            self.growth_demand = delta_growth * f.growth_rate

    def control_growth(self, growth_offer = 0., lesion=None):
        """ Reduce surface of the rings in formation down to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Minimum between the surface available on the leaf for the ring
            to grow (cm2) and 'growth_demand'
        """
        f = lesion.fungus

        # Add surface
        self.surface += growth_offer

        if self.surface==0. and growth_offer == 0.:
            # Turn off the ring
            self.disable()
        elif self.delta_age_ring == None:
            # Turn off growth activity if no delta age left before growth completion
            self.disable_growth()
        
    def incubating(self, lesion=None, **kwds):
        """ Set the status of the ring to CHLOROTIC when needed.
        
        Only the first ring can be incubating. It must wait 220DD to become CHLOROTIC
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        """
        f = lesion.fungus
        time_to_chlorosis = f.degree_days_to_chlorosis
        assert self.is_incubating(fungus=f)
        
        # Compute status transition to necrotic
        if self.age_dday >= time_to_chlorosis:
            self.status = f.CHLOROTIC
            self.chlorotic(lesion=lesion)
    
    def chlorotic(self, lesion=None, **kwds):
        """ Set the status of the ring to NECROTIC when needed.
        
        Each ring entering in the CHLOROTIC stage must wait 110 DD 
        to be NECROTIC.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        """
        f = lesion.fungus
        if lesion.surface_dead==0. and self == lesion.rings[0]:
            time_to_necrosis = f.degree_days_to_chlorosis + f.degree_days_to_necrosis
        else:
            time_to_necrosis = f.degree_days_to_necrosis
        assert self.is_chlorotic(fungus=f)

        # Compute status transition to necrotic
        if self.age_dday >= time_to_necrosis:
            self.status = f.NECROTIC
            self.necrotic(lesion=lesion)
            
    def necrotic(self, lesion=None, **kwds):
        """ Set the status of the ring to SPORULATING when needed.
        
        Each ring entering in the CHLOROTIC stage must wait ??? DD 
        to be SPORULATING.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        """
        f = lesion.fungus
        if lesion.surface_dead==0. and self == lesion.rings[0]:
            time_to_spo = (f.degree_days_to_chlorosis + 
                           f.degree_days_to_necrosis + 
                           f.degree_days_to_sporulation)
        else:
            time_to_spo = (f.degree_days_to_necrosis + 
                           f.degree_days_to_sporulation)
        assert self.is_necrotic(fungus=f)
        
        # Compute status transition to sporulating
        if self.age_dday >= time_to_spo:
            self.status = f.SPORULATING

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

        .. Todo:: Enhance the code to parametrize the code with a function.
        """
        f = lesion.fungus
        assert self.is_sporulating(fungus=f)
        
        # Count dispersal events
        # if lesion.first_rain_hour:
            # self.cumul_rain_event += 1
        
        # Empty the ring after 3 rain events
        # if self.cumul_rain_event >= f.rain_events_to_empty:
            # self.empty(lesion)
            
        # self.update_stock(lesion=lesion)
        
    def update_stock(self, lesion=None):
        """ Update the stock of spores on the ring.
        
        Fill the stock of spores according to production rate and surface of the ring.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        """
        f = lesion.fungus

        production = self.surface * f.production_rate           
        nb_spores_produced = int(round(production))
        self.stock_spores = nb_spores_produced
        # Note: explains the +1 above :
        # At initiation stock_spores=-1. It is only filled once.
        # The aim is to avoid refilling it when it falls back to 0.

    def emission(self, leaf=None, lesion=None, **kwds):
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
        f = lesion.fungus
        emissions = []
        stock_available = int(self.stock_spores*2/3.)
        
        # TODO : improve below
        nb_DU_emitted = int(leaf.rain_intensity * self.surface * 1000)
        nb_spores_by_DU = []
        for DU in range(nb_DU_emitted):
            if stock_available > 0.:
                nb_spores = min(randint(5,100), stock_available)
                nb_spores_by_DU.append(nb_spores)
                stock_available -= nb_spores
                # Update stock_spores
                self.stock_spores -= nb_spores
                
        # Get rid of DUs without spores
        nb_DU_emitted = len(nb_spores_by_DU)
        
        # Empty stock_spores
        if self.stock_spores < 100:
            self.empty(lesion)

        # Return emissions
        emissions = [SeptoriaDU(nb_spores = nb_spores_by_DU[i], status='emitted')
                                for i in range(nb_DU_emitted)]
        
        return emissions
           
    def is_stock_available(self, lesion=None):
        """ Check if the stock of DU can be emitted.
        
        DU is free for liberation if :
            - there are DUs in stock_du
            - ring is sporulating
            - relative humidity is above a treshold
            - it is the first hour of rain
            
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        
        Returns
        -------
        True or False :
            Availability of the stock
        """
        f = lesion.fungus
        
        if (self.stock_spores and self.status == f.SPORULATING and 
            lesion.first_rain_hour):
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
        """
        f = lesion.fungus
        self.status = f.EMPTY
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
        self.is_active = False
        
    def disable_growth(self):
        """ Shut down lesion growth activity (turn it to False).
        
        Parameters
        ----------
            None
        """
        self.growth_is_active = False
        self.growth_demand = 0.
        self.delta_growth = 0.
        self.delta_age_ring = None
        
    def stage(self, lesion=None):
        """ Send the ring toward the proper function according to its status.
        
        Parameters
        ----------
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
            
        Returns
        -------
        None
            Call the method corresponding to ring status  
        """
        f = lesion.fungus
        
        # Cumulation of surface
        if self.is_in_formation(fungus=f):
            self.in_formation(lesion=lesion)
        
        # Ageing
        if self.is_incubating(fungus=f):
            self.incubating(lesion=lesion)
        elif self.is_chlorotic(fungus=f):
            self.chlorotic(lesion=lesion)
        elif self.is_necrotic(fungus=f):
            self.necrotic(lesion=lesion)
        
        # Update stock
        # if self.is_sporulating(fungus=f):
            # NOTE : update stock moved in growth control
            # self.sporulating(lesion=lesion)
        elif self.is_empty(fungus=f):
            self.empty(lesion=lesion)
        elif self.is_dead(fungus=f):
            self.dead(lesion=lesion)