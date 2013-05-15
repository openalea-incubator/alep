""" Classes of dispersal unit, lesion and ring specific of wheat septoria.

"""
# Imports #########################################################################
from alinea.alep.fungal_objects import *
# The following import would provoke a circular reference
# "from alinea.alep.septoria import SeptoriaDU"
# --> Moved in the method 'SeptoriaExchangingRings.emission()'
from alinea.alep.septoria import Disease as _Disease, SeptoriaParameters as _SeptoriaParameters
from random import randint, seed
from math import floor, ceil
import numpy as np

seed(1)
            
# Lesion ##########################################################################
class SeptoriaExchangingRings(Lesion):
    """ Septoria Lesion implemented with rings that exchange surfaces. """

    def __init__(self, nb_spores=None, position=None):
        """ Initialize the lesion of septoria. 
        
        Parameters
        ----------
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        position: non defined
            Position of the dispersal unit on the phyto-element
        """
        super(SeptoriaExchangingRings, self).__init__(nb_spores=nb_spores, position=position)
        f = self.fungus
        # Particular case of first ring
        self.first_ring = FirstRing(lesion=self)        
        # Other rings in same array
        self.surface_rings = np.array([])
        # Surfaces in each state
        self.surface_inc = 0.
        self.surface_chlo = 0.
        self.surface_nec = 0.
        self.surface_spo = 0.
        # Surface alive on the lesion
        self.surface_alive = 0.
        # Surface of disabled rings
        self.surface_dead = 0.
        # Position of senescence the time step before (Useful in case of senescence
        # to compute the time left for growth before senescence occur)
        self.old_position_senescence = None
        # Degree days before senescence if 'self.is_senescent'
        self.ddday_before_senescence = None
        # Age to reach to create another ring
        self.delta_age_ring = 0.
        # Age of the center of the lesion (degree days)
        self.age_dday = 0.
        # Delta degree days during time step
        self.ddday = 0.
        # Growth activity of the lesion
        self.growth_is_active = True
        # Growth demand of the lesion (cm2)
        self.growth_demand = None
        # Is first hour of rain
        self.first_rain_hour = False
        # Stock of spores
        self.stock_spores = None

    def update(self, dt, leaf=None):
        """ Update the status of the lesion and create a new growth ring if needed.
                
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        f = self.fungus
        # Compute delta degree days in dt
        # self.compute_delta_ddays(dt, leaf)
        self.compute_delta_ddays_from_weather(leaf)
        ddday = self.ddday
        
        # If senescence, compute length of growth period before senescence during time step
        if self.is_senescent:
            self.compute_time_before_senescence(ddday=ddday, leaf=leaf)
            ddday = self.ddday_before_senescence

        # Update the age of the lesion
        self.age_dday += ddday
        
        # Compute growth demand
        self.update_growth_demand()
        
        # Manage first ring
        self.first_ring.update(lesion=self)
        # Manage the other rings
        if len(self.surface_rings)>0:
            if not self.growth_is_active:
                # If growth is over, suppress empty surfaces corresponding to young rings
                self.surface_rings = self.surface_rings[self.surface_rings.nonzero()]
            if self.age_dday-f.degree_days_to_chlorosis > self.delta_age_ring:
                # diff = self.delta_age_ring - (self.age_dday - f.degree_days_to_chlorosis)
                # nb_rings = ceil(float(diff)/f.delta_age_ring)
                self.surface_rings = np.append(self.surface_rings, 0.)
                self.delta_age_ring += f.delta_age_ring
            self.exchange_surfaces()

        # Update stock of spores
        if self.status==f.SPORULATING:
            self.update_stock()
            
            # Manage rain perception
            if leaf.rain_intensity > 0. and leaf.relative_humidity >= f.rh_min:
                self.first_rain_hour = True if not self.first_rain_hour else False
            else:
                self.first_rain_hour = False
            
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
            ddday = max(0, sum((temp_list - f.basis_for_dday))/24.)
        else:
            ddday = 0.
        # Save variable
        self.ddday = ddday

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
    
    def update_growth_demand(self):
        """ Update the growth demand of the lesion according to its current growth rate.
        
        Growth demand is a simple product between growth rate (cm2/degree days) and 
        a delta degree days. 
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        ddday = self.ddday
        
        # Update growth rate
        self.update_growth_rate()  
        r = self.current_growth_rate
        # Compute demand
        demand = min(r * ddday, f.Smax - self.surface)
        self.growth_demand = demand        
    
    def update_growth_rate(self):
        """ Update the growth rate of the lesion in cm2/degree days.
        
        Growth rate is low during incubation and faster after this stage.
        If the growth is between the stages, then growth rate is the mean
        between the lower and the faster according to the time spent in each
        stage.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        ddday = self.ddday
        age_dday = self.age_dday
        time_to_chlo = f.degree_days_to_chlorosis
        
        if age_dday < time_to_chlo: 
            r = f.Smin / time_to_chlo
        elif (age_dday - ddday) < time_to_chlo:
            r1 = f.Smin / time_to_chlo
            r2 = f.growth_rate
            diff1 = time_to_chlo - (age_dday - ddday)
            diff2 = age_dday - time_to_chlo
            r = (diff1*r1 + diff2*r2)/ddday
        else:
            r = f.growth_rate

        self.current_growth_rate = r
        
    def exchange_surfaces(self):
        """ Compute the exchanges of surfaces between rings.
        
        At each time step a fraction of surface pass from each age class 
        (i.e. ring) to the following.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        ddday = self.ddday
        dring = f.delta_age_ring
        s = self.surface_rings
        
        ds = s * ddday / dring
        # Works only if ddday < 20 dd
        
        # Compute the exchanges of surfaces between rings 
        s[:-1] -= ds[:-1]
        s[1:] += ds[:-1]
    
    def control_growth(self, growth_offer=0.):
        """ Reduce surface of the rings up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Surface available on the leaf for the ring to grow (cm2)
        """
        if self.growth_is_active:
            f = self.fungus
            # Growth offer is added to surface alive
            self.surface_alive += growth_offer 
            # TODO : Improve because duplication of information
            
            # Share growth offer between rings
            if len(self.surface_rings)==0.:
                if self.first_ring.growth_is_active:
                    # Compute growth offer for first ring
                    if self.first_ring.surface + growth_offer < f.Smin:
                        self.first_ring.grow(growth_offer)
                    else:
                        go_first_ring = f.Smin - self.first_ring.surface
                        self.first_ring.grow(go_first_ring)
                        go_other_ring = growth_offer - go_first_ring
                        self.surface_rings = np.append(self.surface_rings, go_other_ring)
                        self.delta_age_ring = f.delta_age_ring
            else:
                self.surface_rings[0]+=growth_offer
                
            if growth_offer < self.growth_demand:
                self.disable_growth()
    
    def update_stock(self):
        """ Update the stock of spores on the lesion.
        
        Fill the stock of spores according to production rate and surface of the rings.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        surface_spo_before = self.surface_spo
        
        # First time first ring sporulates
        if not self.stock_spores:
            self.surface_spo += self.first_ring.surface
        
        # For the other rings
        if len(self.surface_rings)>0:
            self.compute_sporulating_surface()

        # Inputs of the stock
        surface_spo = self.surface_spo
        delta_surface_spo = max(0, surface_spo - surface_spo_before)
        try:
            self.stock_spores += delta_surface_spo * f.production_rate
        except:
            self.stock_spores = delta_surface_spo * f.production_rate

    def compute_sporulating_surface(self):
        """ Compute only the sporulating surface.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        dring = f.delta_age_ring
        s = self.surface_rings
        time_to_spo = f.degree_days_to_necrosis + f.degree_days_to_sporulation
        
        diff = self.delta_age_ring - time_to_spo
        nb_full_rings = floor(diff/dring)
        portion_last_ring = (diff%dring)/dring
        
        if nb_full_rings>0:
            surface_spo = sum(s[-nb_full_rings:]) + portion_last_ring*s[-(nb_full_rings+1)]
        else:
            surface_spo = portion_last_ring*s[-(nb_full_rings+1)]
        self.surface_spo = surface_spo + self.first_ring.surface    
        
    def compute_all_surfaces(self):
        """ Compute all the surfaces in different states of the lesion.
        
        This method sums the surfaces in age classes (i.e. rings or part
        of rings) corresponding to the organization in states of the lesion.
        Example: Sum the rings surfaces between 'time_to_spo' and 'time_to_nec'
        in order to find the necrotic surface of the lesion.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        status = self.status
        dring = f.delta_age_ring
        self_dring = self.delta_age_ring
        s = np.copy(self.surface_rings)
        time_to_nec = f.degree_days_to_necrosis
        time_to_spo = f.degree_days_to_necrosis + f.degree_days_to_sporulation
        
        # Initiation
        surface_inc = 0.
        surface_chlo = 0.
        surface_nec = 0.
        surface_spo = 0.
        
        if status == f.INCUBATING:
            surface_inc = self.first_ring.surface
        
        elif status == f.CHLOROTIC:
            surface_chlo = self.surface_alive
        
        elif status == f.NECROTIC:
            # Compute necrotic surface
            diff = self_dring - time_to_nec
            nb_full_rings = floor(diff/dring)
            portion_last_ring = (diff%dring)/dring
            if nb_full_rings>0:
                surface_nec = sum(s[-nb_full_rings:]) + portion_last_ring*s[-(nb_full_rings+1)]
                s[-nb_full_rings:] = 0.
                s[-(nb_full_rings+1)] *= (1-portion_last_ring)
            else:
                surface_nec = portion_last_ring*s[-1]
                s[-1] *= (1-portion_last_ring)
            surface_nec += self.first_ring.surface
            # Compute chlorotic surface
            surface_chlo = sum(s)
        
        elif status == f.SPORULATING:
            # Compute sporulating surface
            diff = self_dring - time_to_spo
            nb_full_rings = floor(diff/dring)
            portion_last_ring = (diff%dring)/dring
            if nb_full_rings>0:
                surface_spo = sum(s[-nb_full_rings:]) + portion_last_ring*s[-(nb_full_rings+1)]
                s[-nb_full_rings:] = 0.
                s[-(nb_full_rings+1)] *= (1-portion_last_ring)
                self_dring -= dring * nb_full_rings
                s = s[s.nonzero()]
            else:
                surface_spo = portion_last_ring*s[-(nb_full_rings+1)]
                s[-1] *= (1-portion_last_ring)
            surface_spo += self.first_ring.surface
            # Compute necrotic surface
            diff = self_dring - time_to_nec
            nb_full_rings = floor(diff/dring)
            portion_last_ring = (diff%dring)/dring
            if nb_full_rings>0:
                surface_nec = sum(s[-nb_full_rings:]) + portion_last_ring*s[-(nb_full_rings+1)]
                s[-nb_full_rings:] = 0.
                s[-(nb_full_rings+1)] *= (1-portion_last_ring)
            else:
                surface_nec = portion_last_ring*s[-(nb_full_rings+1)]
                s[-1] *= (1-portion_last_ring)
            # Compute chlorotic surface
            surface_chlo = sum(s)

        # Save variables
        self.surface_inc = surface_inc
        self.surface_chlo = surface_chlo
        self.surface_nec = surface_nec
        self.surface_spo = surface_spo

    def emission(self, leaf=None):
        """ Create a list of dispersal units emitted by the lesion.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        
        .. Todo:: Implement a real formalism.
        """
        # Import to break circular reference
        from alinea.alep.septoria import SeptoriaDU
        # TODO : Improve ?
        
        if self.is_stock_available(leaf):
            f = self.fungus
            emissions = []
            stock_available = int(self.stock_spores*2/3.)
            
            # TODO : improve below
            nb_DU_emitted = int(leaf.rain_intensity * self.surface_spo * 1000)
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
            if self.stock_spores < 1000:
                self.stock_spores = 0.

            # Return emissions
            SeptoriaDU.fungus = f
            emissions = [SeptoriaDU(nb_spores = nb_spores_by_DU[i], status='emitted')
                                    for i in range(nb_DU_emitted)]
                                    
            return emissions
        else:
            return []
    
    def is_stock_available(self, leaf):
        """ Check if the stock of DU can be emitted.
        
        DU is free for liberation if :
            - there are DUs in stock_du
            - ring is sporulating
            - relative humidity is above a treshold
            - it is the first hour of rain
            
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
       
        Returns
        -------
        True or False :
            Availability of the stock
        
        """
        f = self.fungus
        
        if (self.stock_spores>0. and self.status == f.SPORULATING and 
            leaf.relative_humidity >= f.rh_min and
            self.first_rain_hour):
            return True
        else:
            return False
    
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
    
    def senescence_response(self):
        """ Compute surface alive and surface dead.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        ddday = self.ddday
        ddday_sen = self.ddday_before_senescence
        time_to_nec = f.degree_days_to_necrosis
        dring = f.delta_age_ring
        s = self.surface_rings
        
        # Stop growth
        self.disable_growth()
        
        # Kill all surfaces under necrosis at senescence occurence
        if self.status <= f.CHLOROTIC:
            self.disable()
            self.first_ring.disable()
            surface_dead = self.surface_alive
        else:
            self.compute_all_surfaces()
            surface_dead = self.surface_chlo
            # Update the list of rings
            diff = self.delta_age_ring - time_to_nec
            nb_full_rings = floor(diff/dring)
            portion_last_ring = (diff%dring)/dring
            s[-(nb_full_rings+1)] *= portion_last_ring
            s = s[-(nb_full_rings+1):]
        
        # Update 'surface_alive' and 'surface_dead'
        self.surface_alive -= surface_dead if self.surface_alive > 0. else 0.
        self.surface_dead = surface_dead
        
        # Complete the age of the lesion up to the end of time step
        self.ddday = ddday - ddday_sen
        self.age_dday += ddday - ddday_sen
        
        # Manage first ring
        self.first_ring.update(lesion=self)
        # Manage the other rings
        if len(self.surface_rings)>0:
            if self.age_dday-f.degree_days_to_chlorosis > self.delta_age_ring:
                self.surface_rings = np.append(self.surface_rings, 0.)
                self.delta_age_ring += f.delta_age_ring
            self.exchange_surfaces()

        # Update stock of spores
        if self.status==f.SPORULATING:
            self.update_stock()
        
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
        return self.surface_dead + self.surface_alive
            
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
        if self.first_ring:
            return self.first_ring.status
            
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
        if self.first_ring:
            self.first_ring.status = value
        
# Rings ###########################################################################
class FirstRing(Ring):
    """ First ring of a lesion of septoria.
    """
    def __init__(self, lesion):
        """
        """
        super(FirstRing, self).__init__()
        f = lesion.fungus
        # Status of the ring
        self.status = f.INCUBATING
        # Surface of the ring
        self.surface = 0.
        # Age of the ring
        self.age_dday = 0.
        # Activity of the ring
        self.is_active = True
        # Growth activity of the ring
        self.growth_is_active = True

    def is_incubating(self, fungus):
        return self.status == fungus.INCUBATING

    def is_chlorotic(self, fungus):
        return self.status == fungus.CHLOROTIC

    def is_necrotic(self, fungus):
        return self.status == fungus.NECROTIC
    
    def update(self, lesion=None):
        """
        """
        f = lesion.fungus
        # Ageing of the ring
        self.age_dday = lesion.age_dday
        if self.is_incubating(fungus=f):
            self.incubating(lesion=lesion)
        elif self.is_chlorotic(fungus=f):
            self.chlorotic(lesion=lesion)
        elif self.is_necrotic(fungus=f):
            self.necrotic(lesion=lesion)
            
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
        time_to_necrosis = f.degree_days_to_chlorosis + f.degree_days_to_necrosis

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
        time_to_spo = (f.degree_days_to_chlorosis + 
                           f.degree_days_to_necrosis + 
                           f.degree_days_to_sporulation)
       
        # Compute status transition to sporulating
        if self.age_dday >= time_to_spo:
            self.status = f.SPORULATING
    
    def grow(self, surface):
        self.surface += surface
        assert self.surface <= 0.03

    def disable_growth(self):
        """ Shut down ring growth activity (turn it to False)
        
        Parameters
        ----------
            None
        """
        self.growth_is_active = False
    
    def disable(self):
        """ Disable all activities of the ring.
        
        Set the activity of the lesion to False and its growth demand to 0.
        
        Parameters
        ----------
            None
        """
        self.is_active = False
        self.growth_demand = 0.

class Parameters(_SeptoriaParameters):
    def __init__(self,**kwds):
        _SeptoriaParameters.__init__(self,SeptoriaExchangingRings,**kwds)

class Disease(_Disease):
    @classmethod
    def parameters(cls, **kwds):
        return Parameters(**kwds)
    
    @classmethod
    def lesion(cls, **kwds):
        SeptoriaExchangingRings.fungus=cls.parameters(**kwds)
        return SeptoriaExchangingRings
 
