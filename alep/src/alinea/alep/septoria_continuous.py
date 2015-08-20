""" Class of lesion specific of wheat septoria, with a continuous model.

"""
# Imports #########################################################################
from alinea.alep.fungal_objects import *
# The following import would provoke a circular reference
# "from alinea.alep.septoria import SeptoriaDU"
# --> Moved in the method 'ContinuousSeptoria.emission()'
from alinea.alep.septoria import Disease as _Disease, SeptoriaParameters as _SeptoriaParameters
from random import randint, seed
import numpy as np
seed(1)

# Lesion ##########################################################################
class ContinuousSeptoria(Lesion):
    """ Septoria Lesion implemented as a continuous model. """

    def __init__(self, nb_spores=None, position=None):
        """ Initialize the lesion of septoria. 
        
        Parameters
        ----------
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        position: non defined
            Position of the dispersal unit on the phyto-element
        """
        super(ContinuousSeptoria, self).__init__(nb_spores=nb_spores, position=position)
        # Status of the lesion
        self.status = self.fungus.INCUBATING
        # Surface alive of the lesion
        self.surface_alive = 0.
        # Surfaces in each state
        self.surface_inc = 0.
        self.surface_chlo = 0.
        self.surface_nec = 0.
        self.surface_spo = 0.
        # Surface sporulating the time step before
        self.surface_spo_before = 0.
        # Age of the center of the lesion (degree days)
        self.age_tt = 0.
        # Delta age for each time step
        self.ddday = 0.
        # Position of senescence the time step before (Useful in case of senescence
        # to compute the time left for growth before senescence occur)
        self.old_position_senescence = None
        # Degree days before senescence if 'self.is_senescent'
        self.ddday_before_senescence = None
        # Surface killed by senescence
        self.surface_dead = 0.
        # Growth activity of the lesion
        self.growth_is_active = True
        # Growth rate as a function of lesion status (cm2/degree days)
        self.current_growth_rate = 0.        
        # Growth demand of the lesion (cm2)
        self.growth_demand = 0.
        # Is first hour of rain
        self.first_rain_hour = False
        # Stock of spores (number of spores)
        self.stock_spores = 0.
        # List of DUs emitted
        self.emissions = []

    def update(self, dt=1., leaf=None):
        """ Update the status, age and growth demand of the lesion.

        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. area, green area, healthy area,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        f = self.fungus

        # Disable the lesion if all is sporulating and stock is empty
        if (self.age_tt > 0. and
            self.surface_spo != 0. and
            self.surface_spo == self.surface_alive and
            self.stock_spores == 0.):
            self.disable()
            return
        
        # Compute delta degree days in dt
        # TODO : modify if list of temp coming from weather data
        self.compute_delta_ddays(dt, leaf)
        ddday = self.ddday

        # If senescence, compute length of growth period before senescence during time step
        if self.is_senescent:
            self.compute_time_before_senescence(ddday=ddday, leaf=leaf)
            ddday = self.ddday_before_senescence
            
        # Update the age of the lesion
        self.age_tt += ddday
        
        # Update the status of the lesion
        self.update_status()
        
        # Update the growth variables of the lesion
        if self.growth_is_active:
            # Update growth demand
            self.update_growth_demand()
        else:
            self.current_growth_rate = 0.
            self.growth_demand = 0.        

        # Update the perception of rain by the lesion
        if self.status == f.SPORULATING:
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
            A leaf sector with properties (e.g. area, green area, healthy area,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        f = self.fungus
        # Calculation
        if dt != 0.:
            ddday = max(0,(leaf.temp - f.basis_for_dday))*(dt/24.)
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
            A leaf sector with properties (e.g. area, green area, healthy area,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        old_pos = self.old_position_senescence
        new_pos = leaf.position_senescence
        speed = abs(new_pos-old_pos)/ddday if ddday > 0. else 0.
        
        self.ddday_before_senescence = abs(self.position[0]-old_pos)/speed if speed >0. else 0.
    
    def update_status(self):
        """ Find the status of a lesion of septoria according to its age in degree days.
        
        Before 220 DD, a lesion of septoria is INCUBATING.
        ______ 330 DD ________________________ CHLOROTIC.
        ______ 350 DD ________________________ NECROTIC.
        After, it is SPORULATING if not EMPTY or DEAD.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        status: int
            Status of the lesion
        """
        f = self.fungus
        age_tt = self.age_tt
        status = [f.SPORULATING, f.INCUBATING, f.CHLOROTIC, f.NECROTIC]
        times = [0,f.degree_days_to_chlorosis, f.degree_days_to_necrosis, f.degree_days_to_sporulation]
        times = np.cumsum(times)
        
        self.status = status[np.argmin(times<=age_tt)] 

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
        age_tt = self.age_tt
        time_to_chlo = f.degree_days_to_chlorosis
        
        if age_tt < time_to_chlo: 
            r = f.Smin / time_to_chlo
        elif (age_tt - ddday) < time_to_chlo:
            r1 = f.Smin / time_to_chlo
            r2 = f.growth_rate
            diff1 = time_to_chlo - (age_tt - ddday)
            diff2 = age_tt - time_to_chlo
            r = (diff1*r1 + diff2*r2)/ddday
        else:
            r = f.growth_rate

        self.current_growth_rate = r
    
    def compute_all_surfaces(self):
        """ Compute all the surfaces in different states of the lesion.
        
        Before chlorosis, a lesion of septoria grows slowly at a rate 'r1'.
        After, it grows at a higher rate 'r2'. 
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        surface_alive = self.surface_alive
        status = self.status
        age_tt = self.age_tt
        r = f.growth_rate
        Smin = f.Smin
        
        time_to_spo = f.degree_days_to_chlorosis + f.degree_days_to_necrosis + f.degree_days_to_sporulation
        time_to_nec = f.degree_days_to_chlorosis + f.degree_days_to_necrosis
        time_to_chlo = f.degree_days_to_chlorosis
        
        # Initiation
        surface_inc = 0.
        surface_chlo = 0.
        surface_nec = 0.
        surface_spo = 0.
        
        # Compute surfaces
        if status == f.INCUBATING:
            surface_inc = surface_alive
            
        elif status == f.CHLOROTIC:
            surface_chlo = surface_alive
            
        elif status == f.NECROTIC:
            delta_age_nec = age_tt - time_to_nec
            # Potential surface in necrosis if no interruption
            pot_surface_nec = Smin + r*delta_age_nec
            if surface_alive < pot_surface_nec:
                # All the lesion is necrotic
                surface_nec = surface_alive
            else:
                # There are necrotic and chlorotic surfaces
                surface_nec = pot_surface_nec
                surface_chlo = surface_alive - surface_nec

        elif status == f.SPORULATING:
            delta_age_spo = age_tt - time_to_spo
            # Potential surface in necrosis if no interruption
            pot_surface_spo = Smin + r*delta_age_spo
            if surface_alive < pot_surface_spo:
                # All the lesion is necrotic
                surface_spo = surface_alive
            else:
                # There are at least sporulating and necrotic surfaces
                surface_spo = pot_surface_spo
                # Potential surface in necrosis if no interruption
                pot_surface_nec = r*f.degree_days_to_necrosis
                if (surface_alive - surface_spo) < pot_surface_nec:
                    # All the rest is necrotic
                    surface_nec = surface_alive - surface_spo
                else:
                    # There are at sporulating, necrotic and chlorotic surfaces
                    surface_nec = pot_surface_nec
                    surface_chlo = surface_alive - surface_spo - surface_nec
        
        # Save variables
        self.surface_inc = surface_inc
        self.surface_chlo = surface_chlo
        self.surface_nec = surface_nec
        self.surface_spo = surface_spo
    
    def compute_sporulating_surface(self):
        """ Compute only the sporulating surface.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        surface_alive = self.surface_alive
        status = self.status
        age_tt = self.age_tt
        r = f.growth_rate
        Smin = f.Smin
        time_to_spo = f.degree_days_to_chlorosis + f.degree_days_to_necrosis + f.degree_days_to_sporulation
        
        # Initiation
        surface_spo = 0.
        
        if self.status == f.SPORULATING:
            delta_age_spo = age_tt - time_to_spo
            # Potential surface in necrosis if no interruption
            pot_surface_spo = Smin + r*delta_age_spo
            if surface_alive < pot_surface_spo:
                # All the lesion is necrotic
                surface_spo = surface_alive
            else:
                # There are at least sporulating and necrotic surfaces
                surface_spo = pot_surface_spo
                
        # Save sporulating surface
        self.surface_spo = surface_spo
    
    def control_growth(self, growth_offer = 0.):
        """ update surface of the lesion according to its growth demand
            and available surface on the leaf.
        
        Parameters
        ----------
        growth_offer: float
            Minimum between 'growth_demand' and the surface available on
            the leaf for the lesion to grow (cm2)
        """
        f = self.fungus
        growth_demand = self.growth_demand

        if self.growth_is_active:
            # Growth offer is added to surface alive
            self.surface_alive += growth_offer
            
            # Check if any interruption of growth:
            if growth_offer < growth_demand or self.surface == f.Smax:
                self.disable_growth()

        # Update the production of spores of the lesion
        if self.status == f.SPORULATING:
            self.update_stock()
    
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
        
        # Complete the age of the lesion up to the end of time step
        # Note : Age was stopped for update at the time of senescence occurence
        self.age_tt += ddday - ddday_sen
        
        # Stop growth
        self.disable_growth()
        
        # Kill all surfaces under necrosis at senescence occurence
        if self.status <= f.CHLOROTIC:
            self.disable()
            surface_dead = self.surface_alive
        else:
            self.compute_all_surfaces()
            surface_dead = self.surface_chlo
            # temp
            # assert self.surface_alive - surface_dead == self.surface_nec + self.surface_spo # temp
        
        # Update 'surface_alive' and 'surface_dead'
        self.surface_alive -= surface_dead if self.surface_alive > 0. else 0.
        self.surface_dead += surface_dead
        
    def update_stock(self):
        """ Update the stock of spores produced by the lesion.
        
        Each time a new surface is sporulating, it produces a given
        amount of spores by unit of surface : parameter 'production_rate'.
        
        Parameters
        ----------
        None        
        """
        f = self.fungus
        
        # Compute sporulating surface
        self.compute_sporulating_surface()
        
        surface_spo_before = self.surface_spo_before
        surface_spo = self.surface_spo
        
        # Inputs of the stock
        delta_surface_spo = max(0, surface_spo - surface_spo_before)
        self.stock_spores += delta_surface_spo * f.production_rate
        
        # Save 'surface_spo_before' for next time step
        self.surface_spo_before = surface_spo

    def is_stock_available(self, leaf):
        """ Check if the stock of DU can be emitted.
        
        DU is free for liberation if :
            - there are DUs in stock_du
            - lesion is sporulating
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
    
    def emission(self, leaf = None):
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
        
        f = self.fungus        
        if self.is_stock_available(leaf):
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
            
    @property
    def surface(self):
        """ Compute the surface of the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        surface: float
            Surface of the lesion
        """
        return self.surface_alive + self.surface_dead
    
    @property
    def necrotic_area(self):
        """ Compute the necrotic area of the lesion.
        
        Necrotic area is composed by surfaces in state:
            - NECROTIC
            - SPORULATING
        
        Parameters
        ----------
            None
            
        Returns
        -------
        status: int
            Status of the lesion
        """
        self.compute_all_surfaces()
        return self.surface_nec + self.surface_spo
    
class Parameters(_SeptoriaParameters):
    def __init__(self,**kwds):
        _SeptoriaParameters.__init__(self,ContinuousSeptoria,**kwds)

class Disease(_Disease):
    @classmethod
    def parameters(cls, **kwds):
        return Parameters(**kwds)
    
    @classmethod
    def lesion(cls, **kwds):
        ContinuousSeptoria.fungus=cls.parameters(**kwds)
        return ContinuousSeptoria