""" Classes of dispersal unit, lesion and ring specific of wheat septoria.

"""
# Imports #########################################################################
from alinea.alep.cycle2 import *
from random import random, randint
from math import floor, ceil
import numpy as np

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

        """
        super(Septoria, self).__init__(nb_spores=nb_spores, position=position)
        self.rings = []
        ring = SeptoriaRing(lesion = self, status = self.fungus.INCUBATING, dt=1.)
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

        """
        assert self.active
        f = self.fungus
        
        # remove non active rings
        self.surface_dead += sum(ring.surface for ring in self.rings if not ring.active)
        self.rings = [ring for ring in self.rings if ring.active] 

        # Disable lesion when no ring left
        if not self.rings:
            self.disable()
            return

        # Delta degree days
        ddday = max(0, (leaf.temp - f.basis_for_dday)/(24./dt))
        
        # Compute current Dt: width of a ring in degree days
        if len(self.rings) == 1 and self.surface_dead == 0.:
            Dt = f.degree_days_to_chlorosis
        else:
            Dt = f.Dt
        
        # Update every active ring of the lesion except the last because it can generate a new ring
        for ring in self.rings[:-1]:
            ring.update(dt=dt, ddday=ddday, leaf=leaf, lesion=self, **kwds)
        
        # Manage the ring in formation and create a new ring when needed
        if self.can_form_new_ring(Dt, ddday):
            # A new ring will be added
            # 1. Compute age_dday of the two last rings and their ddday for growth demand
            # ( + Add leaf in new_ring properties)
            old_ring = self.rings[-1]
            old_ring.ddday = Dt - old_ring.age_dday
            old_ring.age_dday += ddday
            old_ring.need_last_growth = True
            
            new_ring = SeptoriaRing(lesion=self, status=f.CHLOROTIC, dt=dt)
            new_ring.leaf = leaf
            new_ring.age_dday = ddday - old_ring.ddday
            new_ring.ddday = new_ring.age_dday
            
            # 2. Compute their growth_demand with their new ddday
            old_ring.stage(dt=dt, lesion=self)
            new_ring.stage(dt=dt, lesion=self)
            assert new_ring.is_in_formation(f)
            
            # 3. Update the list of rings of the lesion
            self.rings[-1] = old_ring
            self.rings.append(new_ring)
        else:
            self.rings[-1].update(dt=dt, ddday=ddday, leaf=leaf, lesion=self, **kwds)          
        
    def can_form_new_ring(self, Dt, ddday):
        """ Check if the lesion can form a new ring.
        
        A new ring is created when the physiologic age is reached.
        
        Parameters
        ----------
        Dt : int
            Time step in degree days to create a new ring
        ddday: float
            Number of degree days in time step
        
        """
        
        if self.surface < self.fungus.Smax:
            if self.rings[-1].age_dday + ddday < Dt :
                return False
            else:
                return True
        else:
            return False
        
        # Keep it in mind --> To be integrated later

        # if( leaf.temp < self.fungus.basis_for_dday) ):
            # return False
        # else:
            # return True
        
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
           if ring.is_stock_available(leaf, lesion = self):
                ring_emissions = ring.emission(leaf, lesion=self)
                if ring_emissions:
                    emissions += ring_emissions
        
        return emissions
    
    def control_growth(self, growth_offer = 0.):
        """ Reduce surface of the last ring up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Surface available on the leaf for the ring to grow (cm2)
        
        """
        fungus = self.fungus

        if self.rings:
            forming_rings = [r for r in self.rings if r.is_in_formation(fungus)]
            if len(forming_rings) == 2:
                # A new ring is emerging at this time step: The last one on the list still has
                # a growth demand, the new one has a growth demand two but will be served after
                forming_rings[0].control_growth(lesion=self, 
                                                growth_offer=max(growth_offer, forming_rings[0].growth_demand))
                forming_rings[0].need_last_growth = False                                
                forming_rings[1].control_growth(lesion=self, 
                                                growth_offer=max(0., growth_offer - forming_rings[0].growth_demand))           
            else:
                self.rings[-1].control_growth(lesion=self, growth_offer=growth_offer)
    
    def is_dead(self):
        """ Update the status of all the rings to 'DEAD' if the lesion is dead.
        
        Parameters
        ----------
            None
        
        """
        return all(ring.is_dead() for ring in self.rings)
        
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
            return sum(r.growth_demand for r in forming_rings)
        else:
            return 0.
  
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
    def stock(self):
        """ Compute the stock of spores on the lesion.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        surface: float
            Surface of the whole lesion (cm2)
        """
        stock = sum(ring.stock for ring in self.rings)
        return stock
    
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

class ContinuousSeptoria(Lesion):
    """ Septoria Lesion implemented as a continuous model. """

    def __init__(self, nb_spores=None, position=None):
        """ Initialize the lesion of septoria. 
        
        Parameters
        ----------
        position: non defined
            Position of the dispersal unit on the phyto-element
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        
        """
        super(ContinuousSeptoria, self).__init__(nb_spores=nb_spores, position=position)
        # Status of the lesion
        self.status = self.fungus.INCUBATING
        # Surface of the lesion
        self.surface = 0.
        # Surfaces in each state
        self.surface_inc = 0.
        self.surface_chlo = 0.
        self.surface_nec = 0.
        self.surface_spo = 0.
        # Surface sporulating the time step before
        self.surface_spo_before = 0.
        # Age of the center of the lesion (degree days)
        self.age_dday = 0.
        # Delta age for each time step
        self.ddday = 0.
        # Growth rate as a function of lesion status (cm2/degree days)
        self.current_growth_rate = 0.        
        # Growth demand of the lesion (cm2)
        self.growth_demand = 0.
        # List of continuous sequences of time (degree days) : [growth stop growth stop ...]
        self.sequences = [0]
        # Stock of spores (number of spores)
        self.stock = 0.
        # Is first hour of rain
        self.first_rain_hour = False
        # List of DUs emitted
        self.emissions = []

    def update(self, dt=1., leaf=None, **kwds):
        """ Update the status, surface and stock of spores of the lesion.

        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)

        """       
        f = self.fungus
        
        # Delta degree days
        ddday = max(0,(leaf.temp - f.basis_for_dday)/(24./dt))
        self.ddday = ddday
        # Update growth rate
        self.update_growth_rate(ddday)        
        
        self.age_dday += ddday
        
        # Update growth demand
        self.update_growth_demand(ddday)
        # Update the status of the lesion
        self.update_status()
        # Update the surfaces of the lesion in different states
        self.update_surfaces()  
        # Update the production of the lesion
        if self.status == f.SPORULATING:
            if leaf.rain_intensity > 0. and leaf.relative_humidity >= f.rh_min:
                self.first_rain_hour = True if not self.first_rain_hour else False
            else:
                self.first_rain_hour = False
            
            self.update_stock()   

            # Disable the lesion if all is sporulating and stock is empty
            if (self.age_dday > 0. and
                self.surface_inc == 0. and
                self.surface_chlo == 0. and
                self.surface_nec == 0. and
                self.stock == 0.):
                self.disable()
        
    def update_growth_rate(self, ddday):
        """ Update the growth rate of the lesion in cm2/degree days.
        
        Growth rate is low during incubation and faster after this stage.
        If the transition is between the stages, then growth rate is a 
        mean of the lower and the faster.
        
        Parameters
        ----------
        ddday: float
            Delta degree days.
        """
        f = self.fungus
        age_dday = self.age_dday
        time_to_chlo = f.degree_days_to_chlorosis
       
        if age_dday < time_to_chlo: 
            r1 = f.Smin / time_to_chlo
            if (age_dday + ddday) < time_to_chlo:
                r = r1
            else:
                r2 = f.growth_rate
                diff1 = time_to_chlo - age_dday
                diff2 = age_dday + ddday - time_to_chlo
                r = (diff1*r1 + diff2*r2)/ddday
        else:
            r = f.growth_rate
        
        # If no surface is incubating or chlorotic, growth stops
        # TODO : Biological hypothesis to confirm
        if self.surface_inc == self.surface_chlo == 0.:
            r == 0.
        
        self.current_growth_rate = r
    
    def update_growth_demand(self, ddday=0.):
        """ Update the growth demand of the lesion according to its current growth rate.
        
        Growth demand is a simple product between growth rate (cm2/degree days) and 
        a delta degree days. Lesion growth is limited by the parameter Smax.
        
        Parameters
        ----------
        ddday: float
            Delta degree days.
        
        """
        Smax = self.fungus.Smax
        demand = self.current_growth_rate * ddday
        
        potential_surface = self.surface + demand
        if potential_surface > Smax :
            demand = Smax - self.surface
        
        self.growth_demand = demand

    def update_status(self):
        """ Update the status of the lesion.
        
        Before 220 DD, a lesion of septoria is INCUBATING.
        ______ 330 DD ________________________ CHLOROTIC.
        ______ 350 DD ________________________ NECROTIC.
        After, it is SPORULATING if not EMPTY or DEAD.
        
        Parameters
        ----------
            None
        
        TODO: add another time for competition
        Manage the time for each state when competition
        """
        f = self.fungus

        status = [f.SPORULATING, f.INCUBATING, f.CHLOROTIC, f.NECROTIC]
        times = [0,f.degree_days_to_chlorosis, f.degree_days_to_necrosis, f.degree_days_to_sporulation]
        times = np.cumsum(times)
        
        status = status[np.argmin(times<=self.age_dday)] 
        
        self.status = status
    
    def update_surfaces(self):
        """ Update the surface of the lesion.
        
        Before chlorosis, a lesion of septoria grows slowly at a rate 'r1'.
        After, it grows at a higher rate 'r2'. 
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        surface = self.surface
        status = self.status
        age_dday = self.age_dday
        sequences = self.sequences
        growths = sequences[0::2]
        stops = sequences[1::2]
        r = self.current_growth_rate
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
            surface_inc = surface
            
        elif status == f.CHLOROTIC:
            surface_chlo = surface
            
        elif status == f.NECROTIC:
            if sum(growths) < time_to_chlo:
                # The whole surface is the same age and is necrotic (plateau of the curve)
                surface_nec = surface
            elif len(sequences) == 1.:
                # No interruption --> Easy calculation :
                delta_age_nec = sequences[0] - time_to_nec
                surface_nec = Smin + r * delta_age_nec
                # The remaining surface is chlorotic
                surface_chlo = surface - surface_nec
                
                # Temp:
                surface_chlo = r * delta_age_chlo
                assert delta_age_chlo == time_to_nec
            else:
                # At least Smin is necrotic. Remove time to Smin from list of sequences
                ind = np.where(np.cumsum(growths) > time_to_chlo)[0][0]
                sequences[ind] = sum(growths[0:ind+1]) - time_to_chlo
                sequences = sequences[ind*2:]
                                
                # Calculation of the delta time where necrosis have grown
                # delta_age_nec excludes times of interruption
                diff = age_dday - time_to_nec
                ind = np.where(np.cumsum(sequences) > diff)[0][0]
                delta_age_nec = sum(sequences[0:ind+1:2])
                
                # Deduction of necrotic surface
                surface_nec = Smin + r * delta_age_nec
                # The remaining surface is chlorotic
                surface_chlo = surface - surface_nec
                
        elif status == f.SPORULATING:
            if sum(growths) < time_to_chlo:
                # The whole surface is the same age and is sporulating (plateau of the curve)
                surface_spo = self.surface
            else:
                if len(sequences) == 1.:
                    # No interruption --> Easy calculation :
                    age_spo = sequences[0] - time_to_spo
                    surface_spo = Smin + r * age_spo
                    
                    age_nec = time_to_spo - time_to_nec
                    surface_nec = r * age_nec
                    
                    age_chlo = time_to_nec
                    surface_chlo = r * age_chlo
                else:
                    # At least Smin is sporulating. Remove time to Smin from list of sequences
                    ind = np.where(np.cumsum(growths) > time_to_chlo)[0][0]
                    
                    # Reduction of the list of sequences
                    sequences[ind] = sum(growths[0:ind+1]) - time_to_chlo
                    sequences = sequences[ind*2:]
                                    
                    # Calculation of the delta time where sporulation have grown
                    # delta_age_spo excludes times of interruption
                    diff = age_dday - time_to_spo
                    ind = np.where(np.cumsum(sequences) > diff)[0][0]
                    delta_age_spo = sum(sequences[0:ind+1:2])
                    
                    # Deduction of sporulating surface
                    surface_spo = Smin + r * delta_age_spo
                    
                    # Remove time to delta_age_spo from list of sequences 
                    sequences = sequences[ind:]
                    if ind%2 == 0.:
                        is_growing = True
                    else:
                        is_growing = False 
                    
                    # Calculation of the delta time where necrosis have grown
                    # delta_age_nec excludes times of interruption
                    diff = f.degree_days_to_necrosis
                    ind = np.where(np.cumsum(sequences) > diff)[0][0]
                    if is_growing:
                        delta_age_spo = sum(sequences[0:ind+1:2])
                    else:
                        delta_age_spo = sum(sequences[1:ind+1:2])
                        
                    # The remaining surface is chlorotic
                    surface_chlo = surface - surface_spo - surface_nec
        
        # Save variables
        self.surface_inc = surface_inc
        self.surface_chlo = surface_chlo
        self.surface_nec = surface_nec
        self.surface_spo = surface_spo
        
    def control_growth(self, growth_offer = 0.):
        """ Reduce surface of the lesion up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Minimum between the surface available on the leaf for the lesion
            to grow (cm2) and 'growth_demand'

        """
        f = self.fungus
        age_dday = self.age_dday
        ddday = self.ddday
        sequences = self.sequences
        growth_demand = self.growth_demand
        time_to_chlo = f.degree_days_to_chlorosis
        r1 = f.Smin / time_to_chlo
        r2 = f.growth_rate
        
        # Update total surface of the lesion
        self.surface += growth_offer
        
        # Age before last growth
        age_before = age_dday - ddday
        
        if growth_offer == growth_demand:
            # No limitation of growth
            if self.was_growing_before():
                sequences[-1] += ddday
            else:
                sequences.append(ddday)
                
        elif growth_offer > 0.:
            # The lesion can grow but is limited
            if growth_offer < r1*(time_to_chlo - age_before):
                # Competition occured before reaching chlorosis
                # Growth rate has not changed yet and equals r1
                time_of_growth = growth_offer/r1
            else:
                # Competition occured after reaching chlorosis
                if age_before > time_to_chlo:
                    # Growth rate has not changed in the last time step and equals r2
                    time_of_growth = growth_offer/r2
                else:
                    # Growth rate has changed in the last time step
                    time_with_r1 = time_to_chlo - age_before
                    time_with_r2 = (growth_offer - r1*time_with_r1)/r2
                    time_of_growth = time_with_r1 + time_with_r2
                    assert time_of_growth < ddday
            
            time_of_stops = ddday - time_of_growth
            
            # Manage the list of sequences
            if self.was_growing_before():
                sequences[-1] += time_of_growth
                sequences.append(time_of_stops)
            else:
                sequences.append(time_of_growth)
                sequences.append(time_of_stops)    
        
        else:
            # The lesion is stopped
            if self.was_growing_before():
                sequences.append(ddday)
            else:
                sequences[-1] += ddday
        
    def was_growing_before(self):
        """ Check if the lesion was in a phase of growth or stop before competition occured.
        
        The lesion saves a list of sequences of alternative growth periods and stop periods.
        If the length of the list is pair, then the lesion was is a phase of growth. Otherwise
        it was stopped.
        
        Parameters
        ----------
            None
            
        Returns
        -------
        True or False
            True if the lesion was growing before, False otherwise.
        
        """
        if len(self.sequences)%2 == 0.:
            return True
        else:
            return False
                
    def update_stock(self):
        """ Update the stock of spores produced by the lesion.
        
        Each time a new surface is sporulating, it produces a given
        amount of spores by unit of surface : parameter 'production_rate'.
        
        Parameters
        ----------
        None        
        """
        f = self.fungus
        surface_spo_before = self.surface_spo_before
        surface_spo = self.surface_spo
        
        # Inputs of the stock
        delta_surface_spo = max(0, surface_spo - surface_spo_before)
        self.stock += delta_surface_spo * f.production_rate
        
        # Update 'surface_spo_before' for next time step
        self.surface_spo_before = surface_spo

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
        
        if (self.stock and self.status == f.SPORULATING and 
            leaf.relative_humidity >= f.rh_min and
            self.first_rain_hour):
            return True
        else:
            return False
    
    def emission(self, leaf = None):
        """ Create a list of dispersal units emitted by the ring.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        
        .. Todo:: Implement a real formalism.
        """
        if self.is_stock_available(leaf):
            f = self.fungus
            emissions = []
            stock_available = self.stock/3.
            
            # TODO : improve below
            nb_DU_emitted = int(leaf.rain_intensity * self.surface * 100)
            nb_spores_by_DU = []
            for DU in range(nb_DU_emitted):
                if stock_available > 0.:
                    nb_spores = min(randint(5,100), stock_available)
                    nb_spores_by_DU.append(nb_spores)
                    stock_available -= nb_spores
            
            # Get rid of DUs without spores
            nb_DU_emitted = len(nb_spores_by_DU)
            
            # Update stock
            self.stock -= stock_available
            
            if self.stock < 100:
                self.stock = 0.

            # Return emissions
            emissions = [SeptoriaDU(nb_spores = nb_spores_by_DU[i], status='emitted')
                                    for i in range(nb_DU_emitted)]
            return emissions
        else:
            return []
    
    # @property
    # def status(self):
        # """ Compute the status of the lesion.
        
        # Parameters
        # ----------
            # None
            
        # Returns
        # -------
        # status: int
            # Status of the lesion
        # """
        # return self.status

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
        self.need_last_growth = False
        
        # TODO : integrate the notion of epsilon ??
        # if not lesion.rings:
            # self.surface = lesion.fungus.epsilon # TODO : see later with growth_demand
        # else:
            # self.surface = 0.
        
        self.surface = 0.
        
        self.age = 0.
        self.age_dday = 0.
        self.ddday = 0.
        self.cumul_rain_event = 0.
        self.first_rain_hour = False
        self.rain_before = False
        self.stock = -1

    def is_in_formation(self, fungus):
        """ Can keep growing!!! """
        ok = self.is_incubating(fungus)
        if not ok:
            ok = self.is_chlorotic(fungus) and (self.age_dday <= fungus.Dt or self.need_last_growth)
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

    def update(self, dt, ddday, leaf, lesion=None):
        """ Update the status of the ring.
        
        * Cumulate the age of the ring.
        * Assign leaf data to the ring in order to access it in the methods.
        * Call the property 'stage' depending on the current status of the ring.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        ddday: float
            Delta degree days during 'dt'
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
            
        Returns
        -------
            None
          
        """
        f = lesion.fungus
        self.leaf = leaf
        
        # MANAGE TIME SCALE
        self.age += dt
        self.age_dday += ddday
        self.ddday = ddday
        
        # Check if it is the first hour of rain
        if leaf.rain_intensity > 0. and leaf.relative_humidity >= f.rh_min:
            self.first_rain_hour = True if not self.first_rain_hour else False
        else:
            self.first_rain_hour = False
        
        # Compute status of the ring
        self.stage(dt=dt, lesion=lesion)
    
    def can_form_new_ring(self, Dt, lesion):
        """ Check if the lesion can form a new ring.
        
        A new ring is created when the physiologic age is reached.
        
        Parameters
        ----------
        Dt : int
            Time step in degree days to create a new ring
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        """
        
        if lesion.surface < lesion.fungus.Smax:
            if lesion.rings[-1].age_dday < Dt :
                return False
            else:
                return True
        else:
            return False
    
    def in_formation(self, lesion=None, **kwds):
        """ Cumulate surface while Dt in the state is not reached.        
        
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
        ddday = self.ddday
       
        if self.is_incubating(fungus=f):
            # First ring in incubation
            self.growth_demand = f.Smin * ddday / f.degree_days_to_chlorosis
        elif (self.is_chlorotic(fungus=f) and 
              self.age_dday >= f.degree_days_to_chlorosis and
              (self.age_dday-ddday) < f.degree_days_to_chlorosis):
            # Second ring when it emerges
            if ddday != 0. :
                diff = f.degree_days_to_chlorosis - (self.age_dday-ddday)
                self.growth_demand =( (diff*f.Smin/f.degree_days_to_chlorosis +
                                       (ddday-diff)*f.growth_rate) / ddday)
            else:
                self.growth_demand = 0.
        else:
            # Second ring after emergence and all other rings
            size_before_Smax = f.Smax - lesion.surface
            self.growth_demand = min(size_before_Smax, f.growth_rate * ddday)
        
        # print('growth demand %f' % self.growth_demand)
        
    def control_growth(self, lesion=None, growth_offer = 0.):
        """ Reduce surface of the last ring up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Minimum between the surface available on the leaf for the ring
            to grow (cm2) and 'growth_demand'

        """
        self.surface += growth_offer
        
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
        fungus = lesion.fungus
        assert self.is_chlorotic(fungus)
        
        self.growth_demand = 0.
        
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

        .. Todo:: Enhance the code to parametrize the code with a function.
        """
        fungus = lesion.fungus
        leaf = self.leaf

        assert self.is_sporulating(fungus)
        
        # (float) : rain intensity on the leaf sector during the time step (in mm/h).
        rain_intensity = leaf.rain_intensity 
        # (float) : relative humidity on the leaf sector during the time step (in %).
        relative_humidity = leaf.relative_humidity 
        
        # Count dispersal events
        if self.first_rain_hour:
            self.cumul_rain_event += 1
        
        if self.cumul_rain_event >= fungus.rain_events_to_empty:
            self.empty(lesion) # The ring is empty.
        
        # Fill the stock of spores according to production rate and surface of the ring
        if self.stock==-1.:
            production = self.surface * fungus.production_rate
            # nb_spores_produced = int(floor(production))
            # if proba(1 - (ceil(production) - production)):
                # nb_spores_produced += 1
            
            nb_spores_produced = int(round(production))
            
            self.stock += 1 + nb_spores_produced
                    
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
        stock_available = self.stock
        
        # TODO : improve below
        nb_DU_emitted = int(leaf.rain_intensity * self.surface * 100)
        nb_spores_by_DU = []
        for DU in range(nb_DU_emitted):
            if stock_available > 0.:
                nb_spores = min(randint(5,100), stock_available)
                nb_spores_by_DU.append(nb_spores)
                stock_available -= nb_spores
                
        # Get rid of DUs without spores
        nb_DU_emitted = len(nb_spores_by_DU)
        
        # Update stock
        self.stock -= stock_available
        if self.stock < 1.:
            self.status == fungus.EMPTY

        # Return emissions
        emissions = [SeptoriaDU(nb_spores = nb_spores_by_DU[i], status='emitted')
                                for i in range(nb_DU_emitted)]
        
        return emissions
    
    def is_stock_available(self, leaf, lesion=None):
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
        lesion : Lesion instantiation
            The lesion carrying the ring, with properties 
            (e.g. fungus parameters, surface, status, age, rings, etc.)
        
        Returns
        -------
        True or False :
            Availability of the stock
        
        """
        fungus = lesion.fungus
        
        # print(self.first_rain_hour)
        
        if (self.stock and self.status == fungus.SPORULATING and 
            leaf.relative_humidity >= fungus.rh_min and
            self.first_rain_hour):
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
        
    def stage(self, dt, lesion=None):
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
        f = lesion.fungus
        
        # Particular case of first ring:
        # Need last growth when status is modified. Is still in formation 
        if (self.is_incubating(f) and
            f.degree_days_to_chlorosis <= self.age_dday < f.degree_days_to_chlorosis + f.degree_days_to_necrosis):
            self.status = f.CHLOROTIC
            self.need_last_growth = True
    
        if self.is_in_formation(fungus=f):
            return self.in_formation(lesion=lesion)
        elif self.is_chlorotic(fungus=f):
            return self.chlorotic(lesion=lesion)
        elif self.is_necrotic(fungus=f):
            return self.necrotic(lesion=lesion)
        elif self.is_sporulating(fungus=f):
            return self.sporulating(lesion=lesion)
        elif self.is_empty(fungus=f):
            return self.empty(lesion=lesion)
        else:
            return self.dead(lesion=lesion)

# Fungus parameters (e.g. .ini): config of the fungus #############################
class SeptoriaParameters(Parameters):
    def __init__(self,
                 INCUBATING = 0,
                 CHLOROTIC = 1,
                 NECROTIC = 2,
                 SPORULATING = 3,
                 EMPTY = 4,
                 DEAD = 5,
                 Dt = 20,
                 basis_for_dday = -2.,
                 temp_min = 10.,
                 temp_max = 30.,
                 wd_min = 10.,
                 loss_rate = 1./120,
                 degree_days_to_chlorosis = 220.,
                 degree_days_to_necrosis = 110.,
                 degree_days_to_sporulation = 20.,
                 epsilon = 0.001,
                 Smin = 0.03,
                 Smax = 0.3,
                 growth_rate = 0.0006,
                 rh_min = 85.,
                 rain_events_to_empty = 3,
                 production_rate = 100000,
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
        self.INCUBATING = INCUBATING
        self.CHLOROTIC = CHLOROTIC
        self.NECROTIC = NECROTIC
        self.SPORULATING = SPORULATING
        self.EMPTY = EMPTY
        self.DEAD = DEAD
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

    # def __call__(self, nb_spores = None, position = None):
        # if Septoria.fungus is None:
            # Septoria.fungus = self
        # if SeptoriaDU.fungus is None:
            # SeptoriaDU.fungus = self
        # return Septoria(nb_spores=nb_spores, position=position)
        
    def __call__(self, nb_spores = None, position = None):
        if ContinuousSeptoria.fungus is None:
            ContinuousSeptoria.fungus = self
        if SeptoriaDU.fungus is None:
            SeptoriaDU.fungus = self
        return ContinuousSeptoria(nb_spores=nb_spores, position=position)

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
