""" Classes of lesion and first ring specific of wheat septoria.

"""
# Imports #########################################################################
from alinea.alep.fungal_objects import *
# The following import would provoke a circular reference
# "from alinea.alep.septoria import SeptoriaDU"
# --> Moved in the method 'SeptoriaExchangingRings.emission()'
from alinea.alep.septoria import Disease as _Disease, SeptoriaParameters as _SeptoriaParameters
from random import randint, seed, random
from math import floor, ceil
import numpy as np
  
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
        # Counter of calculation for senescence
        self.can_compute_senescence = True
        # Age to reach to create another ring
        self.delta_age_ring = 0.
        # Age of the center of the lesion (degree days)
        self.age_tt = 0.
        # Delta degree days during time step
        self.ddday = 0.
        # Is first hour of rain
        self.first_rain_hour = False
        # Stock of spores
        self.stock_spores = None
        
        # Temporary
        self.surface_empty = 0.
        self.hist_delta_spo = []
        self.hist_stock = []
        self.nb_spores_emitted = 0.
        self.hist_inc = []
        self.hist_chlo = []
        self.hist_nec = []
        self.hist_spo = []
        self.hist_spo2 = []
        self.hist_spo3 = []
        self.hist_empty = []
        self.hist_surf = []
        self.previous_surface = 0.

    def update(self, dt, leaf=None):
        """ Update the status of the lesion and create a new growth ring if needed.
                
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. area, green area, healthy area,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        f = self.fungus
        # Compute delta degree days in dt
        self.compute_delta_ddays(dt, leaf)
        ddday = self.ddday
 
        # If senescence, compute length of growth period before senescence during time step
        # The condition on 'can_compute_senescence' allows the calculation to be made only once
        if self.is_senescent and self.can_compute_senescence==True:
            self.compute_time_before_senescence(leaf=leaf)
            ddday = self.ddday_before_senescence

        # if self.surface<0:
        if any(self.surface_rings<0):
            import pdb
            pdb.set_trace()
        assert self.surface >= 0, 'surface must be positive'
            
        # Update the age of the lesion
        self.age_tt += ddday
        
        # Temp
        if self.ddday<self.ddday_before_senescence:
            import pdb
            pdb.set_trace()
        
        try:
            if self.age_tt<self.previous_age:
                import pdb
                pdb.set_trace()
        except:
            pass
        
        
        if self.growth_is_active:
            assert self.age_tt - f.degree_days_to_chlorosis - self.delta_age_ring <= 20.

        # Compute growth demand
        if self.growth_is_active:
            self.update_growth_demand()
        
        # Manage first ring
        self.first_ring.update(lesion=self)
        
        # Manage the other rings
        if len(self.surface_rings)>0:
            if not self.growth_is_active:
                # If growth is over, suppress empty surfaces corresponding to young rings
                if self.surface_rings[0] == 0.:
                    ind = self.surface_rings.nonzero()[0][0]
                    self.surface_rings = self.surface_rings[ind:]
            # if self.age_tt - f.degree_days_to_chlorosis >= self.delta_age_ring - f.delta_age_ring:
            if round(self.age_tt - f.degree_days_to_chlorosis, 6) > round(self.delta_age_ring, 6):
                # Create a new ring when needed
                self.surface_rings = np.append(self.surface_rings, 0.)
                self.delta_age_ring += f.delta_age_ring
                
                # Attempt to create several rings in one big time step
                # Note : Failed so far...
                # diff = self.delta_age_ring - (self.age_tt - f.degree_days_to_chlorosis)
                # nb_rings = ceil(float(diff)/f.delta_age_ring)
            
            # Exchange the surfaces between the rings
            self.exchange_surfaces()

        # Update stock of spores
        if self.status==f.SPORULATING:
            self.update_stock()
            
            # Manage rain perception
            if leaf.rain_intensity > 0. and leaf.relative_humidity >= f.rh_min:
                self.first_rain_hour = True if not self.first_rain_hour else False
            else:
                self.first_rain_hour = False
        
        self.old_position_senescence = leaf.position_senescence
        
        # Temporary
        if self.status!=f.SPORULATING:
            self.hist_stock.append(0.)
            self.hist_delta_spo.append(0.)
        else:
            self.hist_stock.append(self.stock_spores)
        
        self.previous_surface=self.surface
        self.hist_surf.append(self.surface)
        
        
        try:
            self.hist_leaf_sen.append(leaf.position_senescence)
            self.hist_leaf_id.append(leaf._vid)
            self.hist_leaf_area.append(leaf.area)
        except:
            self.hist_leaf_sen=[leaf.position_senescence]
            self.hist_leaf_id = [leaf._vid]
            self.hist_leaf_area = [leaf.area]
        
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
            ddday = max(0,(leaf.temp - f.basis_for_dday*dt)/(24./dt))
        else:
            ddday = 0.
        # Save variable
        self.ddday = ddday

    def compute_time_before_senescence(self, leaf=None):
        """ Compute length of growth period before senescence during time step.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
        """
        old_pos = self.old_position_senescence
        new_pos = leaf.position_senescence
        ddday = self.ddday
        speed = (old_pos - new_pos)/ddday if ddday > 0. else 0.
        
        self.ddday_before_senescence = (old_pos-self.position[0])/speed if speed >0. else 0.
        
        # Temporary
        self.hist_leaf_id.append(leaf._vid)
        self.hist_leaf_area.append(leaf.area)
        
        if speed<0.:
            import pdb
            pdb.set_trace()
        
        if round(ddday,6)<round(self.ddday_before_senescence,6):
            import pdb
            pdb.set_trace()
        #
        
        # Reset self.old_position_senescence
        self.can_compute_senescence = False
    
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
               
        # ds = s * ddday / dring
        # Works only if ddday < 20 dd
        
        r = self.current_growth_rate
        ds = np.zeros(len(s))
        if self.surface_empty == 0. and self.surface_rings[-1]==0.:
            for ind in range(len(s)):
                ds[ind] = min(s[ind], r*(self.age_tt%dring))
        else:
            for ind in range(len(s)):
                ds[ind] = min(s[ind], r*ddday)

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
                    if self.first_ring.surface + growth_offer <= f.Smin:
                        self.first_ring.grow(self, growth_offer)
                    else:
                        # Compute sharing between first ring and following rings
                        go_first_ring = f.Smin - self.first_ring.surface
                        self.first_ring.grow(self, go_first_ring)
                        go_other_ring = growth_offer - go_first_ring
                        self.surface_rings = np.append(self.surface_rings, go_other_ring)
                        self.delta_age_ring = f.delta_age_ring
            else:
                self.surface_rings[0]+=growth_offer
                if self.surface == f.Smax:
                    self.disable_growth()                    
                
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
        
        # Temporary
        self.hist_delta_spo.append(delta_surface_spo)
            
    def compute_sporulating_surface(self):
        """ Compute only the sporulating surface.
        
        Parameters
        ----------
            None
        """
        f = self.fungus
        if self.is_sporulating():
            width_ring = f.delta_age_ring
            delta_ring = self.delta_age_ring
            age_tt = self.age_tt
            s = self.surface_rings
            time_to_chlo = f.degree_days_to_chlorosis
            time_to_spo = f.degree_days_to_necrosis + f.degree_days_to_sporulation
            
            diff = delta_ring - time_to_spo
            nb_full_rings = floor(diff/width_ring)
            if len(s)>0:
                if nb_full_rings>0:
                    portion_sporulating = (diff%width_ring)/width_ring
                    if portion_sporulating > 0:
                        if len(s)>nb_full_rings:
                            surface_spo = sum(s[-nb_full_rings:]) + portion_sporulating*s[-(nb_full_rings+1)]
                        else:
                            nb_full_rings = len(s)
                            surface_spo = sum(s[-nb_full_rings:])
                    else:
                        surface_spo = sum(s[-nb_full_rings:])
                    delta_ring -= width_ring * nb_full_rings
                else:
                    # portion_sporulating = (age_tt-time_to_spo)/(age_tt+width_ring-delta_ring)
                    portion_sporulating = (age_tt-time_to_chlo-time_to_spo)/width_ring
                    surface_spo = portion_sporulating*s[-(nb_full_rings+1)]
            # self.surface_spo = surface_spo + self.first_ring.surface
            else:
                surface_spo = 0.
            # Temporary
            self.surface_spo = surface_spo + self.first_ring.surface - self.surface_empty
        else:
            self.surface_spo = 0.

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
        width_ring = f.delta_age_ring
        delta_ring = self.delta_age_ring
        s = np.copy(self.surface_rings)
        time_to_chlo = f.degree_days_to_chlorosis
        time_to_nec = f.degree_days_to_necrosis
        time_to_spo = f.degree_days_to_necrosis + f.degree_days_to_sporulation
        age_tt = self.age_tt

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
            diff = delta_ring - time_to_nec
            nb_full_rings = floor(diff/width_ring)
            if len(s)>0:
                if nb_full_rings>0:
                    portion_necrotic = (diff%width_ring)/width_ring
                    if portion_necrotic > 0:
                        if len(s)>nb_full_rings:
                            surface_nec = sum(s[-nb_full_rings:]) + portion_necrotic*s[-(nb_full_rings+1)]
                            s[-(nb_full_rings+1)] *= (1-portion_necrotic)
                        else:
                            nb_full_rings = len(s)
                            surface_nec = sum(s[-nb_full_rings:])
                    else:
                        surface_nec = sum(s[-nb_full_rings:])
                    s[-nb_full_rings:] = 0.
                else:
                    # portion_necrotic = (age_tt-time_to_nec)/(age_tt+width_ring-delta_ring)
                    portion_necrotic = (age_tt-time_to_chlo-time_to_nec)/width_ring
                    surface_nec = portion_necrotic*s[-1]
                    s[-1] *= (1-portion_necrotic)
            else:
                surface_nec=0.
            surface_nec += self.first_ring.surface
            # Compute chlorotic surface
            surface_chlo = sum(s)
        
        elif status == f.SPORULATING:
            # Compute sporulating surface
            diff = delta_ring - time_to_spo
            nb_full_rings = floor(diff/width_ring)
            # Temp
            dring = delta_ring

            if len(s)>0:
                if nb_full_rings>0:
                    portion_sporulating = (diff%width_ring)/width_ring
                    if portion_sporulating > 0:
                        if len(s)>nb_full_rings:
                            surface_spo = sum(s[-nb_full_rings:]) + portion_sporulating*s[-(nb_full_rings+1)]
                            s[-(nb_full_rings+1)] *= (1-portion_sporulating)
                        else:
                            nb_full_rings = len(s)
                            surface_spo = sum(s[-nb_full_rings:])
                    else:
                        surface_spo = sum(s[-nb_full_rings:])
                    s[-nb_full_rings:] = 0.
                    delta_ring -= width_ring * nb_full_rings
                    s = s[:-nb_full_rings]
                else:
                    # portion_sporulating = (age_tt-time_to_spo)/(age_tt+width_ring-delta_ring)
                    portion_sporulating = (age_tt-time_to_chlo-time_to_spo)/width_ring
                    surface_spo = portion_sporulating*s[-1]
                    
                    # Temp
                    try:
                        if round(surface_spo+self.first_ring.surface,6)<round(self.previous_spo,6):
                            import pdb
                            pdb.set_trace()
                    except:
                        pass
                        
                    
                    s[-1] *= (1-portion_sporulating)
            else:
                surface_spo=0.
            # Temporary
            portion_sporulating = 0.
            surface_spo += self.first_ring.surface
            
            # # if surface_spo<self.surface_empty:
            # try:
                # if round(self.previous_spo, 6)>round(surface_spo, 6):
                # # if self.previous_spo>surface_spo:
                    # import pdb
                    # pdb.set_trace()
            # except:
                # pass
                    
            # Temporary
            # try:
                # if round(surface_spo,4)<round(self.previous_spo,4):
                    # import pdb
                    # pdb.set_trace()
            # except:
                # pass
            
            self.previous_rings = self.surface_rings
            self.previous_first = self.first_ring.surface
            self.previous_age = self.age_tt
            self.previous_dring = dring
            self.previous_spo = surface_spo
            self.previous_full = nb_full_rings
            self.previous_portion = portion_sporulating

            # Compute necrotic surface
            diff = delta_ring - time_to_nec
            nb_full_rings = floor(diff/width_ring)
            if len(s)>0:
                if nb_full_rings>0:
                    portion_necrotic = (diff%width_ring)/width_ring
                    if portion_necrotic > 0:
                        if len(s)>nb_full_rings:
                            surface_nec = sum(s[-nb_full_rings:]) + portion_necrotic*s[-(nb_full_rings+1)]
                            s[-(nb_full_rings+1)] *= (1-portion_necrotic)
                        else:
                            nb_full_rings = len(s)
                            surface_nec = sum(s[-nb_full_rings:])
                    else:
                        surface_nec = sum(s[-nb_full_rings:])
                    s[-nb_full_rings:] = 0.
                else:
                    # portion_necrotic = (age_tt-time_to_nec)/(age_tt+width_ring-delta_ring)
                    portion_necrotic = (age_tt-time_to_chlo-time_to_nec)/width_ring
                    surface_nec = portion_necrotic*s[-1]
                    s[-1] *= (1-portion_necrotic)
            # Compute chlorotic surface
            surface_chlo = sum(s)

        # Save variables
        self.surface_inc = surface_inc
        self.surface_chlo = surface_chlo
        self.surface_nec = surface_nec
        # self.surface_spo = surface_spo
        # Temporary
        # if self.surface_empty<0.:
        if surface_spo<0.:
            import pdb
            pdb.set_trace()
        
        self.surface_spo = surface_spo - self.surface_empty
        
        if self.surface_spo<0.:
            import pdb
            pdb.set_trace()
        
        self.hist_inc.append(surface_inc)
        self.hist_chlo.append(surface_chlo)
        self.hist_nec.append(surface_nec)
        self.hist_spo.append(surface_spo)
        self.hist_spo2.append(self.surface_spo)
        self.hist_empty.append(self.surface_empty)

    def emission(self, leaf=None):
        """ Create a list of dispersal units emitted by the lesion.
        
        Parameters
        ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. area, green area, healthy area,
            senescence, rain intensity, wetness, temperature, lesions)
        
        .. Todo:: Implement a real formalism.
        """
        # Import to break circular reference
        from alinea.alep.septoria import SeptoriaDU
        # TODO : Improve ?
        
        import pdb
        pdb.set_trace()
        print('emission is computed with an external model now')
        
        if self.is_stock_available(leaf):
            f = self.fungus
            emissions = []
            stock_available = int(self.stock_spores*2/3.)
            
            # Temporary
            initial_stock = self.stock_spores
            
            # Temporary REVISION 25/10/2013 G.Garin:
            tot_surf_spo = 0.
            for les in leaf.lesions:
                les.compute_sporulating_surface()
                tot_surf_spo += les.surface_spo
            contribution = self.surface_spo/tot_surf_spo if tot_surf_spo>0. else 0.
            tot_fraction_spo = tot_surf_spo/leaf.area if leaf.area>0. else 0.
            total_DU_emitted = 0.36 * 6.19e7 * tot_fraction_spo * leaf.rain_intensity
            nb_DU_emitted = int(contribution * total_DU_emitted)
            
            # TODO : improve below
            # nb_DU_emitted = int(leaf.rain_intensity * self.surface_spo * 1000)
            nb_spores_by_DU = []
            for DU in range(nb_DU_emitted):
                if stock_available > 0.:
                    # nb_spores = min(randint(5,100), stock_available)
                    nb_spores = min(1, stock_available)
                    nb_spores_by_DU.append(nb_spores)
                    stock_available -= nb_spores
                    # Update stock_spores
                    self.stock_spores -= nb_spores
                    
            # Temporary : Test for calculation of empty surface
            self.update_empty_surface(initial_stock)

            # Get rid of DUs without spores
            nb_DU_emitted = len(nb_spores_by_DU)
            
            # Empty stock_spores
            if self.stock_spores < 1000:
                self.stock_spores = 0.

            # Return emissions
            SeptoriaDU.fungus = f
            emissions = [SeptoriaDU(nb_spores = nb_spores_by_DU[i], status='emitted', 
                            position=self.position) for i in range(nb_DU_emitted)]

            return emissions
        else:
            return []
    
    def reduce_stock(self, nb_spores_emitted):
        """ Reduce the stock of spores after emission.
        
        Parameters
        ----------
        nb_spores_emitted: int
            Number of spores emitted
        """
        self.stock_spores -= nb_spores_emitted
        if self.stock_spores < self.fungus.threshold_spores:
            self.stock_spores = 0.
    
    def update_empty_surface(self, nb_spores_emitted, initial_stock):  
        """ Update empty surface after emission. 
        
        In this case, the surface is emptied proportionally 
        to the number of spores emitted.
        
        Parameters
        ----------
        nb_spores_emitted: int
            Number of spores emitted
        initial_stock: int
            Number of spores on the lesion before emission
        """
        assert initial_stock>0.
        
        self.compute_all_surfaces()
       
        # Temporary
        surface_spo = self.surface_spo
        
        surface_empty = (nb_spores_emitted/initial_stock) * self.surface_spo

        # Temporary
        if surface_empty>surface_spo:
            import pdb
            pdb.set_trace()
            
        self.surface_empty += surface_empty

        # Temporary
        self.nb_spores_emitted += nb_spores_emitted
        
    def is_sporulating(self):
        """ Check if the lesion is sporulating.
        """
        return self.status == self.fungus.SPORULATING
    
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
            A leaf sector with properties (e.g. area, green area, healthy area,
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
    
    def become_senescent(self):
        """ The lesion will become senescent during this time step.
        """
        # Turn on 'is_senescent'
        self.is_senescent = True

    def senescence_response(self):
        """ Compute surface alive and surface dead.
        
        Parameters
        ----------
            None
        """
        # Runs only for the first time step of senescence occurence
        if self.ddday_before_senescence!=None:
            
            f = self.fungus
            ddday = self.ddday
            ddday_sen = self.ddday_before_senescence
            time_to_nec = f.degree_days_to_necrosis
            dring = f.delta_age_ring
            s = np.copy(self.surface_rings)
            
            # Stop growth
            self.disable_growth()
            
            # Kill all surfaces under necrosis at senescence occurence
            if self.status <= f.CHLOROTIC:
                self.disable()
                self.first_ring.disable()
                self.surface_dead = self.surface_alive
                self.surface_alive = 0.
                self.surface_rings = np.array([])
            else:
                # Update the list of rings
                if len(s)>0 and self.surface_dead==0.:
                    diff = self.delta_age_ring - time_to_nec
                    nb_full_rings = floor(diff/dring)
                    portion_last_ring = (diff%dring)/dring
                    if len(s)>nb_full_rings:
                        surface_dead = sum(s[:-(nb_full_rings+1)])+(1-portion_last_ring)*s[-(nb_full_rings+1)]
                        s[:-(nb_full_rings+1)] = 0.
                        s[-(nb_full_rings+1)] *= portion_last_ring
                    else:
                        surface_dead = 0.
                    self.surface_rings = s
                    self.surface_alive = sum(self.surface_rings) + self.first_ring.surface
                    self.surface_dead = surface_dead
            
            if any(self.surface_rings<0):
                import pdb
                pdb.set_trace()
            
            # Complete the age of the lesion up to the end of time step
            # self.ddday = ddday - ddday_sen
            self.ddday = ddday_sen
            self.age_tt += ddday - ddday_sen
            
            # Manage first ring
            self.first_ring.update(lesion=self)
            # Manage the other rings
            if len(self.surface_rings)>0:
                if self.age_tt-f.degree_days_to_chlorosis > self.delta_age_ring:
                    self.surface_rings = np.append(self.surface_rings, 0.)
                    self.delta_age_ring += f.delta_age_ring
                self.exchange_surfaces()
                
            # Update stock of spores
            if self.status==f.SPORULATING:
                self.update_stock()
            
            # Reset self.ddday_before_senescence
            self.ddday_before_senescence=None

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
        return self.surface_nec + self.surface_spo + self.surface_empty
    
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
            
# First ring ######################################################################
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
        self.age_tt = 0.
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
        self.age_tt = lesion.age_tt
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
        if self.age_tt >= time_to_chlorosis:
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
        if self.age_tt >= time_to_necrosis:
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
        if self.age_tt >= time_to_spo:
            self.status = f.SPORULATING
    
    def grow(self, lesion, surface):
        f = lesion.fungus
        self.surface += surface
        if self.surface == f.Smin:
            self.disable_growth()
        assert self.surface <= f.Smin

    def disable_growth(self):
        """ Shut down ring growth activity (turn it to False)
        
        Parameters
        ----------
            None
        """
        self.growth_is_active = False
        self.growth_demand = 0.
    
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
        _SeptoriaParameters.__init__(self, model=SeptoriaExchangingRings, **kwds)
        
class Disease(_Disease):
    @classmethod
    def parameters(cls, **kwds):
        return Parameters(**kwds)
    
    @classmethod
    def lesion(cls, **kwds):
        SeptoriaExchangingRings.fungus=cls.parameters(**kwds)
        return SeptoriaExchangingRings
 
