""" Class of lesion of wheat septoria. 
    Progress from one state to another is computed with physiological age.
"""

# Imports #########################################################################
from alinea.alep.fungal_objects import *
from alinea.alep.septoria import SeptoriaDU, septoria_parameters, is_iterable
import numpy as np
from math import floor, ceil

# Lesion ##########################################################################
class SeptoriaAgePhysio(Lesion):
    """ Lesion of septoria implemented with growth stages that exchange surfaces
        according to their physiological age.
    """
    def __init__(self, mutable=False):
        """ Initialize the lesion of septoria. 
        
        Parameters
        ----------
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        position: non defined
            Position of the dispersal unit on the phyto-element
        """
        super(SeptoriaAgePhysio, self).__init__(mutable=False)
        # Calculate parameter for emission from other parameters
        self.fungus.density_dus_emitted_max = np.array([self.fungus.density_dus_emitted_ref*(self.fungus.reduction_by_rain**ind) for ind in range(self.fungus.rain_events_to_empty)])
        # Status of the center of the lesion
        self.status = self.fungus.INCUBATING
        # Age of the center of the lesion
        self.ddday=None
        self.age_t = 0.
        self.age_dday = 0.
        self.age_physio = 0.
        # Status of the periphery of the lesion
        self.status_edge = self.fungus.INCUBATING
        # Age of the periphery of the lesion
        self.age_physio_edge = 0.
        # Ratio left to progress in new status when center of the lesion changes status during time step
        self.ratio_left = 0.
        self.ratio_left_edge = 0.
        # Surfaces exchanged
        self.to_necrosis = 0.
        self.to_sporulation = 0.
        # Factor for sharing surfaces in new forming rings (number of rings to fill)
        self.distribution_new_rings = 0.
        # Rings in each state
        self.surface_first_ring = 0.
        self.surfaces_chlo = np.array([])
        self.surfaces_nec = np.array([])
        self.surfaces_spo = np.zeros(self.fungus.rain_events_to_empty)
        self.surface_empty = 0.
        self.surface_senesced = 0.
        self.incubation_completed = False
        # Surface of disabled rings
        self.surface_dead = 0.
        # Stock of spores
        self.stock_spores = None
        self.nb_spores_emitted = None
        # Is first hour of rain
        # self.first_rain_hour = False
        # Counter of calculation for senescence
        # self.can_compute_senescence = True
        self.senescence_response_completed = False
        # Old position senescence
        # self.old_position_senescence = None
        # dt left after senescence
        # self.dt_before_senescence = None
        # self.dt_left_after_senescence = None
    
        # Temporary
        # self.previous_surface = 0.
        # self.hist_age = []
        # self.hist_inc = []
        self.hist_chlo = []
        # self.hist_nec = []
        # self.hist_spo = []
        # self.hist_empty = []
        # self.hist_can_compute_rings = []
        self.hist_nb_rings_chlo = []
        self.hist_progress = []
        self.hist_age_physio = []
    
    def is_incubating(self):
        """ Check if lesion status is incubation. """
        return self.status == self.fungus.INCUBATING
    
    def is_chlorotic(self):
        """ Check if lesion status is chlorosis. """
        return self.status == self.fungus.CHLOROTIC
        
    def is_necrotic(self):
        """ Check if lesion status is necrosis. """
        return self.status == self.fungus.NECROTIC
    
    def is_sporulating(self):
        """ Check if lesion status is sporulation. """
        return self.status == self.fungus.SPORULATING
    
    def update(self, dt, leaf=None):
        """ Update the status of the lesion and create a new growth ring if needed.
                
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. area, green area, healthy area,
            senescence, rain intensity, wetness, temperature, lesions, etc.)
        """
        # Manage senescence
        # if self.is_senescent and self.can_compute_senescence==True:
            # dt = self.compute_time_before_senescence(dt, leaf=leaf)
        # self.old_position_senescence = leaf.position_senescence
        if any([x[0]<leaf.senesced_length for x in self.position]):
            self.senescence_response(leaf.senesced_length)
        
        if self.is_active:
            
            # Compute delta degree days in dt
            self.compute_delta_ddays(dt, leaf)
            
            if self.ddday > 0.:
                # Update age in degree days of the lesion
                self.age_dday += self.ddday        
                # Update growth demand and status
                self.update_status()
            else:
                self.growth_demand = 0.
            
            # Temporary
            # self.hist_age.append(self.age_dday)
            # self.hist_inc.append(self.surface_inc)
            self.hist_chlo.append(self.surface_chlo)
            # self.hist_nec.append(self.surface_nec)
            # self.hist_spo.append(self.surface_spo)
            # self.hist_empty.append(self.surface_empty)
            # self.hist_can_compute_rings.append(self.can_compute_rings())
            self.hist_nb_rings_chlo.append(len(self.surfaces_chlo))
            self.hist_age_physio.append(self.age_physio)
            
            # Temporary      
            # if self.surface_chlo>0. and len(self.surfaces_nec>1.) and self.surfaces_nec[0]==0.:
                # import pdb
                # pdb.set_trace()

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
        if self.ddday!=None:
            f = self.fungus
            # Calculation
            if dt != 0.:
                ddday = sum([max(0,(temp - f.basis_for_dday)*1/24.) for temp in leaf.temperature_sequence])
                if 'global_efficacy' in leaf.properties():
                    ddday *= (1 - max(0, min(1, leaf.global_efficacy['eradicant'])))
            else:
                ddday = 0.
            # Save variable
            self.ddday = ddday
        else:
            self.ddday = 0.
    
    def progress(self, age_threshold=0.):
        """ Compute progress in physiological age according to age_threshold. 
        """
        left = self.ratio_left
        if left==0.:
            # Normal situation
            progress = self.ddday/age_threshold if age_threshold>0. else 0.
        elif left>0.:
            # Passing from one state to another in same time step
            assert left<1.
            progress = left*self.ddday/age_threshold if age_threshold>0. else 0.
            # Reset ratio left and age physio
            self.ratio_left = 0.
        return progress
    
    def control_growth(self, growth_offer=0.):
        """ Reduce surface of the rings up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Surface available on the leaf for the ring to grow (cm2)
        """
        if self.growth_is_active:
            
            # Disable growth in case of competition 
            if round(growth_offer,14) < round(self.growth_demand,14):
                self.disable_growth()
                # TEMPORARY
                self.age_competition = self.age_dday

            f = self.fungus
            # Growth offer is added to surface according to state
            if self.is_incubating():
                if growth_offer<0:
                    import pdb
                    pdb.set_trace()
                self.surface_first_ring += growth_offer
                if round(self.surface_first_ring,14) == round(self.fungus.Smin * self.nb_lesions, 14) and self.incubation_completed==False:
                    self.incubation_completed = True
            else:
                Smin = self._surface_min
                if round(self.surface_first_ring,14) < round(Smin,14):
                    if self.surface_first_ring + growth_offer <= Smin:
                        self.surface_first_ring += growth_offer
                        growth_offer = 0.
                    else:
                        diff = Smin - self.surface_first_ring
                        self.surface_first_ring = Smin
                        growth_offer -= diff
                    if round(self.surface_first_ring,14) == round(self._surface_min, 14) and self.incubation_completed==False:
                        self.incubation_completed = True
                
                nb_full_rings = int(floor(self.distribution_new_rings))
                surf = np.array([])
                for rg in range(nb_full_rings):
                    filling = growth_offer/self.distribution_new_rings
                    growth_offer-=filling
                    surf = np.append(surf, filling)
                surf = np.append(surf, round(growth_offer, 14))
                # Fill surfaces chlo
                if len(surf)<=len(self.surfaces_chlo):
                    self.surfaces_chlo[:len(surf)]+=surf
                else:
                    nb_existing = len(self.surfaces_chlo)
                    nb_to_create = len(surf) - len(self.surfaces_chlo)
                    self.surfaces_chlo += surf[:len(self.surfaces_chlo)]
                    self.surfaces_chlo = np.append(self.surfaces_chlo, surf[len(self.surfaces_chlo):])
                
                if self.is_sporulating() and self.status_edge==f.CHLOROTIC and self.age_physio_edge==0. and len(self.surfaces_chlo)<f.nb_rings_by_state:
                    import pdb
                    pdb.set_trace()

            # Reset distribution in new rings and growth demand
            self.distribution_new_ring = 0.         

            # If lesion has reached max size, disable growth
            if round(self.surface,14) >= round(self._surface_max,14):
                self.disable_growth()  

            self.growth_demand = 0.
    
    def incubation(self):
        """ Compute growth demand and physiological age progress to chlorosis.
        """
        f = self.fungus
        time_to_chlo = f.degree_days_to_chlorosis
        # Compute progress in incubation
        progress = self.progress(age_threshold=time_to_chlo)
        self.age_physio += progress
        # Compute growth demand
        if self.age_physio<1.:
            if self.surface_first_ring<0:
                import pdb
                pdb.set_trace()
            if self.growth_is_active:
                self.growth_demand = progress * f.Smin * self.nb_lesions
        else:
            if self.growth_is_active:
                Smin = self._surface_min
                self.growth_demand = round(Smin - self.surface_alive,14)
                    
            if self.growth_demand < 0:
                import pdb
                pdb.set_trace()
            
            # Change status
            self.ratio_left = round((self.age_physio - 1.)/progress,14)
            self.change_status()
            self.change_status_edge()
            self.reset_age_physio()
            self.chlorosis()

    def chlorosis(self):
        """ Compute growth demand and physiological age progress to necrosis.
        """
        f = self.fungus
        # Compute progress in chlorosis
        time_to_nec = f.degree_days_to_necrosis
        progress = self.progress(age_threshold=time_to_nec)
        
        # Temp debug
        self.hist_progress.append(progress)
        # self.hist_can_compute_rings.append(self.can_compute_rings())
                
        # Compute growth demand
        if self.growth_is_active:
            # Note : '+=' because might be added to growth demand in incubation 
            # if transition in same time step; added to zero otherwise.
            self.growth_demand += progress * time_to_nec * f.growth_rate * self.nb_lesions
            # Limit growth to size max
            Smax = self._surface_max
            if self.surface + self.growth_demand >= Smax:
                self.growth_demand = Smax - self.surface
                
        # Compute exchanges of surfaces
        if self.incubation_completed:
            age_physio = self.age_physio
            age_edge = self.age_physio_edge
            default_nb_rings = f.nb_rings_by_state
            # 'rings' is the division of the stage in classes of age
            # 'width' is the width in physiologic age of each class
            (rings, width) = np.linspace(0,1, default_nb_rings+1, retstep=True)
            # Calculate the number of rings in which chlorosis input will be shared
            self.distribution_new_rings = progress/width        
            if self.is_chlorotic() and age_physio==0.:
                pass
                # (No exchange of surface the first time step the lesion enters stage)
            elif len(self.surfaces_chlo)>0:
                # Calculate exchanges of surfaces between rings 
                # 1. Reduce the superior limit of ring ages if age_physio in chlorosis
                if self.is_chlorotic():
                    rings = rings[:ceil(age_physio/width)+1]
                    rings[-1] = age_physio
                # 2. Reduce the inferior limit of ring ages if age_edge in chlorosis
                if self.status_edge==f.CHLOROTIC:
                    rings = rings[floor(age_edge/width):]
                    rings[0] = age_edge
                
                # 3. Get the beginnings and the ends of age classes
                begs = rings[:-1]
                ends = rings[1:]
                # 4. Apply progress to the beginnings and the ends of age classes
                begs_prog = begs + progress
                ends_prog = ends + progress
                # 5. Find ends of new classes in which surfaces will be distributed after progress
                new_ends = np.arange(max(0.1, width*ceil(begs_prog[0]/width)), width*(ceil(ends_prog[-1]/width)+1), width)
                new_ends = np.round(new_ends, 14)
                
                # 6. Loop over this classes to calculate the new distribution in each class
                new_surf = np.zeros(len(new_ends))              
                for j in range(len(new_ends)):
                    # Find beginnings and ends of rings in new class
                    ind_begs = np.where((new_ends[j]-width <= begs_prog) * (begs_prog < new_ends[j]))[0]
                    ind_ends = np.where((new_ends[j]-width <= ends_prog) * (ends_prog < new_ends[j]))[0]
                    indexes = np.unique(np.append(ind_begs, ind_ends))
                    # Calculate ratio of old class in new class
                    for ind in indexes:
                        new_surf[j] += round((min(ends_prog[ind], new_ends[j])-
                                              max(begs_prog[ind],new_ends[j]-width))*
                                              self.surfaces_chlo[ind]/(ends[ind]-begs[ind]), 14)
                # 7. Update surfaces and calculate what passes to necrosis
                self.surfaces_chlo = np.extract(new_ends<=1, new_surf)
                self.to_necrosis = sum(np.extract(new_ends>1, new_surf))
        
        # Ageing of the periphery of the lesion if growth has been stopped
        if self.status_edge==f.CHLOROTIC and not self.growth_is_active :
            if self.age_physio_edge+progress < 1.:
                self.age_physio_edge += progress
            else:
                diff = self.age_physio_edge + progress - 1.
                self.ratio_left_edge = diff/progress
                self.change_status_edge()
                self.reset_age_physio_edge()
        
        # Ageing of the center of the lesion
        if self.is_chlorotic():
            if self.age_physio+progress < 1.:
                self.age_physio += progress
            else:
                diff = self.age_physio + progress - 1.
                self.ratio_left = diff/progress
                self.change_status()
                self.reset_age_physio()
                self.necrosis()
                    
    def necrosis(self):
        """ Compute physiological age progress to sporulation.
        """
        f = self.fungus
        # Compute progress in necrosis
        time_to_spo = f.degree_days_to_sporulation
        progress = self.progress(age_threshold=time_to_spo)
        
        # Compute exchanges of surfaces
        if self.incubation_completed:
            age_physio = self.age_physio
            age_edge = self.age_physio_edge
            default_nb_rings = f.nb_rings_by_state
            (rings, width) = np.linspace(0,1, default_nb_rings+1, retstep=True)
            if self.is_necrotic() and age_physio==0.:
                pass
            elif len(self.surfaces_nec)>0:
                if self.is_necrotic():
                    rings = rings[:ceil(age_physio/width)+1]
                    rings[-1] = self.age_physio
                if self.status_edge==f.NECROTIC:
                    rings = rings[floor(age_edge/width):]
                    rings[0] = self.age_physio_edge
                begs = rings[:-1]
                ends = rings[1:]
                begs_prog = begs + progress
                ends_prog = ends + progress          
                # new_classes = np.arange(width*ceil(begs_prog[0]/width),
                                        # width*(ceil(ends_prog[-1]/width)+1),
                                        # width)
                # new_classes = np.round(new_classes, 14)
                new_ends = np.arange(max(0.1, width*ceil(begs_prog[0]/width)), width*(ceil(ends_prog[-1]/width)+1), width)
                new_ends = np.round(new_ends, 14)
                
                new_surf = np.zeros(len(new_ends))
                for j in range(len(new_ends)):
                    # Find beginnings and ends of rings in class
                    ind_begs = np.where((new_ends[j]-width <= begs_prog) * (begs_prog < new_ends[j]))[0]
                    ind_ends = np.where((new_ends[j]-width <= ends_prog) * (ends_prog < new_ends[j]))[0]
                    indexes = np.unique(np.append(ind_begs, ind_ends))
                    for ind in indexes:
                        new_surf[j] += round((min(ends_prog[ind], new_ends[j])-
                                              max(begs_prog[ind],new_ends[j]-width))*
                                              self.surfaces_nec[ind]/(ends[ind]-begs[ind]), 14)
                                    
                # Get what passes to next status
                self.surfaces_nec = np.extract(new_ends<=1, new_surf)
                self.to_sporulation = sum(np.extract(new_ends>1, new_surf))
            
            # Filling of new rings
            nb_full_rings = int(floor(progress/width))
            surf = np.array([])
            for rg in range(nb_full_rings):
                filling = self.to_necrosis*width/progress if progress > 0. else 0.
                self.to_necrosis -= filling
                surf = np.append(surf, filling)
            surf = np.append(surf, self.to_necrosis)
            if len(surf)<=len(self.surfaces_nec):
                self.surfaces_nec[:len(surf)]+=surf
            else:
                nb_existing = len(self.surfaces_nec)
                nb_to_create = len(surf) - len(self.surfaces_nec)
                self.surfaces_nec += surf[:len(self.surfaces_nec)]
                self.surfaces_nec = np.append(self.surfaces_nec, surf[len(self.surfaces_nec):])
            self.to_necrosis = 0.
        
        # Ageing of the periphery of the lesion if growth has been stopped       
        if self.status_edge==f.NECROTIC and not self.growth_is_active:
            if self.age_physio_edge==0. and self.ratio_left_edge>0.:
                self.age_physio_edge += progress*self.ratio_left_edge
                self.ratio_left_edge = 0.
            elif self.age_physio_edge+progress < 1.:
                self.age_physio_edge += progress
            else:
                diff = self.age_physio_edge + progress - 1.
                self.ratio_left_edge = diff/progress
                self.change_status_edge()
    
        # Ageing of the center of the lesion
        if self.is_necrotic():
            if self.age_physio+progress < 1.:
                self.age_physio += progress
            else:
                diff = self.age_physio + progress - 1.
                self.ratio_left = diff/progress
                self.change_status()
                self.reset_age_physio()
                self.sporulation()
     
    def sporulation(self):
        """ Compute production of spores. """
        f = self.fungus
        # First time first ring sporulates
        # if self.stock_spores==None:
            # self.stock_spores = self.surface_first_ring * f.production_rate
            # self.surface_spo += self.surface_first_ring
        # self.stock_spores += self.to_sporulation * f.production_rate
        # self.surface_spo += self.to_sporulation
        # self.to_sporulation = 0.
        
        if self.stock_spores==None:
            self.stock_spores = self.surface_first_ring * f.production_rate
            self.surfaces_spo[0] += self.surface_first_ring
        self.stock_spores += self.to_sporulation * f.production_rate
        self.surfaces_spo[0] += self.to_sporulation
        self.to_sporulation = 0.
    
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
        # f = self.fungus
        # return (self.stock_spores>0. and self.is_sporulating() and 
                # leaf.relative_humidity >= f.rh_min and self.first_rain_hour)
        return self.is_sporulating()

    def emission(self, density_DU_emitted):
        """ Return number of DUs emitted by the lesion. """
        if density_DU_emitted>0:
            f = self.fungus
            # print 'density DU popDrops %f' % density_DU_emitted
            du_factor = [float(density_DU_emitted)/dmax if dmax>0. else 0. for dmax in f.density_dus_emitted_max]
            du_factor[du_factor>1.]=1.
            delta_spo = self.surfaces_spo*du_factor
            self.surfaces_spo -= delta_spo
            self.surfaces_spo[1:]+=delta_spo[:f.rain_events_to_empty-1]
            if delta_spo[-1]>0.:
                self.surface_empty += delta_spo[-1]
            if self.status_edge == f.SPORULATING and round(self.surface_spo,14) <= round(f.threshold_spo * self.nb_lesions,14):
                # Everything becomes empty and the lesion is disabled
                self.surface_empty += sum(self.surfaces_spo)
                self.surfaces_spo = np.zeros(len(self.surfaces_spo))
                self.change_status()
                self.disable()
            return [f.dispersal_unit() for i in range(int(sum(delta_spo*[min(x,density_DU_emitted) for x in f.density_dus_emitted_max])))]
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
        
        In this case, the surface is emptied proportionally  to the number of spores emitted.
        
        Parameters
        ----------
        nb_spores_emitted: int
            Number of spores emitted
        initial_stock: int
            Number of spores on the lesion before emission
        """
        assert initial_stock>0.
        f = self.fungus
        new_surface_empty = (nb_spores_emitted/initial_stock) * self.surface_spo
        self.surface_empty += new_surface_empty
        self.surface_spo = max(0., self.surface_spo-new_surface_empty)
        self.nb_spores_emitted += nb_spores_emitted
        if self.status_edge == f.SPORULATING and self.surface_spo == 0.:
            self.change_status()
            self.disable()
    
    def compute_time_before_senescence(self, dt=1., leaf=None):
        """ Compute portion of time step before senescence. """
        old_pos = self.old_position_senescence
        new_pos = leaf.position_senescence
        speed = (old_pos - new_pos)/dt if dt > 0. else 0.
        new_dt = (old_pos-self.position[0])/speed if speed >0. else 0.
        self.dt_before_senescence = dt
        self.dt_left_after_senescence = dt - new_dt
        return new_dt
    
    def senescence_response(self, senesced_length):
        """ Compute surface alive and surface dead after senescence. """
        if not self.is_senescent:
            self.become_senescent()
            
        if not self.senescence_response_completed:
        
            nb_sen = len(filter(lambda x: x[0]<senesced_length, self.position))
            ratio_sen = float(nb_sen)/self.nb_lesions
            
            # Save surface of senesced lesions in cohort
            self.surface_senesced = self.surface * ratio_sen
            
            # Reduce surfaces alive
            f = self.fungus
            age_switch = f.age_physio_switch_senescence
            age_physio = self.age_physio
            age_edge = self.age_physio_edge
            if self.status < f.CHLOROTIC:
                self.surface_dead += round(self.surface_first_ring * ratio_sen,14)
                self.surface_first_ring = round(self.surface_first_ring * (1-ratio_sen),14)
            elif self.is_chlorotic() and age_switch >= age_physio:
                self.surface_dead += round(self.surface_first_ring * ratio_sen + sum(self.surfaces_chlo * ratio_sen),14)
                self.surface_first_ring = round(self.surface_first_ring * (1-ratio_sen),14)
                self.surfaces_chlo *= (1-ratio_sen)
                self.surfaces_chlo = self.surfaces_chlo[self.surfaces_chlo>0.]
            elif self.status_edge == f.CHLOROTIC and age_switch > age_edge:
                default_nb_rings = f.nb_rings_by_state
                (rings, width) = np.linspace(0,1, default_nb_rings+1, retstep=True)
                if self.is_chlorotic():
                    rings = rings[:ceil(age_physio/width)+1]
                    rings[-1] = age_physio
                rings = rings[floor(age_edge/width):]
                rings[0] = age_edge
                if age_switch>rings[0]:
                    self.surface_dead += sum(np.extract(rings<age_switch-width, self.surfaces_chlo * ratio_sen))
                    self.surface_dead += self.surfaces_chlo[np.where(rings<age_switch)[0][-1]]*ratio_sen*(round(age_switch%width, 14)/width if 
                                                                            (round(age_switch%width, 14)/width)>0. else 1.)                    
                    self.surfaces_chlo[:np.where(rings<age_switch)[0][-1]] *= (1 - ratio_sen)
                    self.surfaces_chlo[np.where(rings<age_switch)[0][-1]]*= (1 - ratio_sen)*((1-round(age_switch%width, 14)/width) if 
                                                                            (1-round(age_switch%width, 14)/width)>0. else 1.)
                    if nb_sen==self.nb_lesions:
                        if age_switch==1 or len(self.surfaces_chlo>0.)==0:
                            self.age_physio_edge = 0.
                            self.change_status_edge()
                        else:
                            self.age_physio_edge = age_switch
                    elif any([x==0. for x in self.surfaces_chlo]):
                        rings = rings[self.surfaces_chlo>0.]
                        if len(rings)>0:
                            self.age_physio_edge = rings[0]
                        else:
                            self.age_physio_edge = 0.
                            self.change_status_edge()
                    self.surfaces_chlo = self.surfaces_chlo[self.surfaces_chlo>0.]
                else:
                    self.senescence_response_completed==True

            self.position = filter(lambda x: x[0]>senesced_length, self.position)

            # Stop growth when last lesion of cohort is reached
            if self.nb_lesions==0.:
                if self.is_incubating() or (self.is_chlorotic() and age_switch >= age_physio):
                    self.senescence_response_completed==True
                    self.disable()
                else:
                    self.disable_growth()
                    self.senescence_response_completed==True
            
        
        # if self.can_compute_senescence:
            
            # # Reduce number of lesions in cohort
            
            # # Stop growth when last lesion of cohort is reached
            # if len(self.nb_lesions)==1.:
                # self.can_compute_senescence = False
                # self.disable_growth()
        
            # # Calculate surface dead according to status
            # f = self.fungus
            # age_switch = f.age_physio_switch_senescence
            # age_physio = self.age_physio
            # age_edge = self.age_physio_edge
                        
            # if self.status < f.CHLOROTIC:
                # self.surface_dead = self.surface_alive
                # if self.surface_dead<0.:
                    # import pdb
                    # pdb.set_trace()
                # self.surface_first_ring = 0.
                # self.disable()
            # elif self.is_chlorotic() and age_switch >= age_physio:
                # self.surface_dead = self.surface_alive
                # self.surface_first_ring = 0.
                # self.surfaces_chlo = np.array([])
                # self.disable()
            # elif self.status_edge == f.CHLOROTIC and age_switch > age_edge:
                # default_nb_rings = f.nb_rings_by_state
                # (rings, width) = np.linspace(0,1, default_nb_rings+1, retstep=True)
                # if self.is_chlorotic():
                    # rings = rings[:ceil(age_physio/width)+1]
                    # rings[-1] = age_physio
                # rings = rings[floor(age_edge/width):]
                # rings[0] = age_edge
                # self.surface_dead = sum(np.extract(rings<age_switch-width, self.surfaces_chlo))
                # self.surfaces_chlo = self.surfaces_chlo[np.where(rings<age_switch)[0][-1]:]
                # self.surface_dead += self.surfaces_chlo[0]*round(age_switch%width, 14)/width
                # self.surfaces_chlo[0]*=(1-round(age_switch%width, 14)/width)
                # if self.surfaces_chlo[0]==0.:
                    # self.surfaces_chlo = self.surfaces_chlo[1:]
                # if age_switch==1:
                    # self.age_physio_edge = 0.
                    # self.change_status_edge()
                # else:
                    # self.age_physio_edge = age_switch
        
            # # Complete the evolution of the lesion up to the end of time step
            # if self.is_active:
                # self.age_dday += self.dt_left_after_senescence*self.ddday/self.dt_before_senescence
                # self.update_status()

    def update_status(self):
        """ Update growth demand and status. """        
        # A REORDONNER
        f = self.fungus
        status = self.status
        if status <= f.INCUBATING:
            self.incubation()
        elif status <= f.CHLOROTIC:
            self.chlorosis()
        elif status <= f.NECROTIC:
            if self.status_edge <= f.CHLOROTIC:
                self.chlorosis()
            self.necrosis()
        else:
            if self.status_edge <= f.CHLOROTIC:
                self.chlorosis()
            if self.status_edge <= f.NECROTIC:
                self.necrosis()
            self.sporulation()
        
    def change_status(self):
        """ Passes the status of the center of the lesion from one status to the next. """
        self.status += 1
    
    def change_status_edge(self):
        """ Passes the status of the edge of the lesion from one status to the next. """
        self.status_edge +=1
    
    def reset_age_physio(self):
        """ Turn age physio to 0. """
        self.age_physio = 0.
    
    def reset_age_physio_edge(self):
        """ Turn age physio to 0. """
        self.age_physio_edge = 0.
    
    # def can_compute_rings(self):
        # """ Allow ring formation only if particular first ring has reached Smin. """
        # if self.status_edge<=self.fungus.CHLOROTIC:
            # return (len(self.surfaces_chlo)>0. and round(self.surface_first_ring, 14) == round(self._surface_min, 14))
        # elif self.status_edge<=self.fungus.NECROTIC:
            # return (len(self.surfaces_nec)>0. and round(self.surface_first_ring, 14) == round(self._surface_min, 14))
            
    def compute_all_surfaces(self):
        """ Not needed in this model.
        ..TODO: Remove """
        pass
    
    def set_position(self, position=None):
        """ Set the position of the DU to position given in argument.
        
        Parameters
        ----------
        position: list[x, y] on leaf blade
            Position of the DU.
        """
        if not is_iterable(position[0]):
            self.position = [position]
        else:
            self.position = position
        
    def set_nb_spores(self, nb_spores=0.):
        """ Set the position of the DU to position given in argument.
        
        Parameters
        ----------
        nb_spores: int
            Number of spores in the DU forming the lesion
        """
        self.nb_spores = nb_spores
    
    @property
    def surface_inc(self):
        """ Calculate the surface in incubation. """
        return self.surface_first_ring if self.is_incubating() else 0.

    @property
    def surface_chlo(self):
        """ Calculate the surface in chlorosis. """
        return (sum(self.surfaces_chlo)+self.surface_first_ring) if self.is_chlorotic() else sum(self.surfaces_chlo)

    @property
    def surface_nec(self):
        """ Calculate the surface in necrosis. """
        return (sum(self.surfaces_nec)+self.surface_first_ring) if self.is_necrotic() else sum(self.surfaces_nec)

    @property
    def surface_spo(self):
        """ Calculate the surface in sporulation. """
        return sum(self.surfaces_spo)
    
    @property
    def necrotic_area(self):
        """ calculate surface necrotic + sporulating. """
        return self.surface_nec + self.surface_spo + self.surface_empty
    
    @property
    def surface_alive(self):
        """ Calculate the surface alive of the lesion. """
        return self.surface_inc + self.surface_chlo + self.surface_nec + self.surface_spo + self.surface_empty

    @property
    def _surface_min(self):
        """ Calculate the surface min for a lesion in chlorosis. """
        return self.fungus.Smin * self.nb_lesions
        
    @property
    def _surface_max(self):
        """ Calculate the surface max for a lesion. """
        # return self.fungus.Smax * self.nb_lesions + self.surface_dead
        return self.fungus.Smax * self.nb_lesions + self.surface_senesced
    
    @property
    def surface(self):
        """ Calculate the total surface of the lesion. """
        return self.surface_alive + self.surface_dead
        
    @property
    def nb_lesions(self):
        if self.position is None:
            return None
        else:
            return len(self.position)

class SeptoriaFungus(Fungus):
    def __init__(self, name='septoria', Lesion=SeptoriaAgePhysio, DispersalUnit=SeptoriaDU, parameters=septoria_parameters):
        super(SeptoriaFungus, self).__init__(name=name, Lesion=Lesion, DispersalUnit=DispersalUnit, parameters=parameters)
        
# class Parameters(_SeptoriaParameters):
    # def __init__(self,**kwds):
        # _SeptoriaParameters.__init__(self, model=SeptoriaAgePhysio, **kwds)
        
# class Disease(_Disease):
    # @classmethod
    # def parameters(cls, **kwds):
        # return Parameters(**kwds)
    
    # @classmethod
    # def lesion(cls, **kwds):
        # SeptoriaAgePhysio.fungus=cls.parameters(**kwds)
        # return SeptoriaAgePhysio

import collections        
def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)