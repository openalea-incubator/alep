""" Class of lesion of wheat septoria. 
    Progress from one state to another is computed with physiological age.
"""

# Imports #########################################################################
from alinea.alep.fungal_objects import *
from alinea.alep.septoria import SeptoriaDU, septoria_parameters, is_iterable
import numpy as np
from scipy.interpolate import interp1d
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
        mutable: bool
            False if all instances of the class share the same parameters
        """
        super(SeptoriaAgePhysio, self).__init__(mutable=mutable)
        # Calculate parameter for emission from other parameters
        self.density_dus_emitted_max = np.array([self.fungus.density_dus_emitted_ref*(self.fungus.reduction_by_rain**ind) for ind in range(self.fungus.rain_events_to_empty)])
        # Status of the center of the lesion
        self.status = self.fungus.INCUBATING
        # Age of the center of the lesion
        self.ddday=None
        self.age_tt = 0.
        self.age_physio = 0.
        # Status of the periphery of the lesion
        self.status_edge = self.fungus.INCUBATING
        # Age of the periphery of the lesion
        self.age_physio_edge = 0.
        # Ratio left to progress in new status when centre of the lesion changes status during time step
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
        self.surface_dead = 0.
        # Marker of incubation completion
        self.incubation_completed = False
        # Marker of senescence completion
        self.senescence_response_completed = False
        # Number of lesions in senescence
        self.nb_lesions_sen = 0.
        # Potential surface of the lesion if no competition
        self.potential_surface = 0.

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
    
    def update(self, dt=1., leaf=None):
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
        if any([x[0]<=leaf.senesced_length for x in self.position]):
            self.senescence_response(leaf.senesced_length)

        if self.is_active:
            # Compute delta degree days in dt
            self.compute_delta_ddays(dt, leaf)
            
            if self.ddday > 0.:
                # Update age in degree days of the lesion
                self.age_tt += self.ddday        
                # Update growth demand and status
                self.update_status()
                # Disable growth if edge is necrotic
                self.check_edge()
            else:
                self.growth_demand = 0.
                
            # Update potential surface
            self.potential_surface += self.growth_demand 
        
        if (self.incubation_completed == False and
            round(self.surface_first_ring,14) >= round(self._surface_min, 14)):
            self.incubation_completed = True

    def check_edge(self):
        if (self.status > self.fungus.CHLOROTIC and 
            round(self.surface_chlo, 14) == 0. and 
            round(self.surface_inc, 14) == 0.):
            self.disable_growth()
            if (self.status >= self.fungus.SPORULATING and 
                round(self.surface_nec, 14) == 0. and
                round(self.surface_spo, 14) == 0.):
                    self.disable()
                
    def rh_response(self, rh):
        f = self.fungus
        if rh > f.rh_max:
            return 1
        elif rh > f.rh_min:
            rh_resp = interp1d([f.rh_min, f.rh_max], [0,1])
            return rh_resp(rh).tolist()
        else:
            return 0.
                
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
        props = leaf.properties()
        # Calculation
        if dt != 0.:
            ddday = sum([max(0,(temp - f.basis_for_dday)*1/24.) if temp<=f.temp_max else 0. for temp in props['temperature_sequence']])
            # TEMP ... et moche
            if f.rh_effect==True:
                if f.apply_rh=='all' :
                    rh_resps = map(self.rh_response, props['relative_humidity_sequence'])
                    rh_factor = np.mean(rh_resps)
                    ddday *= rh_factor
                elif f.apply_rh=='chlorosis' and self.status>=f.CHLOROTIC:
                    rh_resps = map(self.rh_response, props['relative_humidity_sequence'])
                    rh_factor = np.mean(rh_resps)
                    ddday *= rh_factor
                elif f.apply_rh=='necrosis' and self.status>=f.NECROTIC:
                    rh_resps = map(self.rh_response, props['relative_humidity_sequence'])
                    rh_factor = np.mean(rh_resps)
                    ddday *= rh_factor
#            if self.is_necrotic() and f.rh_effect==True:
#                rh_resps = map(self.rh_response, props['relative_humidity_sequence'])
#                rh_factor = np.mean(rh_resps)
#                ddday *= rh_factor
#                if leaf.dd_sequence[-1]>880 and leaf.dd_sequence[-1]<1500:
#                    import pdb
#                    pdb.set_trace()
#                    print '-----------------------'
#                    print 'rh_factor %f' %rh_factor
#                    print '-----------------------'
#            if 'global_efficacy' in leaf.properties():
#                ddday *= (1 - max(0, min(1, leaf.global_efficacy['eradicant'])))
        else:
            ddday = 0.
        
        if ddday > f.degree_days_to_chlorosis or ddday > f.degree_days_to_necrosis or ddday > f.degree_days_to_sporulation:
            raise SeptoError('Can not handle a dt > minimum stage duration')
        
        # Save variable
        self.ddday = ddday
    
    def progress(self, age_threshold=0.):
        """ Compute progress in physiological age according to age_threshold. 
        """
        left = self.ratio_left
        if left==0.:
            # Normal situation
            progress = self.ddday/age_threshold if age_threshold>0. else 0.
        elif left>0.:
            # Passing from one state to another in same time step
            progress = left*self.ddday/age_threshold if age_threshold>0. else 0.
            # Reset ratio left and age physio
            self.ratio_left = 0.
        return progress
    
    def control_growth(self, growth_offer=0.):
        """ Reduce surface of the lesion up to available surface on leaf.
        
        Parameters
        ----------
        growth_offer: float
            Surface available on the leaf for the ring to grow (cm2)
        """
        if self.growth_is_active:
            # Growth offer is added to surface according to state
            if self.is_incubating():
                if growth_offer<0:
                    self.surface_dead -= growth_offer                   
                self.surface_first_ring = max(0, self.surface_first_ring + growth_offer)
                if self.surface_first_ring == 0. and self.age_tt>self.ddday:
                    self.disable()
                    return
            else:
#                growth_offer = min(self.surface_inc)
                Smin = self._surface_min
                if self.is_chlorotic() and round(self.surface_first_ring,14) < round(Smin,14):
                    if self.age_tt - self.ddday < self.fungus.degree_days_to_chlorosis:
                        if growth_offer <= 0:
                            self.surface_dead -= growth_offer                   
                            self.surface_first_ring = max(0, self.surface_first_ring + growth_offer)
                            if self.surface_first_ring == 0. and self.age_tt>self.ddday:
                                self.disable()
                                return
                        else:
                            diff = self.fungus.degree_days_to_chlorosis - (self.age_tt - self.ddday)
                            ratio_inc = diff/self.ddday
                            self.surface_first_ring += growth_offer*ratio_inc
                            growth_offer *= (1 - ratio_inc)
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
                    self.surfaces_chlo += surf[:len(self.surfaces_chlo)]
                    self.surfaces_chlo = np.append(self.surfaces_chlo, surf[len(self.surfaces_chlo):])
            
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
            if self.growth_is_active:
                self.growth_demand = progress * f.Smin * self.nb_lesions_non_sen
        else:
            if self.growth_is_active:
                ratio_left = 1 - (self.age_physio-progress)
                self.growth_demand = progress * f.Smin * self.nb_lesions_non_sen *ratio_left

            # Change status
            self.ratio_left = round((self.age_physio - 1.)/progress, 14)
            self.change_status()
            self.change_status_edge()
            self.reset_age_physio()
            self.incubation_completed = True
            self.chlorosis()

    def chlorosis(self):
        """ Compute growth demand and physiological age progress to necrosis.
        """
        f = self.fungus
        # Compute progress in chlorosis
        time_to_nec = f.degree_days_to_necrosis
        progress = self.progress(age_threshold=time_to_nec)

        # Compute growth demand
        if self.growth_is_active:
            # Note : '+=' because might be added to growth demand in incubation 
            # if transition in same time step; added to zero otherwise.
            self.growth_demand += progress * time_to_nec * f.growth_rate * self.nb_lesions_non_sen
            # Limit growth to size max
            Smax = self._surface_max
            if self.surface + self.growth_demand >= Smax:
                self.growth_demand = max(0., Smax - self.surface)

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
                new_ends = np.arange(max(width, width*ceil(begs_prog[0]/width)), width*(ceil(ends_prog[-1]/width)+1), width)
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
                                              self.surfaces_nec[ind]/(ends[ind]-begs[ind]), 14) if round(ends[ind]-begs[ind], 14) > 0. else 0.
                                    
                # Get what passes to next status
                self.surfaces_nec = np.extract(new_ends<=1, new_surf)
                surface_to_next_phase = sum(np.extract(new_ends>1, new_surf))
                self.to_sporulation = f.sporulating_capacity * surface_to_next_phase
                self.surface_empty += (1 - f.sporulating_capacity) * surface_to_next_phase
                
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
        if sum(self.surfaces_spo) == 0:
            self.surfaces_spo[0] += self.surface_first_ring
            self.surface_first_ring = 0.
        self.surfaces_spo[0] += self.to_sporulation
        self.to_sporulation = 0.

    def emission(self, density_DU_emitted):
        """ Return number of DUs emitted by the lesion. """
        if density_DU_emitted>0:       
            f = self.fungus
            density = map(lambda x: min(float(density_DU_emitted),x), self.density_dus_emitted_max)
            du_factor = map(lambda (x,y): x / y if y>0. else 0., zip(density, self.density_dus_emitted_max))
            emissions = map(lambda x: int(x), self.surfaces_spo * density)
            
            delta_spo = self.surfaces_spo*du_factor
            self.surfaces_spo -= delta_spo
            self.surfaces_spo[1:] += delta_spo[:-1]
                       
            if delta_spo[-1]>0.:
                self.surface_empty += delta_spo[-1]
                
            if (self.status_edge == f.SPORULATING and 
                round(self.surface_spo,14)<=round(f.threshold_spo*self.nb_lesions,14)):
                # Everything becomes empty and the lesion is disabled
                self.surface_empty += sum(self.surfaces_spo)
                self.surfaces_spo = np.zeros(len(self.surfaces_spo))
                self.change_status()
                self.disable()
            return sum(emissions)
        else:
            return 0.
    
    def senescence_response(self, senesced_length=0.):
        """ Compute surface alive and surface dead after senescence. """
        if not self.is_senescent:
            self.become_senescent()
            
        if not self.senescence_response_completed:
            nb_sen = len(filter(lambda x: x[0]<=senesced_length, self.position))
            nb_new_sen = nb_sen - self.nb_lesions_sen
            ratio_sen = float(nb_new_sen)/(self.nb_lesions_non_sen)
            
            # Reduce surfaces alive
            f = self.fungus
            age_switch = f.age_physio_switch_senescence
            age_physio = self.age_physio
            age_edge = self.age_physio_edge
            if self.status < f.CHLOROTIC:
                self.surface_dead += self.surface_first_ring * ratio_sen
                self.surface_first_ring = self.surface_first_ring * (1-ratio_sen)
            elif self.is_chlorotic() and age_switch >= age_physio:
                self.surface_dead += self.surface_first_ring * ratio_sen + sum(self.surfaces_chlo * ratio_sen)
                self.surface_first_ring = self.surface_first_ring * (1-ratio_sen)
                self.surfaces_chlo *= (1-ratio_sen)
                self.surfaces_chlo = self.surfaces_chlo[self.surfaces_chlo>0.]
            elif self.status_edge == f.CHLOROTIC and age_switch > age_edge and self.surface_chlo > 0.:
                default_nb_rings = f.nb_rings_by_state
                (rings, width) = np.linspace(0,1, default_nb_rings+1, retstep=True)
                if self.is_chlorotic():
                    rings = rings[:ceil(age_physio/width)+1]
                    rings[-1] = age_physio
                rings = rings[floor(age_edge/width):]
                rings[0] = age_edge
                if age_switch>rings[0] and len(self.surfaces_chlo)>0:
                    ind_cut = np.where(rings<age_switch)[0][-1]
                    ratio_to_dead = round((age_switch - rings[ind_cut])/(rings[ind_cut+1]-rings[ind_cut]), 14)
                    to_dead_on_cut_ring = self.surfaces_chlo[ind_cut]*ratio_sen*ratio_to_dead
                    self.surface_dead += to_dead_on_cut_ring # In the ring which is cut
                    self.surface_dead += sum(self.surfaces_chlo[:ind_cut]*ratio_sen) # In older rings
                    self.surfaces_chlo[:ind_cut] *= (1 - ratio_sen)
                    self.surfaces_chlo[ind_cut] -= to_dead_on_cut_ring
                    
                if nb_new_sen==self.nb_lesions_non_sen:
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

            # Update number of lesions in senescence
            self.nb_lesions_sen += nb_new_sen         

            # Update potential surface
            self.potential_surface *= (1-ratio_sen)

            # Stop growth when last lesion of cohort is reached
            if self.nb_lesions_non_sen==0.:
                if self.is_incubating() or (self.is_chlorotic() and age_switch >= age_physio):
                    self.senescence_response_completed=True
                    self.disable()
                else:
                    self.disable_growth()
                    self.senescence_response_completed=True

    def update_status(self):
        """ Update growth demand and status. """        
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
            
    def disappear(self):
        """ Kill the lesion and pass all surfaces to 0. 
        """
        self.growth_is_active = False
        self.growth_demand = 0.
        self.surface_inc = 0.
        self.surface_chlo = 0.
        self.surface_nec = 0.
        self.surface_spo = 0.
        self.surface_empty = 0.
        self.surface_dead = 0.
        self.disable() 

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
        return self.fungus.Smin * self.nb_lesions_non_sen
        
    @property
    def _surface_max(self):
        """ Calculate the surface max for a lesion. """
        return self.fungus.Smax * self.nb_lesions_non_sen + self.surface_senescent
    
    @property
    def surface(self):
        """ Calculate the total surface of the lesion. """
        return self.surface_alive + self.surface_dead

    @property
    def nb_lesions(self):
        """ Get number of lesions in cohort if group_dus == True.
            Each lesion in cohort has its own position. """
        if self.position is None:
            return None
        elif self.fungus.group_dus == True:
            return len(self.position)
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
    def surface_non_senescent(self):
        """ calculate the surface of the lesion non affected by senescence. """
        if self.is_senescent:
            f = self.fungus
            age_switch = f.age_physio_switch_senescence
            ratio_non_sen = self.nb_lesions_non_sen / self.nb_lesions
            return self.surface_inc + ratio_non_sen*(age_switch*self.surface_chlo +\
                    self.surface_nec + self.surface_spo + self.surface_empty)
        else:
            return self.surface - self.surface_dead
            
    @property
    def surface_senescent(self):
        """ Calculate the surface of the lesion affected by senescence. """
        return self.surface - self.surface_non_senescent

            
class SeptoriaFungus(Fungus):
    def __init__(self, Lesion=SeptoriaAgePhysio, DispersalUnit=SeptoriaDU, parameters=septoria_parameters):
        super(SeptoriaFungus, self).__init__(Lesion=Lesion, DispersalUnit=DispersalUnit, parameters=parameters)

# Useful functions ############################################################
#import collections
#def is_iterable(obj):
#    """ Test if object is iterable """
#    return isinstance(obj, collections.Iterable)
    
class SeptoError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


# Temp
def rh_response(x, rh_max=60, rh_min=35):
    if x > rh_max:
        return 1
    elif x > rh_min:
        f = interp1d([rh_min,rh_max], [0,1])
        return f(x).tolist()
    else:
        return 0.
