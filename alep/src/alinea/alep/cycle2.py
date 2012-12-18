from random  import random
from math import exp
from math import pi
import numpy as np
from pylab import *

nbFract = 10 # Number of fractions into a state

class LesionFactory(object):
    """
    """

    def __init__(self, fungus):
        """ Initialize the lesion. 
        
        :Parameters:
          - `fungus` (function): returns a class of specific parameters for 
          the chosen fungus (e.g. 'septoria()' or 'powderyMildew()').

        """
        self.fungus = fungus    
    
    def instantiate(self):
        """ instantiate a Lesion """
        les = Lesion(self.fungus, nbFract)
        les.initiate()
        return les
        
class Lesion(object):
    """ 
    """
    
    def __init__(self, fungus, nbFract=nbFract):
        """ Initialize the lesion. 
        
        :Parameters:
          - `fungus` (function): returns a class of specific parameters for 
          the chosen fungus (e.g. 'septoria()' or 'powderyMildew()').

        """
        self.surfaces = []
        self.fungus = fungus
        self.nbFract = nbFract
    
     def initiate(self):
        """ Create a new property 'lesions' on g[vid] 
        
        """
        new_surface = self.fungus_factory()
        new_surface.status = EMERGENT
        self.surfaces.append(new_surface)
    
    def growth(iSim, iStep, t, s, dt, Dt, growth_rate):
        """ 
        Parameters
        ==========

        - iSim : int
            current step of simulation
        - iStep : int
            current step of simulation into the state
        - t : float
            time needed to achieve state
        - s : array
            state variable
        - dt : float
            step of the simulation
        - Dt : float
            step for grouping surfaces
        - growth_rate : float
            constant input = growth by time step

        Returns
        =======
        - output_growth_rate: float
            produced growth rate in the new state
        """
        ds = s * (float(dt) / float(Dt))
        
        # If no surface has completed this state
        if iSim < t:
            if iStep == 1:
                # current_Dt is the fraction of the state in which the simulation is currently running
                current_Dt = 0
            else:
                current_Dt = floor(float(iStep)/Dt)

            # Balance between input and output until current_Dt
            s[0:current_Dt+1] -= ds[0:current_Dt+1]
            s[1:current_Dt+1] += ds[0:current_Dt]    
        
        # If this state has been completed by at least 1 surface, same as is 'test_temporal.py'
        else:
            ds = s * (float(dt) / float(Dt))
            s = s - ds
            s[1:] += ds[0:-1]
        
        # In all cases the input is the same
        s[0]+= growth_rate
        
        return s, ds[-1]
    
    def fungus_factory(self):
        """ This factory will search for entry points.
        
        Thus, we will not have to change this file to add a new fungus.
        """
        eval(self.fungus.name)()
            return eval(self.fungus.name)()
        # TODO G.Garin : 18/12/12 : Pas bien d'utiliser eval non ? 

        
##############################################################################

class Septoria(Lesion):
    """ Lesion of Septoria.
    """
    EMERGENT = 0
    INCUBATING = 1
    CHLOROTIC = 2
    NECROTIC = 3
    SPORULATING = 4
    EMPTY = 5
    DEAD = 6

    def __init__(self, status = EMERGENT):
        """ Initialize the lesion. 

        """
        super(Septoria, self).__init__()
        for i in range(sef.nbStates):
            self.surfaces.append(np.zeros(self.nbFract))
        self.growth_rate = np.zeros(self.nbStates)
        self.age_dday = 0.
    
    def update(self, dt, leaf, **kwds):
        """ Update the status of the lesion and make it grow if needed.

        """        
        
        # Update the age of the lesion
        ddday = (leaf.temp - self.basis_for_dday)/dt
        self.age_dday += ddday

        # Update the surfaces in each state
        for state in range(CHLOROTIC, SPORULATING):
            iStep[state] = dt - tEnd[state-1]
            self.surfaces[state], growth_rate[state+1] = growth(self.age_dday, iStep, 
                                                                tEnd[state], 
                                                                self.surfaces[state],
                                                                dt, self.nbFrac, growth_rate[state])


    @property
    def surface(self):
        """ Compute the surface of the lesion.
        """
        surf = sum(self.surfaces)
        return surf
    
    @property
    def status(self):
        """ Compute the status of the lesion.
        """
        if self.rings:
            return self.rings[0].status
    
        
        
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
                 nbStates = 7
                 basis_for_dday = -2,
                 temp_min = 10,
                 temp_max = 30,
                 wd_min = 10,
                 loss_rate = 1./120,
                 degree_days_to_chlorosis = 220,
                 degree_days_to_necrosis = 110,
                 degree_days_to_sporulation = 50,
                 epsilon = 0.001,
                 Smin = 0.03,
                 Smax = 0.3,
                 growth_rate = 0.0006,
                 rh_min = 85,
                 rain_events_to_empty = 3,            
                 *args, **kwds):
        """ Parameters for septoria.
        
        :Parameters:
            - 'nbStates' (int): number of successive states on a lesion
            - `basis_for_dday` (float): basis temperature for the accumulation of degree days (°C)
            - `temp_min` (float): Minimal temperature for infection.
            - `temp_max` (float): Maximal temperature for infection.
            - `wd_min` (float): Minimal wetness duration for infection.
            - `loss_rate` (float): Loss rate of dispersal units in 1 hour.
            - `degree_days_to_chlorosis` (float): Thermal time between emergence 
            and chlorosis (i.e. incubation for the first rings).
            - `degree_days_to_necrosis` (float): Thermal time between chlorosis 
            and necrosis (i.e. incubation for the first rings).
            - `degree_days_to_sporulation` (float): Thermal time between necrosis 
            and sporulation.
            - `epsilon` (float): Initial size of incubating lesion (cm2).
            - `Smin` (float): Initial size of chlorotic lesion (cm2).
            - `Smax` (float): Lesion maximum size (cm2).
            - `growth_rate` (float): Lesion growth rate (cm2.dday-1).
            - `rh_min` (float): Minimal relative humidity for sporulation.
            - `rain_events_to_empty` (int): number of rain events to empty a sporulating ring.
            
        """
        self.name = "Septoria"
        self.nbStates = nbStates
        self.basis_for_dday = basis_for_dday
        self.temp_min = temp_min
        self.temp_max = temp_max
        self.wd_min = wd_min
        self.loss_rate = loss_rate
        self.degree_days_to_chlorosis = degree_days_to_chlorosis
        self.degree_days_to_necrosis = degree_days_to_necrosis
        self.degree_days_to_sporulation = degree_days_to_sporulation
        self.epsilon = epsilon
        self.Smin = Smin
        self.Smax = Smax
        self.growth_rate = growth_rate
        self.rh_min = rh_min
        self.rain_events_to_empty = rain_events_to_empty

def septoria(**kwds):
    return SeptoriaParameters(**kwds)
    
##############################################################################
def proba(p):
    """ p in 0,1 """
    return random() < p