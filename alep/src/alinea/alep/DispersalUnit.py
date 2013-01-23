""" TODO : Description Dispersal_unit """

from random  import random

class DispersalUnit(object):
    """    
    Define a dispersal unit protocol.

    To implement a lesion, you need to implement the following methods:
        - infect(dt, leaf, ...
    
    """
    
    def __init__(self, fungus, nbSpores):
        """ Initialize the dispersal unit. 
        """
        self.fungus = fungus
        self.nbSpores = nbSpores
        
    def infect(self, dt, leaf, **kwds):
    # Si infection réussie :
    # - la dispersal unit disparait
    # - une lésion est instanciée
    
        pass

class SeptoriaDU(DispersalUnit):

    def __init__(self, fungus):
        """ Initialize the dispersal unit of septoria. 
              
        """
        super(SeptoriaDU, self).__init__(fungus=fungus, nbSpores=nbSpores)
        self.cumul_wetness = 0.
    
    def infect(self, dt, leaf, **kwds):
        """ Check if infection is a success.
        
        Return boolean (True) or (False).
        
        """
        
        leaf_wet = leaf.wetness # (boolean): True if the leaf sector is wet during this time step.
        temp = leaf.temp # (float) : mean temperature on the leaf sector during the time step (in °C).
        healthy_surface = leaf.healthy_surface # (float) : healthy surface (=with no lesion) on the leaf sector during the time step (in cm^2).
        
        # TODO: design a new equation : see Magarey (2005)
        if leaf_wet:
            self.cumul_wetness += 1
        elif self.cumul_wetness > 0: 
            assert not leaf_wet
            self.cumul_wetness = 0
            self.status = self.DEAD 
        else:
            assert not leaf_wet
            assert self.cumul_wetness == 0
        
        if (self.fungus.temp_min <= temp <= self.fungus.temp_max) and self.cumul_wetness >= self.fungus.wd_min :
            # TODO : create a function of the number of spores            
            spores_factor = nbSpores / nbSpores # always equals 1 for now
            if proba(spores_factor):
                return True
            else:
                return False
        elif self.cumul_wetness == 0 :
            return False
            # TODO : Proba conditionnelle doit se cumuler.
            if proba(self.fungus.loss_rate): 
                self.nbSpores = 0.
        

class PowderyMildewDU(DispersalUnit):

    def __init__(self, fungus):
        """ Initialize the dispersal unit of septoria. 
        """
        super(PowderyMildewDU, self).__init__(fungus=fungus, nbSpores=nbSpores)
    
    def infect(self, dt, leaf, **kwds):
        """ Check if infection is a success.
        
        Return boolean (True) or (False).
        
        """
        pass
        
####################################################################################################
# Useful functions
def proba(p):
    """ p in 0,1 """
    return random() < p

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
                 temp_min = 10,
                 temp_max = 30,
                 wd_min = 10,
                 loss_rate = 1./120,
           
                 *args, **kwds):
        """ Parameters for septoria.
        
        :Parameters:
            - `temp_min` (float): Minimal temperature for infection.
            - `temp_max` (float): Maximal temperature for infection.
            - `wd_min` (float): Minimal wetness duration for infection.
            - `loss_rate` (float): Loss rate of dispersal units in 1 hour.
            
        """
        self.name = "Septoria"
        self.Dt = Dt
        self.basis_for_dday = basis_for_dday
        self.temp_min = temp_min
        self.temp_max = temp_max
        self.wd_min = wd_min
        self.loss_rate = loss_rate

    def __call__(self):
        return Septoria(fungus=self)

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
        
    def __call__(self):
        return PowderyMildew(fungus=self)

def powdery_mildew(**kwds):
    return PowderyMildewParameters(**kwds)