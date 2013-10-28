""" Classes of dispersal unit, lesion and ring specific of wheat septoria.

"""
# Imports #########################################################################
from alinea.alep.fungal_objects import *
# from alinea.alep import septoria_continuous, septoria_with_rings, septoria_exchanging_rings
#from alinea.alep.septoria_continuous import *
#from alinea.alep.septoria_with_rings import *
#from alinea.alep.septoria_exchanging_rings import *
from openalea.vpltk import plugin
from random import random, randint, seed
from math import floor, ceil
import numpy as np
seed(1)

# Dispersal unit ##################################################################
class SeptoriaDU(DispersalUnit):
    """ Define a dispersal unit specific of septoria.
    
    """
    fungus = None
    def __init__(self, nb_spores=None, position=None, status=None):
        """ Initialize the dispersal unit of septoria.
        
        Parameters
        ----------
        nb_spores: int
            Number of spores aggregated in the dispersal unit
        position: non defined
            Position of the dispersal unit on the phyto-element
        status: str
            'emitted' or 'deposited'
        
        Returns
        -------
            None
        """
        super(SeptoriaDU, self).__init__(nb_spores=nb_spores, position=position, status=status)
        self.cumul_wetness = 0.
            
    def infect(self, dt, leaf, **kwds):
        """ Compute infection by the dispersal unit of Septoria.
        
        Parameters
        ----------
        dt: int
            Time step of the simulation (in hours)
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. area, green area, healthy area,
            senescence, rain intensity, wetness, temperature, lesions)

        Returns
        -------
            None
        """
        leaf_wet = leaf.wetness # (boolean): True if the leaf sector is wet during this time step.
        temp = leaf.temp # (float) : mean temperature on the leaf sector during the time step (in degree).
        if self.can_infect_at_position:
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
                            # print('l.91 infection septo')
                            # raise Exception('')
                    elif self.cumul_wetness == 0 :
                        # TODO : Proba conditionnelle doit se cumuler.
                        if proba(self.fungus.loss_rate): 
                            self.disable()
        else:
            self.disable()

# Fungus parameters (e.g. .ini): config of the fungus #############################
class SeptoriaParameters(Parameters):
    model = None
    def __init__(self,
                 model=None,
                 INCUBATING = 0,
                 CHLOROTIC = 1,
                 NECROTIC = 2,
                 SPORULATING = 3,
                 EMPTY = 4,
                 DEAD = 5,
                 delta_age_ring = 20.,
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
                 treshold_spores = 1000,
                 *args, **kwds):
        """ Parameters for septoria.
        
        Parameters
        ----------
        model: model of lesion
            Model of lesion among : "ContinuousSeptoria", "SeptoriaExchangingRings", 
            "SeptoriaWithRings"
        delta_age_ring: int
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
        production_rate: float
            Number of spores produced by cm2 of new sporulating surface
        treshold_spores: int
            Number of spores left in stock to consider the lesion empty
        """
        self.name = "septoria"
        self.__class__.model = model
        self.INCUBATING = INCUBATING
        self.CHLOROTIC = CHLOROTIC
        self.NECROTIC = NECROTIC
        self.SPORULATING = SPORULATING
        self.EMPTY = EMPTY
        self.DEAD = DEAD
        self.delta_age_ring = delta_age_ring
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
        # TODO : Improve this parameter. Very Sensitive.
        self.treshold_spores = treshold_spores
        
    def __call__(self, nb_spores=None, position=None):
        model = self.model
        
        if model.fungus is None:
            model.fungus = self
        if SeptoriaDU.fungus is None:
            SeptoriaDU.fungus = self
        return model(nb_spores=nb_spores, position=position)

def septoria(**kwds):
    return SeptoriaParameters(**kwds)

class Disease(object):
    name = 'septoria'

    @classmethod
    def parameters(cls, **kwds):
        return septoria(**kwds)
    
    @classmethod
    def dispersal_unit(cls, **kwds):
        SeptoriaDU.fungus=cls.parameters(**kwds)
        return SeptoriaDU
           
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

# Plugin function #################################################################
def plugin_septoria(model='septoria_exchanging_rings'):
    diseases=plugin.discover('alep.disease')
    try:
        septoria = diseases[model].load()
    except KeyError:
        if model=='septoria_exchanging_rings':
            from alinea.alep.septoria_exchanging_rings import Disease
        elif model=='septoria_continuous':
            from alinea.alep.septoria_exchanging_rings import Disease
        elif model=='septoria_with_rings':
            from alinea.alep.septoria_exchanging_rings import Disease
        septoria=Disease()
    return septoria        