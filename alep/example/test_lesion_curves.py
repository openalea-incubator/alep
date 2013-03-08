""" Draft created to represent the lesion as a curve of age of mycelium
    as a function of its distance from the center of the lesion.
    
"""
from random import random
from random import randint
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

def test_lesion():
    dt = 1
    nb_steps = 400
    leaf = Leaf()
    sl = SeptoriaLesion(fungus = septoria(), nbSpores = 0., position=None)
    for i in range(nb_steps):
        print(i)
        sl.update(dt, leaf)
        
    return sl
    
def test_draw_lesion():
    dt = 1
    nb_steps = 500
    leaf = Leaf()
    sl = SeptoriaLesion(fungus = septoria(), nbSpores = 0., position=None)

    x_values = np.arange(0., 0.3, 0.01)
    y_values = np.arange(0., nb_steps)
    X, Y = np.meshgrid(x_values, y_values)
    
    Z = np.zeros((len(y_values), len(x_values)))
    
    for i in range(len(Y)):
        print(i)
        sl.update(dt, leaf)
        Z[i] = sl.compute_growth_curve(X[i])
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
    plt.show()
    
    return Z, sl

class Leaf(object):
    def __init__(self):
        self.temp = 24.
        
class Lesion(object):
    """ Define the interface of a lesion.
    
    """
    def __init__(self, fungus=None, nbSpores=0., position=None):
        """ Initialize the lesion.
        
        """
        self.fungus = fungus
        self.active = True
        self.nbSpores = nbSpores
        self.position = position
        self.age_physio = 0.
        self.stock_du = []
        self.emissions = []

class SeptoriaLesion(Lesion):
    """ Define the interface of a lesion of septoria.
    
    """
    def __init__(self, fungus, nbSpores, position):
        """ Initialize the lesion of septoria.
        
        """
        super(SeptoriaLesion, self).__init__(fungus=fungus, nbSpores=nbSpores, position=position)
        self.age_dday = 0.
        self.surface_total = 0.
        self.surface_incubating = 0.
        self.surface_chlorotic = 0.
        self.surface_necrotic = 0.
        self.surface_sporulating = 0.
        self.surface_dead = 0.
        self.is_growing = True
        self.count= []
    
    def compute_ddday(self, dt, leaf) :
        """ Compute the number of degree days cumulated in dt.
        
        """
        ddday = (leaf.temp - self.fungus.basis_for_dday)/(24./dt)
        return ddday
    
    def compute_growth_rate(self):
        """ Find growth rate on a piecewise linear curve.
                
        """
        if self.is_growing:
            if self.age_dday < self.fungus.degree_days_to_chlorosis : 
                r = self.fungus.r1
            else:
                r = self.fungus.r2
        else:
            r = self.fungus.r3
            
        return r 
        
    def update(self, dt, leaf):
        """ Update the lesion of septoria (ageing and growth)
        
        """
        # Ageing :
        ddday = self.compute_ddday(dt, leaf)
        self.age_dday += ddday

        # Growth :  
        r = self.compute_growth_rate()
        self.surface_total = min(self.fungus.Smax, self.surface_total + r*ddday)
        self.count.append(self.surface_total)
        self.update_surfaces_in_state()

    def update_surfaces_in_state(self):
        """ Compute the surfaces in each state of the lesion except dead:
            - incubating
            - chlorotic
            - necrotic
            - sporulating
            
        """
        if self.surface_total < self.fungus.Smin:
            # print(self.age_dday)
            # print('incubating')
            self.surface_incubating = self.surface_total
            self.surface_chlorotic = 0.
            self.surface_necrotic = 0.
            self.surface_sporulating = 0.
        else:
            self.surface_incubating = 0.
                        
            # solve first system of linear equations
            surface_max_sporulating = self.surface_total + (self.fungus.degree_days_to_sporulation)*(self.fungus.Smin-self.surface_total)/self.age_dday
                        
            # solve second system of linear equations
            surface_max_necrotic = self.surface_total + (self.fungus.degree_days_to_necrosis)*(self.fungus.Smin-self.surface_total)/self.age_dday

            
            if surface_max_necrotic < self.fungus.Smin :
                print(self.age_dday)
                print('chlorotic')
                self.surface_chlorotic = self.surface_total
                self.surface_necrotic = 0.
                self.surface_sporulating = 0.
            elif surface_max_sporulating < self.fungus.Smin:
                assert(surface_max_necrotic >= self.fungus.Smin)
                print(self.age_dday)
                print('necrotic')
                self.surface_chlorotic = self.surface_total - surface_max_necrotic
                self.surface_necrotic = surface_max_necrotic
                self.surface_sporulating = 0.
            else:
                assert(surface_max_sporulating >= self.fungus.Smin)
                self.surface_chlorotic = self.surface_total - surface_max_necrotic
                self.surface_necrotic = self.surface_chlorotic - surface_max_sporulating
                self.surface_sporulating = surface_max_sporulating
                
    def compute_growth_curve(self, x_values):
        """ Compute growth curve for the lesion in given age.
        
        """
        z = np.zeros(len(x_values))
        j = -1.
        for x in x_values:
            j +=1
            if x < self.fungus.Smin:
                z[j] = self.age_dday
            else:
                z[j] = max(0., -(1/self.fungus.r2)*(x - self.surface_total))
                
        return z
                

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
                 basis_for_dday = -2.,
                 temp_min = 10.,
                 temp_max = 30.,
                 wd_min = 10.,
                 loss_rate = 1./120,
                 degree_days_to_chlorosis = 220.,
                 degree_days_to_necrosis = 330.,
                 degree_days_to_sporulation = 350.,
                 epsilon = 0.001,
                 r1 = 0.000132,
                 r2 = 0.0006,
                 r3 = 0.,
                 Smin = 0.03,
                 Smax = 0.3,
                 growth_rate = 0.0006,
                 rh_min = 85.,
                 rain_events_to_empty = 3,
                 production_rate = 10*0.36*6.19*10**3,
                 *args, **kwds):
        """ Parameters for septoria.
        
        :Parameters:
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
            - `r1`, `r2`, `r3` (float) : Growth rates according to lesion age (cm2.dday-1)
            - `Smin` (float): Initial size of chlorotic lesion (cm2).
            - `Smax` (float): Lesion maximum size (cm2).
            - `rh_min` (float): Minimal relative humidity for sporulation.
            - `rain_events_to_empty` (int): number of rain events to empty a sporulating ring.
            
        """
        self.name = "Septoria"
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
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.rh_min = rh_min
        self.rain_events_to_empty = rain_events_to_empty
        self.production_rate = production_rate
        # TODO : Improve this parameter. 
        # Rapilly would say : RI = 10 mm/h * fDU = 0.36 * pDr = 6.19e3 spores produced by cm2.

    def __call__(self,nbSpores = None, position = None):
        return Septoria(fungus=self,nbSpores = None, position = None)

def septoria(**kwds):
    return SeptoriaParameters(**kwds)