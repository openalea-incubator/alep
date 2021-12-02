""" Test the protocol of communication between adel and alep """

# Imports #########################################################################

import random as rd
import numpy
import pandas
from pylab import *

from alinea.adel.data_samples import adel_two_metamers, adel_one_leaf
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

from alinea.alep import septoria
from alinea.alep.septoria import *
from alinea.alep import powdery_mildew

from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.dispersal_transport import RandomDispersal
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.architecture import set_properties, update_healthy_area

from alinea.alep.protocol import *
from alinea.alep.diseases import get_disease
from alinea.alep.disease_operation import generate_stock_du


from datetime import datetime
import time
from math import ceil, sqrt

# Climate #########################################################################
def set_initial_properties_g(g, 
                             surface_leaf_element=5.,
                             position_senescence=None,
                             label = 'LeafElement'):
    """ Give initial values for plant properties of each LeafElement. 
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    surface: float
        Initial surface of each leaf element
    position_senescence: float
        Position of senescence on blade axis
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy
    """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.surface = surface_leaf_element
        n.healthy_surface = surface_leaf_element # TODO : Manage properly
        n.position_senescence = position_senescence
        
    return g

def update_climate_all(g, wetness=True,
                          temp = 22.,
                          rain_intensity=0.,
                          rain_duration=0.,
                          relative_humidity=85.,
                          label = 'LeafElement'):
    """ Simulate an environmental program.
    
    All leaf elements have the same values for all variables.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    wetness: bool
        True if the leaf element is wet, False otherwise
    temp: float
        Temperature of the leaf element (degrees celsius)
    rain_intensity : float
        Rain intensity on the leaf element (mm/h)
    rain_duration: float
        Rain duration (in hours)
    relative_humidity : float
        Relative humidity on the leaf element (percent)
    label: str
        Label of the part of the MTG concerned by the calculation ('LeafElement')
     
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy
    """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.wetness = wetness
        n.temp = temp
        n.rain_intensity = rain_intensity
        n.rain_duration = rain_duration
        n.relative_humidity = relative_humidity
    
    return g

# Fungus ##########################################################################
def distribute_dispersal_units(g, nb_dus=1, model="septoria"):
    """ Distribute new dispersal units on g. 
    
    Call the method 'initiate' from the protocol with dispersal units.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    nb_dus: int
        Number of dispersal units to put on g
        
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy
    """
    dispersal_units = generate_stock_du(nb_dus, model)

    inoculator = RandomInoculation()
    initiate(g, dispersal_units, inoculator)
    
    return g
    
# def distribute_lesions(g, nb_lesions=1, model="SeptoriaWithRings"):
#     """ Distribute new lesions on g.
#
#     Call the method 'initiate' from the protocol with lesions.
#
#     Parameters
#     ----------
#     g: MTG
#         MTG representing the canopy
#     nb_lesions: int
#         Number of lesions to put on g
#     model: str
#         Type of model of septoria lesion
#
#     Returns
#     -------
#     g: MTG
#         Updated MTG representing the canopy
#     """
#     fungus = septoria.plugin_septoria(model=model)
#     models = ({"SeptoriaExchangingRings":SeptoriaExchangingRings,
#                     "SeptoriaWithRings":SeptoriaWithRings,
#                     "ContinuousSeptoria":ContinuousSeptoria})
#     if model in models:
#         models[model].fungus = fungus
#         lesions = [models[model](nb_spores=rd.randint(1,100)) for i in range(nb_lesions)]
#
#     inoculator = RandomInoculation()
#     initiate(g, lesions, inoculator)
#
#     return g
    
# Call models for disease #########################################################


def inoculator():            
    """ Instantiate the class RandomInoculation().
    
    """
    inoculator = RandomInoculation()
    return inoculator

def dispersor():            
    """ Instantiate the class RandomInoculation().
    
    """
    dispersor = RandomDispersal()
    return dispersor
    
def washor():            
    """ Instantiate the class RandomInoculation().
    
    """
    washor = RapillyWashing()
    return washor

# Display #########################################################################
def temp_plot3D(g):
    """ plot g """
    scene = plot3d(g)
    Viewer.display(scene)

def plot_DU(g):
    """ plot the plant with elements carrying dispersal units in yellow """
    green = (0,180,0)
    yellow = (247, 220, 17)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'dispersal_units' in n.properties() and n.dispersal_units:
            n.color = yellow
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)
    
def plot_lesions(g):
    """ plot the plant with infected elements in red """
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)

def plot_DU_deposited(g):
    """ plot the plant with elements carrying deposited dispersal units in yellow """
    green = (0,180,0)
    yellow = (247, 220, 17)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        ind_deposited = 0.
        if 'dispersal_units' in n.properties():
            for du in n.dispersal_units:
                if du.status == 'deposited':
                    ind_deposited += 1
            
        if ind_deposited > 0 :
            n.color = yellow
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)
    
def plot_lesions_after_DU(g):
    """ plot the plant with :
        - elements carrying dispersal units in yellow
        - infected elements in red
        
    """
    green = (0,180,0)
    yellow = (247, 220, 17)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            n.color = red
        else :
            if 'dispersal_units' in n.properties():
                n.color = yellow
            else : 
                n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)

def plot_lesions_in_state(g, state):
    """ plot the plant with elements carrying lesions in given state in red """
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            for les in n.lesions:
                if les.status == state:
                    n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)

def count_lesions_in_state(g, state):
    """ count lesions of the mtg in give state.
    
    """
    lesions = g.property('lesions')
    return sum(1 for l in lesions.values() for lesion in l if lesion.status == state)
    
def count_DU(g):
    """ count DU of the mtg.
    
    """
    dispersal_units = g.property('dispersal_units')
    return sum(len(du) for du in dispersal_units.values())

def count_lesions(g):
    """ count lesions of the mtg.
    
    """
    lesions = g.property('lesions')
    return sum(len(l) for l in lesions.values())
    
class DisplayLesions(object):
    """ Print the ID of Leaf Elements where new lesions appear. """
    
    def __init__(self):
        self.old_lesions = []
        
    def print_new_lesions(self, g):
        lesions = g.property('lesions')
        
        for vid, l in lesions.items():
            for lesion in l:
                if not lesion in self.old_lesions:
                    self.old_lesions.append(lesion)
                    print(('New Lesions on : ' + g.label(vid) + ' %d' % vid))
        
    def print_all_lesions(self, g):
        from pprint import pprint
        lesions = g.property('lesions')
        ldict = {}
        for vid, l in lesions.items():
            ldict[vid] = len(l)
        pprint(ldict)
        # print('You can find lesions on LeafElements : ' + str(llist).strip('[]'))
        
    def print_lesion_surfaces(self, g):
        from pprint import pprint
        lesions = g.property('lesions')
        ldict = {}
        for vid, l in lesions.items():
            for lesion in l:
                if vid not in ldict:
                    ldict[vid] = 0
                ldict[vid] += lesion.surface
        print(('\n' + 'Sum of lesion surfaces by leaf element : ' + '\n'))             
        pprint(ldict)
    
# Tests ###########################################################################    
def test_initiate():
    """ Check if 'initiate' from 'protocol.py' deposits dispersal units on the MTG.
    
    """
    g = adel_two_metamers()
    set_initial_properties_g(g, surface_leaf_element=5.)
    fungus = get_disease('septoria')
    stock = generate_stock_du(100, fungus)
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    plot_DU(g)
    return g
    
def test_infect():
    """ Check if 'infect' from 'protocol.py' leads to infection by dispersal units on the MTG.

    """
    g = adel_two_metamers()
    set_initial_properties_g(g, surface_leaf_element=5.,
                   temperature_sequence = [18] * 24,
                    relative_humidity_sequence = [85] * 24)
    fungus = get_disease('septoria')
    stock = generate_stock_du(100, fungus)
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    dt = 1
    nb_steps = 100
    plot_DU(g)
    for i in range(nb_steps):
        update_climate_all(g)
        infect(g, dt)
            
    plot_lesions(g)
    return g
       
def test_update():
    """ Check if 'update' from 'protocol.py' provokes the growth of a lesion instantiated on the MTG.

    """
    g = adel_two_metamers()
    set_properties(g, area=20., green_area=20.)
    fungus = get_disease('septoria')
    stock = generate_stock_du(100, fungus)
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    controler = NoPriorityGrowthControl()
    
    dt = 1
    nb_steps = 750
    nb_les_inc = numpy.zeros(nb_steps)
    # nb_les_chlo = numpy.array([0. for i in range(nb_steps)])
    # nb_les_nec = numpy.array([0. for i in range(nb_steps)])
    nb_les_spo = numpy.zeros(nb_steps)
    # nb_les_empty = numpy.array([0. for i in range(nb_steps)])
    nb_les = 0.
    for i in range(nb_steps):

        ts = time.clock()
        #print('time step %d' % i)
        
        update_climate_all(g)
            
        #grow(g)
        update_healthy_area(g, label = 'LeafElement')
        infect(g, dt)
        update(g,dt, growth_control_model=controler)
        # control_growth(g, controler)
        
        # Count of lesions :
        nb_les_inc[i] = count_lesions_in_state(g, state = 0)
        # nb_les_chlo[i] = count_lesions_in_state(g, state = 1)
        # nb_les_nec[i] = count_lesions_in_state(g, state = 2)
        nb_les_spo[i] = count_lesions_in_state(g, state = 3)
        # nb_les_empty[i] = count_lesions_in_state(g, state = 4)

        te = time.clock()
        #print "time ", i, " : ", te-ts
                
    # Display results
    plot(nb_les_inc)
    plot(nb_les_spo)
    ylabel('Nombre de lesions dans cet etat sur le MTG')
    xlabel('Pas de temps de simulation')
    ylim([0, 120])
    
    # displayer = DisplayLesions()
    # displayer.print_all_lesions(g)
  
    # plot_lesions(g)
    return g

def test_disperse():
    """ Check if 'disperse' from 'protocol.py' disperse new dispersal units on the MTG.

    """
    g = adel_two_metamers()
    set_initial_properties_g(g, surface_leaf_element=5.,
                   temperature_sequence = [18] * 24,
                    relative_humidity_sequence = [85] * 24)
    fungus = get_disease('septoria')
    stock = generate_stock_du(100, fungus)
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    controler = NoPriorityGrowthControl()
    
    dt = 1
    nb_steps = 750
    nb_les = numpy.array([0 for i in range(nb_steps)])
    for i in range(nb_steps):
        print(('time step %d' % i))
        
        # Update climate and force rain occurences       
        if i>400 and i%100 == 0:
            global_rain_intensity = 4.
        else:
            global_rain_intensity = 0.
        update_climate_all(g, wetness=True, temp=22., rain_intensity = global_rain_intensity*0.75)
        
        # grow(g)
        infect(g, dt)
        update(g,dt, growth_control_model=controler)
        # control_growth(g, controler)
        
        if global_rain_intensity != 0.:
            scene = plot3d(g)
            disperse(g, scene, dispersor(), "Septoria")
               
        # Count of lesions :
        nb_les[i] = count_lesions(g)
        
        # Display results
        plot_lesions(g)
        
        # displayer = DisplayLesions()
        # displayer.print_new_lesions(g)
    
    plot(nb_les)
    ylabel('Nombre de lesions sur le MTG')
    xlabel('Pas de temps de simulation')
    show()
    
    return g

def test_washing():
    """ Check if 'washing' from 'protocol.py' washes dispersal units out of the MTG.

    """
    g = adel_two_metamers()
    set_initial_properties_g(g, surface_leaf_element=5.)
    fungus = get_disease('septoria')
    stock = generate_stock_du(100, fungus)
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    washor = RapillyWashing()
    
    dt = 1
    nb_steps = 100
    nb_DU = numpy.array([0. for i in range(nb_steps)])
    for i in range(nb_steps):
        # Update climate and force rain occurences       
        if i>2 and i%5 == 0 or (i-1)%5 == 0:
            global_rain_intensity = 4.
        else:
            global_rain_intensity = 0.
        update_climate_all(g, wetness=True, temp=22., rain_intensity = global_rain_intensity*0.75)
                
        # Compute washing 
        # (needs to be done even if no rain to update variables in the washing model)
        wash(g, washor, global_rain_intensity, DU_status='deposited')
           
        # Count of DU :
        nb_DU[i] = count_DU(g)

        # Display results
        if global_rain_intensity != 0. :
            print('\n')
            print('   _ _ _')
            print('  (pluie)')
            print(' (_ _ _ _)')
            print('   |  |')
            print('    |  | ')
            print('\n')
            print(('Sur le MTG il y a %d DU actives en tout' % nb_DU[i]))

    # Display results
    plot(nb_DU)
    ylim([0, 120])
    ylabel('Nombre de DU sur le MTG')
    xlabel('Pas de temps de simulation')
    show()
    
    return g

def test_growth_control():
    g = adel_two_metamers()
    set_initial_properties_g(g, surface_leaf_element=5.,
                   temperature_sequence = [18] * 24,
                    relative_humidity_sequence = [85] * 24)
    fungus = get_disease('septoria')
    stock = generate_stock_du(100, fungus)
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    controler = NoPriorityGrowthControl()
    
    # compute infection:
    nb_steps_inf = 11  
    for i in range(nb_steps_inf):
        print(('time step %d' % i))       
        update_climate_all(g)           
        #grow(g)
        infect(g, dt=1.)
    
    lesions = g.property('lesions')
    if lesions:
        # compute competition
        dt = 50
        nb_steps = 500
        dates = list(range(0,nb_steps,dt))
        sum_surface = numpy.array([0. for i in dates])
        for i in dates:
            print(('time step %d' % i))       
            update_climate_all(g)
                
            #grow(g)      
            update(g,dt, growth_control_model=controler)
            # control_growth(g, controler)
            
            vids = [v for v in g if g.label(v).startswith("LeafElement")]
            count_surf = 0.
            for v in vids:
                leaf = g.node(v)
                count_surf += leaf.healthy_surface
                # if 'lesions' in leaf.properties():
                    # count_surf += sum([l.surface for l in leaf.lesions])
                    # print('leaf element %d - ' % v + 'Healthy surface : %f'  % leaf.healthy_surface) 
            
            index = i/dt
            sum_surface[index] = count_surf
            
        # Display results:
        plot(dates, sum_surface)
        ylabel('Surface saine totale sur le MTG')
        xlabel('Pas de temps de simulation')
        show()
        
    return g

def test_all(model="SeptoriaExchangingRings"):
    # Generate a MTG with required properties :
    g = adel_two_metamers()
    set_initial_properties_g(g, surface_leaf_element=5.)
    
    # Deposit first dispersal units on the MTG :
    fungus = get_disease('septoria')
    stock = generate_stock_du(100, fungus)
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()
    position_checker = BiotrophDUProbaModel()
   
    # Prepare the simulation loop
    dt = 1
    nb_steps = 750
    nb_max_les = 0.
    nb_les = 0.
    for i in range(0,nb_steps,dt):
        # Update climate and force rain occurences       
        if i>400 and i%100 == 0:
            global_rain_intensity = 4.
        else:
            global_rain_intensity = 0.
        update_climate_all(g, wetness=True, temp=22., rain_intensity = global_rain_intensity*0.75)
        
        # grow(g)
        infect(g, dt, position_checker)
        update(g,dt, growth_control_model=controler)
        
        if global_rain_intensity != 0.:
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria")            
            wash(g, washor, global_rain_intensity, DU_status='deposited')
        
        # Count how many lesions are simultaneously active on the MTG at maximum charge
        nb_les = count_lesions(g)
        if nb_les > nb_max_les:
            nb_max_les = nb_les
    
    print(('max number lesions %d' % nb_max_les))
    
    return g

# if __name__ == '__main__':
    # g=test_all()
