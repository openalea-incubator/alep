""" Test basic responses of the grapevine powdery mildew model on wheat. """

# Imports #########################################################################
import random as rd
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

from alinea.alep.wheat_examples import adel_mtg, adel_mtg2, adel_one_leaf
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

from alinea.alep.powdery_mildew import *
from alinea.alep.architecture import *
from alinea.alep.disease_operation import *
from alinea.alep.du_position_checker import BiotrophDUProbaModel
from alinea.alep.dispersal_emission import PowderyMildewWindEmission
from alinea.alep.dispersal_transport import RandomDispersal
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.inoculation import RandomInoculation

from alinea.alep.protocol import *

# Display #########################################################################
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
    
# Tests ###########################################################################    
def test_initiate():
    """ Check if 'initiate' from 'protocol.py' deposits dispersal units on the MTG.
    
    """
    # Generate a wheat MTG
    g = adel_mtg2()
    set_initial_properties_g(g)
    # Generate a stock of powdery mildew dispersal units and distribute it on g
    nb_dus_in_stock = 100
    distribute_dispersal_units(g, nb_dus=nb_dus_in_stock)
    # Display dispersal units
    plot_DU(g)
    return g
    
def test_infect():
    """ Check if 'infect' from 'protocol.py' leads to infection by dispersal units on the MTG.

    """
    # Generate a wheat MTG
    g = adel_mtg2()
    set_initial_properties_g(g)
    # Generate a stock of septoria dispersal units and distribute it on g
    nb_dus_in_stock = 100
    distribute_dispersal_units(g, nb_dus=nb_dus_in_stock)
    # Run a loop of simulation
    dt = 1
    nb_steps = 100
    # Display dispersal units ...
    plot_DU(g)
    for i in range(nb_steps):
        update_climate_all(g)
        update_leaf_age(g, dt)
        infect(g, dt)
    # ... then display lesions       
    plot_lesions(g)
    return g
    
def test_update():
    """ Compute lesion growth and draw surface evolution.
    
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    set_properties(g,label = 'LeafElement',
                    surface=5., 
                    healthy_surface=5., 
                    age = 0.,
                    position_senescence=None)    
    # Attach a lesion on g
    nb_lesion = 1
    distribute_disease(g,
                       fungal_object='lesion', 
                       nb_objects=nb_lesion, 
                       disease_model='powdery_mildew',
                       initiation_model=RandomInoculation())
    # Call model of growth control
    controler = NoPriorityGrowthControl()
    # Loop of simulation
    dt = 1
    nb_steps = 1000
    diameter = np.zeros(nb_steps)
    surface = np.zeros(nb_steps)
    surface2 = np.zeros(nb_steps)
    for i in range(0,nb_steps,dt):
        set_properties(g,label = 'LeafElement',
                        age = i+dt,
                        wetness=True,
                        temp=22.,
                        rain_intensity=0.,
                        rain_duration=0.,
                        relative_humidity=85.,
                        wind_speed=0.2)
        update(g, dt, controler)
        
        lesions = g.property('lesions')
        l = lesions.values()[0][0]
        diameter[i] = l.diameter
        surface[i] = l.surface
        assert round(surface[i],6) == round(pi*(diameter[i]**2)/4,6)
    
    # Display results
    fig = plt.figure()
    
    ax1 = fig.add_subplot(2,1,1)
    plot(range(0,nb_steps,dt), diameter)
    ylabel('Lesion diameter')
    
    ax2 = fig.add_subplot(2,1,2)
    plot(range(0,nb_steps,dt), surface, color='r')
    ylabel('Lesion surface')
    xlabel('Simulation time step')
       
    fig.subplots_adjust(hspace=1)
    
    plt.show()
        
    return g
    
def test_growth_control():
    """ Check if 'control_growth' from 'protocol.py' limits the lesion growth
        up to available surface on the leaf.
    
    Generate a wheat MTG and deposit sumultaneously 1000 lesions on a leaf element.
    Run a loop to compute update. Check that the healthy surface of the leaf
    decreases up to 0. Check that lesion growth is stopped after this point.
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    initial_leaf_surface = 5.
    set_properties(g,label = 'LeafElement',
                    surface=initial_leaf_surface, 
                    healthy_surface=initial_leaf_surface, 
                    age = 0.,
                    position_senescence=None)  
    total_initial_surface = sum(g.node(v).surface for v in g if g.label(v).startswith('LeafElement'))
    # Distribute lesions on MTG
    nb_lesion = 1000
    distribute_disease(g,
                       fungal_object='lesion', 
                       nb_objects=nb_lesion, 
                       disease_model='powdery_mildew',
                       initiation_model=RandomInoculation())
    # Call the model to check if DUs can infect where they are
    position_checker = BiotrophDUProbaModel()
    # Call models of growth control and senescence
    controler = NoPriorityGrowthControl()
    # Loop of simulation
    dt = 1
    nb_steps = 150
    surface = np.zeros(nb_steps)
    # Healthy surface the day before
    healthy_surface_before = []
    for i in range(0,nb_steps,dt):
        set_properties(g,label = 'LeafElement',
                        age = i+dt,
                        wetness=True,
                        temp=22.,
                        rain_intensity=0.,
                        rain_duration=0.,
                        relative_humidity=85.,
                        wind_speed=0.2)
        infect(g, dt)
        update(g, dt, controler)
        
        # Find the value of interest on the MTG (total healthy surface of the leaf)
        lesions = g.property('lesions')
        healthy_surfaces = g.property('healthy_surface')
        labels = g.property('label')
        
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith('LeafElement')]
            leaf_healthy_surface = sum(healthy_surfaces[lf] for lf in leaf)
        
        # Check that healthy surface + lesion surface = initial healthy surface            
        if lesions:
            leaf_lesion_surface = sum(l.surface for les in lesions.itervalues() for l in les)
            assert round(leaf_healthy_surface,6)+round(leaf_lesion_surface,6)==round(total_initial_surface,6)
        
        # Check that healthy surface decreases after emergence of the first lesion
        if lesions and not healthy_surface_before:
            # Emergence of the first lesion
            healthy_surface_before.append(leaf_healthy_surface)
        elif lesions and healthy_surface_before[0] > 0.:
            # Check that healthy surface decreases
            assert leaf_healthy_surface < healthy_surface_before
            # Update of 'healthy_surface_before' for next step
            healthy_surface_before[0] = leaf_healthy_surface
        elif lesions and healthy_surface_before[0] == 0.:
            # Check that healthy surface stays null
            assert leaf_healthy_surface == 0.
            # Update of 'healthy_surface_before' for next step
            healthy_surface_before[0] = leaf_healthy_surface
         
        # Output
        surface[i] = leaf_healthy_surface
        
    # Display results
    plot(range(0,nb_steps,dt), surface)
    ylabel('Leaf healthy surface (cm2)')
    xlabel('Hours of simulation')
    show()
    return g

def test_disperse():
    """ Check if 'disperse' from 'protocol.py' disperse new 
        dispersal units on the MTG.
    
    Generate a wheat MTG and distribute lesions randomly on leaf elements.
    Run a loop to compute update and dispersal with a constant wind.
    Check that the number of lesions on the MTG increases.
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    set_properties(g,label = 'LeafElement',
                    surface=5., 
                    healthy_surface=5., 
                    age = 0.,
                    position_senescence=None)
    # Distribute lesions on MTG
    nb_lesion = 10
    distribute_disease(g,
                       fungal_object='lesion', 
                       nb_objects=nb_lesion, 
                       disease_model='powdery_mildew',
                       initiation_model=RandomInoculation())
    # Call the model to check if DUs can infect where they are
    position_checker = BiotrophDUProbaModel()
    # Call a model of growth control and a model of dispersal
    controler = NoPriorityGrowthControl()
    emitter = PowderyMildewWindEmission()
    transporter = RandomDispersal()
    # Loop of simulation
    dt = 1
    nb_steps = 1000
    for i in range(0,nb_steps,dt):
        print(i)
        set_properties(g,label = 'LeafElement',
                        age = i+dt,
                        wetness=True,
                        temp=22.,
                        rain_intensity=0.,
                        rain_duration=0.,
                        relative_humidity=85.,
                        wind_speed=0.5)
        
        # Run protocols
        infect(g, dt, position_checker)
        update(g, dt, controler)
        
        # Count objects on the MTG before dispersal event
        lesions = g.property('lesions')
        lesions = ([l for les in lesions.values() for l in les 
                    if l.status==l.fungus.SPORULATING and l.production_is_active])
        dispersal_units = g.property('dispersal_units')
        total_stock_spores_before = sum(l.stock_spores for l in lesions)
        total_DUs_before = sum(len(du) for du in dispersal_units.itervalues())
        
        # Dispersal event
        scene = plot3d(g)
        disperse(g, emitter, transporter, "powdery_mildew")
        
        # Count objects on the MTG after dispersal event
        lesions = g.property('lesions')
        lesions = ([l for les in lesions.values() for l in les 
                    if l.status==l.fungus.SPORULATING and l.production_is_active])
        dispersal_units = g.property('dispersal_units')
        total_stock_spores_after = sum(l.stock_spores for l in lesions)
        total_DUs_after = sum(len(du) for du in dispersal_units.itervalues())
        
        if lesions:
            # Check that stocks of spores on lesions decrease
            assert total_stock_spores_after <= total_stock_spores_before
            # Check that new DUs are deposited on the MTG
            assert total_DUs_after >= total_DUs_before
            
         # Find the value of interest on the MTG (total healthy surface of the leaf)
        lesions = g.property('lesions')
        healthy_surfaces = g.property('healthy_surface')
        labels = g.property('label')
        
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith('LeafElement')]
            leaf_healthy_surface = sum(healthy_surfaces[lf] for lf in leaf)
        print(leaf_healthy_surface)
                
    return g
