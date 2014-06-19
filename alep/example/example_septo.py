""" Run model of septoria """

# General useful imports
import random as rd
import numpy as np
import pandas
import matplotlib.pyplot as plt

# Imports for wheat
from alinea.alep.wheat import adel_one_leaf, initialize_stand
from alinea.alep.architecture import set_properties, update_healthy_area

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.disease_outputs import LeafInspector
from alinea.alep.disease_outputs import plot_severity_by_leaf

def example_one_leaf(dispersal_units=[plugin_septoria().dispersal_unit() for i in range(10)], group_dus=False, nb_steps=1000, leaf_area=None):
    # Initiate wheat MTG with one leaf and good conditions for disease
    g = adel_one_leaf()
    set_properties(g, label='LeafElement', wetness=True, temp=22., rain_intensity=0.)
    if leaf_area is not None:
        set_properties(g, label='LeafElement', area=leaf_area, green_area=leaf_area)

    # Initiate disease
    leaf = g.node(10)
    for du in dispersal_units:
        du.set_position([random.random() * leaf.length, 0])
    if group_dus:
        pos = [du.position[0] for du in dispersal_units]
        dispersal_units[0].fungus.group_dus = True
        dispersal_units[0].set_position(pos)
        dispersal_units = [dispersal_units[0]]
    leaf.dispersal_units = dispersal_units
    infection_controler = BiotrophDUPositionModel()
    growth_controler = NoPriorityGrowthControl()

    # Initiate output recorder
    inspector = LeafInspector(g, blade_id=8)

    # Simulation loop
    for i in range(nb_steps):
        # Update wheat healthy area
        update_healthy_area(g, label = 'LeafElement')
        
        # Update dispersal units and lesions
        infect(g, 1, infection_controler, label='LeafElement')
        update(g, 1, growth_controler, senescence_model=None, label='LeafElement')
        
        # Get outputs
        inspector.update_variables(g)
        inspector.update_du_variables(g)
    
    return g, inspector
    
def plot_lesion_states():
    g, insp = example_one_leaf(dispersal_units=[plugin_septoria().dispersal_unit(nb_rings_by_state=1)], nb_steps=1000)
    outputs = pandas.DataFrame({'disease area': insp.leaf_disease_area,
                                'surface incubating': insp.surface_inc,
                                'surface chlorotic': insp.surface_chlo,
                                'surface necrotic': insp.surface_nec,
                                'surface sporulating': insp.surface_spo})
    ax = outputs.plot()
    ax.set_xlabel('Degree days')
    ax.set_ylabel('Surface (in cm2)')


def example_plant(dispersal_units=[plugin_septoria().dispersal_unit() for i in range(10)], group_dus=False, nb_steps=1000):
    # Initiate wheat MTG with one leaf and good conditions for disease
    g, wheat, domain_area, domain = initialize_stand(age=1600., length=0.1, width=0.1,
            sowing_density=150, plant_density=150, inter_row=0.12, nsect=1, seed=3)
    set_properties(g, label='LeafElement', wetness=True, temp=22.)
    rain_data = pandas.DataFrame({'rain':[1.,1.,1.,1.]})
    
    # Initiate disease
    source_leaf = g.node(82)
    for du in dispersal_units:
        du.set_position([random.random() * source_leaf.length, 0])
    if group_dus:
        pos = [du.position[0] for du in dispersal_units]
        dispersal_units[0].fungus.group_dus = True
        dispersal_units[0].set_position(pos)
        dispersal_units = [dispersal_units[0]]
    source_leaf.dispersal_units = dispersal_units
    infection_controler = BiotrophDUPositionModel()
    growth_controler = NoPriorityGrowthControl()
    emitter = PopDropsEmission(domain=domain)
    transporter = PopDropsTransport(domain=domain, domain_area=domain_area)
    
    # Initiate output recorder
    inspectors = {}
    ind=4
    for blade in [80, 88, 96]:
            ind -= 1
            inspectors['F%d' % ind] = LeafInspector(g, blade_id=blade)
    
    # Simulation loop
    for i in range(nb_steps):
        # Update wheat healthy area
        update_healthy_area(g, label = 'LeafElement')
        
        # Update dispersal units and lesions
        infect(g, 1, infection_controler, label='LeafElement')
        update(g, 1, growth_controler, senescence_model=None, label='LeafElement')
        
        # Disperse
        if i%100==0.:
            g = disperse(g, emitter, transporter, "septoria", label='LeafElement', weather_data=rain_data)
        
        # Get outputs
        for inspector in inspectors.itervalues():
            inspector.update_variables(g)
            inspector.update_du_variables(g)

    return g, inspectors

def plot_plant_outputs(call='example_plant(dispersal_units=[plugin_septoria().dispersal_unit() for i in range(10)], group_dus=True, nb_steps=1000)'):
    exec('g, inspectors = ' + call)
    disease_area = pandas.DataFrame({k:v.leaf_disease_area for k,v in inspectors.iteritems()})
    surf_nec = pandas.DataFrame({k:v.surface_total_nec for k,v in inspectors.iteritems()})
    fig, axs = plt.subplots(1,2)
    disease_area.plot(ax=axs[0])
    surf_nec.plot(ax=axs[1])
    
g, inspectors10 = example_plant(dispersal_units=[plugin_septoria().dispersal_unit() for i in range(10)], group_dus=True, nb_steps=1000)
disease_area10 = pandas.DataFrame({k:v.leaf_disease_area for k,v in inspectors10.iteritems()})
surf_nec10 = pandas.DataFrame({k:v.surface_total_nec for k,v in inspectors10.iteritems()})

g, inspectors1 = example_plant(dispersal_units=[plugin_septoria().dispersal_unit(nb_rings_by_state=1) for i in range(10)], group_dus=True, nb_steps=1000)
disease_area1 = pandas.DataFrame({k:v.leaf_disease_area for k,v in inspectors10.iteritems()})
surf_nec1 = pandas.DataFrame({k:v.surface_total_nec for k,v in inspectors10.iteritems()})
fig, axs = plt.subplots(1,2)
disease_area10.plot(ax=axs[0])
surf_nec10.plot(ax=axs[1])
disease_area1.plot(ax=axs[0])
surf_nec1.plot(ax=axs[1])

    
def stat_profiler(call='example_plant(dispersal_units=[plugin_septoria().dispersal_unit() for i in range(10)], group_dus=False, nb_steps=1000)'):
    import cProfile
    import pstats
    cProfile.run(call, 'restats')
    p = pstats.Stats('restats')
    p.sort_stats('time').print_stats(10)