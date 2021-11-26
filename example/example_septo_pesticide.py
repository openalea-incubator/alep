""" Examples highlighting the response of septoria to pesticides """

# General useful imports
import random as rd
import numpy as np
import pandas
import matplotlib.pyplot as plt
plt.ion()

# Imports for wheat
from alinea.alep.wheat import adel_one_leaf, initialize_stand
from alinea.alep.architecture import set_properties, update_healthy_area
from alinea.astk.plant_interface import grow_canopy

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septoria import plugin_septoria
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport
from alinea.alep.growth_control import PriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.disease_outputs import SeptoRecorder
from alinea.alep.disease_outputs import plot_severity_by_leaf
septoria = plugin_septoria()

# Septo 3D
from alinea.septo3d.cycle.alep_objects import Septo3DFungus
s3d_fungus = Septo3DFungus()

def example_one_leaf(fungus = septoria,
                     nb_dus=10,
                     group_dus = False,
                     infection_controler = BiotrophDUPositionModel(),
                     growth_controler = PriorityGrowthControl(),
                     nb_steps = 1000, 
                     leaf_area = None,
                     eradicant = 0.,
                     protectant = 0.):
    # Initiate wheat MTG with one leaf and good conditions for disease
    g = adel_one_leaf()
    set_properties(g, label='LeafElement', wetness=True, 
                   temperature_sequence=np.array([22.]), rain_intensity=0., 
                   global_efficacy={'eradicant':eradicant, 'protectant':protectant})
    if leaf_area is not None:
        set_properties(g, label='LeafElement', area=leaf_area, green_area=leaf_area)

    # Initiate disease
    leaf = g.node(10)
    dispersal_units = [fungus.dispersal_unit() for i in range(nb_dus)]
    for du in dispersal_units:
        du.set_position([random.random() * leaf.length, 0])
    if group_dus:
        pos = [du.position[0] for du in dispersal_units]
        dispersal_units[0].fungus.group_dus = True
        dispersal_units[0].set_position(pos)
        dispersal_units = [dispersal_units[0]]
    leaf.dispersal_units = dispersal_units

    # Initiate output recorder
    recorder = SeptoRecorder(vids=[10], group_dus=group_dus)
    recorder.record(g,0)
    
    # Simulation loop
    for i in range(nb_steps):
        # Update wheat healthy area
        update_healthy_area(g, label = 'LeafElement')
        
        # Update dispersal units and lesions
        infect(g, 1, infection_controler, label='LeafElement')
        update(g, 1, growth_controler, senescence_model=None, label='LeafElement')
        
        # Get outputs
        recorder.record(g, date=i+1)
    recorder.get_complete_dataframe()
    
    # return g, inspector
    return g, recorder
    
def plot_lesion_states(fungus = septoria,
                       nb_dus=10,
                       group_dus = False,
                       infection_controler = BiotrophDUPositionModel(),
                       growth_controler = PriorityGrowthControl(), 
                       nb_steps = 1000, 
                       leaf_area = 5.,
                       eradicant = 0.,
                       protectant = 0.):
    g, recorder = example_one_leaf(dispersal_units = dispersal_units, 
                                   group_dus = group_dus, 
                                   nb_steps = nb_steps, 
                                   leaf_area = leaf_area,
                                   eradicant = eradicant,
                                   protectant = protectant)
    outputs = pandas.DataFrame({'disease area': recorder.data.leaf_disease_area,
                                'surface incubating': recorder.data.surface_inc,
                                'surface chlorotic': recorder.data.surface_chlo,
                                'surface necrotic': recorder.data.surface_nec,
                                'surface sporulating': recorder.data.surface_spo})
    ax = outputs.plot()
    ax.set_xlabel('Degree days')
    ax.set_ylabel('Surface (in cm2)')