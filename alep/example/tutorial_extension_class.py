""" Tutorial for powdery mildew"""
from alinea.alep.wheat_examples import adel_mtg, adel_mtg2, adel_one_leaf
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
from alinea.alep.architecture import *
from alinea.alep.disease_operation import *
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal import RandomDispersal
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.protocol import *

from openalea.vpltk import plugin

# Define a plant or canopy
# g = adel_one_leaf()
g = adel_mtg2()

# Add missing properties needed for the simulation
# The simulation requires the following properties on leaf elements:
#   - 'surface': total surface of the leaf element (in cm2)
#   - 'healthy_surface': surface of the leaf element without lesion or senescence (in cm2)
#   - 'age': age of the leaf (in hours)
#   - 'position_senescence': position of the senescence on blade axis
set_properties(g,label = 'LeafElement',
               surface=5, healthy_surface=5, 
               age = 0., position_senescence=None)
               
# discover the disease implementation for powdery mildew
diseases=plugin.discover('alep.disease')
powdery_mildew = diseases['powdery_mildew'].load()

# Modify a part of the lesion
LesionKlass = powdery_mildew.lesion()

############################################

def update_diameter_max_2(lesion, leaf=None):
    """ Compute maximum diameter of the lesion according to its age.

    Parameters
    ----------
        leaf: Leaf sector node of an MTG 
            A leaf sector with properties (e.g. healthy surface,
            senescence, rain intensity, wetness, temperature, lesions)
    """
    self = lesion
    # Parameters for the calculation
    f = self.fungus
    diameter_max = f.diameter_max
    diameter_min = f.diameter_min
    leaf_age_effect = f.leaf_age_effect
        
    # Maximum diameter of the lesion according to leaf age
    kmax = diameter_min + (diameter_max - diameter_min) * exp(-leaf_age_effect * leaf.age)+2
    self.diameter_max = kmax

LesionKlass.update_diameter_max = update_diameter_max_2

############################################
a = b = c = update_diameter_max_2

def f(LesionKlass, method=None):
    if method:
        LesionKlass.method = method


# Create a pool of dispersal units (DU)
nb_du = 2
dispersal_units = generate_stock_du(nb_du, disease=powdery_mildew)

# Distribute the DU 
inoculator = RandomInoculation()
initiate(g, dispersal_units, inoculator)

# Simulation
controler = NoPriorityGrowthControl()
dispersor = RandomDispersal()
dt = 1
nb_steps = 1000
for i in range(0,nb_steps,dt):
    print i
    set_properties(g,label = 'LeafElement',
                    wetness=True,
                    temp=22.,
                    rain_intensity=0.,
                    rain_duration=0.,
                    relative_humidity=85.,
                    wind_speed=0.2,
                    age=i+dt)
    infect(g, dt)
    update(g, dt, controler)
    
    scene = plot3d(g)
    disperse(g, scene, dispersor, "powdery_mildew")
    
plot_lesions(g)
