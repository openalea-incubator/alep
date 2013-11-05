""" Run a simulation to compare strategies of dispersal. """

# Imports for visualization
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer

# Imports for selection of source leaf
from alinea.alep.architecture import get_leaves
from openalea.plantgl import all as pgl
from collections import OrderedDict
import collections

# Imports for wheat
from alinea.astk.caribu_interface import *
from alinea.alep.wheat import initialize_stand
from alinea.alep.architecture import set_properties,set_property_on_each_id, get_leaves

# Imports for disease
from alinea.alep.fungal_objects import DispersalUnit, Lesion
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.dispersal_transport import PowderyMildewWindDispersal, SeptoriaRainDispersal
from alinea.alep.disease_outputs import count_dispersal_units_by_leaf
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.alep.protocol import disperse

# Useful functions ############################################################
def update_plot(g, leaf_source):
    # Compute severity by leaf
    nb_dus_by_leaf = count_dispersal_units_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'nb_dus', nb_dus_by_leaf, label = 'LeafElement')
    
    # Visualization
    g = alep_colormap(g, 'nb_dus', cmap=green_yellow_red(levels=1000),
                      lognorm=False, zero_to_one=False, vmax=10)

    leaves = get_leaves(g, label='LeafElement')
    # pos_sen = g.property('position_senescence')
    for id in g:
        if not id in nb_dus_by_leaf:
            g.node(id).color = (0,0,0)

    g.node(leaf_source).color = (230, 62, 218)
    scene = plot3d(g)
    Viewer.display(scene)
    return scene

def periodise_canopy(g, domain):
    geometries = g.property('geometry')
    shapes = [geom2shape(k,v) for k,v in geometries.iteritems()]                                          
    cs = CaribuScene(scene=shapes, pattern=domain)
    cs.runPeriodise()
    shapes = cs.generate_scene()
    newgeom = dict([(s.id,s.geometry) for s in shapes])
    geometries.update(newgeom)
    return g

def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)

def get_source_leaf(g, position='center'):
    tesselator = pgl.Tesselator()
    bbc = pgl.BBoxComputer(tesselator)
    leaves = get_leaves(g, label='LeafElement')
    centroids = g.property('centroid')
    geometries = g.property('geometry')    
    targets = list(leaf for leaf in leaves if leaf in geometries.iterkeys())
    for vid in targets:
        if is_iterable(geometries[vid]):
            bbc.process(pgl.Scene(geometries[vid]))
        else:
            bbc.process(pgl.Scene([pgl.Shape(geometries[vid])]))
        center = bbc.result.getCenter()
        centroids[vid] = center
    zmax = max(centroids.items(), key=lambda x:x[1][2])[1][2]
    distances = {vid:pgl.norm(centroids[vid]-(0,0,zmax/2.)) for vid in centroids}
    if position=='center':
        return min(distances.items(), key=lambda x:x[1])[0]
    elif position=='border':
        return max(distances.items(), key=lambda x:x[1])[0]
    
# Dispersal ################################################################### 

# Initialize a wheat canopy
# g, wheat, domain_area, domain = initialize_stand(age=1500., length=1,
                                                 # width=1, sowing_density=150,
                                                 # plant_density=150, inter_row=0.12)
                                       
g, wheat, domain_area, domain = initialize_stand(age=1500., length=1,
                                                width=1, sowing_density=250,
                                                plant_density=250, inter_row=0.12)

# periodise_canopy(g, domain)
                                                
# Create a undetermined lesion emmitting a stock of dispersal units
class DummyLesion(Lesion):
    """ Undetermined lesion with only a method emission."""
    def __init__(self, nb_spores=None, position=None):
        super(DummyLesion, self).__init__(nb_spores=nb_spores, position=position)
        class params():
            def __init__(self, name="dummy"):
                self.name = name
        self.fungus = params()
        # Stock of spores needs to be positive
        # TODO : change this condition to a method 'can_emit'
        self.stock_spores = 1
        
    def is_sporulating(self):
        return True

class DummyEmission():        
    def get_dispersal_units(self, g, fungus_name="dummy", label='LeafElement'):
        DU={}
        lesions = {k:[l for l in les if l.fungus.name is fungus_name and l.is_sporulating()] 
                    for k, les in g.property('lesions').iteritems()} 
        for vid, l in lesions.iteritems():
            for lesion in l:
                emissions = [DispersalUnit(nb_spores=lesion.nb_spores, status='emitted',
                             position=lesion.position) for i in range(int(1e4))]
                try:
                    DU[vid] += emissions
                except:
                    DU[vid] = emissions
        return DU

        
leaf_id = get_source_leaf(g, position='center')
leaf = g.node(leaf_id)
leaf.color = (0,0,180)
leaf.lesions = [DummyLesion(position=[0.5,0])]
emitter = DummyEmission()

transporter = Septo3DSplash(reference_surface=domain_area)
transporter2 = PowderyMildewWindDispersal(cid=0.04, a0=45., k_beer=0.65, label='LeafElement')
transporter3 = SeptoriaRainDispersal()

# Define the wind direction
wind_direction = (-1, -0.5, 0)
set_properties(g,label = 'LeafElement', wind_direction=wind_direction)

disperse(g, emitter, transporter3, "dummy", label='LeafElement')

# scene = plot3d(g)
# Viewer.display(scene)

# plot_dispersal_units(g)
update_plot(g, leaf_source=leaf_id)

# Check proportion intercepted : must be low
