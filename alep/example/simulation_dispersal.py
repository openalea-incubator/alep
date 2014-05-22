""" Run a simulation to compare strategies of dispersal. """

# Imports for visualization
# from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer

# Imports for selection of source leaf
from alinea.alep.architecture import get_leaves
from openalea.plantgl import all as pgl
from collections import OrderedDict
import collections

# Imports for wheat
from alinea.astk.caribu_interface import *
from alinea.alep.wheat import initialize_stand
from alinea.astk.plantgl_utils import get_lai
from alinea.alep.architecture import set_properties,set_property_on_each_id, get_leaves

# Imports for disease
from alinea.alep.fungal_objects import DispersalUnit, Lesion
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.dispersal_transport import PowderyMildewWindDispersal, SeptoriaRainDispersal
from alinea.alep.disease_outputs import count_dispersal_units_by_leaf, count_dispersal_units
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.alep.protocol import disperse

def plot3d(g, 
               leaf_material = None,
               stem_material = None,
               soil_material = None,
               colors = None,
               transparencies = None):
    """
    Returns a plantgl scene from an mtg.
    """
    
    Material = pgl.Material
    Color3 = pgl.Color3
    Shape = pgl.Shape
    Scene = pgl.Scene
    
    if colors is None:
        if leaf_material is None:
            leaf_material = Material(Color3(0,180,0))
        if stem_material is None:
            stem_material = Material(Color3(0,130,0))
        if soil_material is None:
            soil_material = Material(Color3(170, 85,0))
        colors = g.property('color')

    transparencies = g.property('transparency')
    
    geometries = g.property('geometry')
    greeness = g.property('is_green')
    labels = g.property('label')
    scene = Scene()

    def geom2shape(vid, mesh, scene):
        shape = None
        if isinstance(mesh, list):
            for m in mesh:
                geom2shape(vid, m, scene)
            return
        if mesh is None:
            return
        if isinstance(mesh, Shape):
            shape = mesh
            mesh = mesh.geometry
        label = labels.get(vid)
        is_green = greeness.get(vid)
        if colors:
            if transparencies==None:
                shape = Shape(mesh, Material(Color3(* colors.get(vid, [0,0,0]) )))
            else:
                shape = Shape(mesh, Material(Color3(* colors.get(vid, [0,0,0]) ), transparency=transparencies.get(vid,0)))
        elif not greeness:
            if not shape:
                shape = Shape(mesh)
        elif label.startswith('Stem') and is_green:
            shape = Shape(mesh, stem_material)
        elif label.startswith('Leaf') and is_green:
            shape = Shape(mesh, leaf_material)
        elif not is_green:
            shape = Shape(mesh, soil_material)
        shape.id = vid
        scene.add(shape)

    for vid, mesh in geometries.iteritems():
        geom2shape(vid, mesh, scene)
    return scene

# Useful functions ############################################################
def update_plot(g, leaf_source):
    # Compute severity by leaf
    nb_dus_by_leaf = count_dispersal_units_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'nb_dus', nb_dus_by_leaf, label = 'LeafElement')
    
    # Visualization
    g = alep_colormap(g, 'nb_dus', cmap=green_yellow_red(levels=1000),
                      lognorm=False, zero_to_one=False, vmax=25)

    leaves = get_leaves(g, label='LeafElement')
    # pos_sen = g.property('position_senescence')
    for id in g:
        if not id in nb_dus_by_leaf:
            g.node(id).color = (255,255,255)
            g.node(id).transparency = 0.9
        else:
            g.node(id).transparency = 0.

    g.node(leaf_source).color = (230, 62, 218)
    g.node(leaf_source).transparency = 0.
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
    
def get_source_leaf_and_max_height(g, position='center', relative_height=2./3):
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
    distances = {vid:pgl.norm(centroids[vid]-(0,0,relative_height*zmax)) for vid in centroids}
    if position=='center':
        return min(distances.items(), key=lambda x:x[1])[0], zmax
    elif position=='border':
        return max(distances.items(), key=lambda x:x[1])[0], zmax

# Dummy lesion with dummy emission ############################################                                               
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

# Dispersal ################################################################### 
def run_dispersal(model=3, position_in_canopy='center', relative_height=2./3, position_on_leaf=0.5, age=1500., seed=3):
    # Initialize a wheat canopy
    g, wheat, domain_area, domain = initialize_stand(age=age, length=1,
                                                width=1, sowing_density=250,
                                                plant_density=250, inter_row=0.12,
                                                seed=seed)
    lai = get_lai(g.property('geometry'), domain_area)/float(1e4)

    # Define the wind direction
    wind_direction = (-1, -0.5, 0)
    set_properties(g,label = 'LeafElement', wind_direction=wind_direction)
    
    # Get source leaf and make it emit DUs
    leaf_id, zmax = get_source_leaf_and_max_height(g, position=position_in_canopy, relative_height=relative_height)
    leaf = g.node(leaf_id)
    leaf.lesions = [DummyLesion(position=[position_on_leaf,0])]
    emitter = DummyEmission()

    if model==1:
        transporter = Septo3DSplash(reference_surface=domain_area)
    elif model==2:
        transporter = PowderyMildewWindDispersal(cid=0.04, a0=45., k_beer=0.65, label='LeafElement')
    elif model>=3:
        transporter = SeptoriaRainDispersal()
    
    disperse(g, emitter, transporter, "dummy", label='LeafElement')
    update_plot(g, leaf_source=leaf_id)
    
    # Check proportion intercepted : must be low
    dus = g.property('dispersal_units')
    nb_deposited = count_dispersal_units(g)
    max_dus = max([len(du) for du in dus.itervalues()])
    min_dus = min([len(du) for du in dus.itervalues()])
    return g, lai, zmax, nb_deposited, max_dus, min_dus
