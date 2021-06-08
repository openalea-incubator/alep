"""" Demonstrate the dispersal algorithms"""
import collections
from alinea.adel.astk_interface import AdelWheat
from alinea.alep.fungus import Lesion, Fungus
from alinea.alep.dispersal_emission import SimpleEmission
from alinea.alep.dispersal_transport import SeptoriaRainDispersal, PowderyMildewWindDispersal, BrownRustDispersal
from alinea.alep.protocol import disperse
from alinea.alep.disease_outputs import plot3d_transparency
import openalea.plantgl.all as pgl
from alinea.alep.architecture import get_leaves

# utilities to place the source leaf
def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)
    
def get_source_leaf_and_max_height(g, position='center', relative_height=2./3):
    tesselator = pgl.Tesselator()
    bbc = pgl.BBoxComputer(tesselator)
    leaves = get_leaves(g, label='LeafElement')
    centroids = g.property('centroid')
    geometries = g.property('geometry')
    targets = list(leaf for leaf in leaves if leaf in iter(geometries.keys()))
    for vid in targets:
        if is_iterable(geometries[vid]):
            bbc.process(pgl.Scene(geometries[vid]))
        else:
            bbc.process(pgl.Scene([pgl.Shape(geometries[vid])]))
        center = bbc.result.getCenter()
        centroids[vid] = center
    zmax = max(list(centroids.items()), key=lambda x:x[1][2])[1][2]
    distances = {vid:pgl.norm(centroids[vid]-(0,0,relative_height*zmax)) for vid in centroids}
    if position=='center':
        return min(list(distances.items()), key=lambda x:x[1])[0], zmax
    elif position=='border':
        return max(list(distances.items()), key=lambda x:x[1])[0], zmax


fungus_name = "lesion_tutorial"

# Create our own emitting lesion for this example
class LesionTutoDispersal(Lesion):
    def __init__(self, mutable=False):
        super(LesionTutoDispersal, self).__init__()
        self.fungus_name = fungus_name

    def emission(*args, **kwds):
        return 10000

fungus = Fungus(Lesion=LesionTutoDispersal, parameters={"name":fungus_name})
source_lesion = fungus.lesion()

# create wheat canopy 
adel= AdelWheat(nplants=20, nsect=7)
g = adel.setup_canopy(1500)
# plot
# adel.plot(g)

# choose source leaf
leaf_ids = [id for id, label in g.property("label").items() if label.startswith("LeafElement")]
source_leaf = g.node(leaf_ids[1])

# Alternatively, pick a leaf based on ist position in the canopy
vid, hmax = get_source_leaf_and_max_height(g)
source_leaf = g.node(vid)

# inoculate this leaf
source_lesion.is_sporulating = True # Required for Popdrops usage
source_leaf.lesions = [source_lesion]

# Setup dispersal models + parameters for dispersion
emitter = SimpleEmission()
#rain 3d (slow)
transporter = SeptoriaRainDispersal()
# wind 3d (callonec model)
#transporter = PowderyMildewWindDispersal(label='LeafElement')
# wind 1D
#transporter = BrownRustDispersal(domain_area=adel.domain_area)
# missing : rain 1D (septo 3d dispersal model)

# Simulate one dispersal event
g = disperse(g, emission_model=emitter, transport_model=transporter, fungus_name=fungus_name)

# Visualize (currently only works for septoria)
# rain/ wind
transporter.plot_distri_3d(g)
# wind 1D
#transporter.plot_layers(g)







