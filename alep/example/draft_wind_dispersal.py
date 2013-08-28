""" Draft preparing implementation of powdery mildew dispersal.
    
    - Generate a fully developed vine canopy.
    - List all leaves in the scene.
    - Choose a source leaf.
    - Make the source leaf emit a stock of dispersal units.
    - Define a wind direction.
    - Sort the only leaves that are in the dispersal cone.
    - Apply the formalism of dispersal unit interception from Calonnec
    et al., 2008 on each sorted leaf.
"""
# Imports ##########################################################
from alinea.alep.vine import Vine
from random import seed, choice, shuffle
from math import degrees, exp
from collections import OrderedDict
from alinea.alep.disease_operation import generate_stock_du
from alinea.alep.disease_outputs import count_dispersal_units_by_leaf
from alinea.alep.architecture import set_property_on_each_id
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer
from openalea.plantgl import all as pgl
from openalea.vpltk import plugin
from time import sleep

# Color leaves with dispersal units ################################
def update_plot(g):
    # Count dispersal units by id & add it as MTG property 
    nb_dus_by_leaf = count_dispersal_units_by_leaf(g, label = 'lf')
    set_property_on_each_id(g, 'nb_dus', nb_dus_by_leaf, label = 'lf')
    
    # Visualization
    g = alep_colormap(g, 'nb_dus', cmap=green_yellow_red(levels=10), lognorm=False, zero_to_one=False)
    labels = g.property('label')
    trunk_ids = [k for k,l in labels.iteritems() 
                if l.startswith('tronc') or l.startswith('en')]
    brown = (100,70,30)
    for id in trunk_ids:
        trunk = g.node(id)
        trunk.color = brown
    scene = plot3d(g)
    Viewer.display(scene)
    return scene

def wind_dispersal(g, source, targets, wind = (1,0,0), cid = .04, a0 = 45.):
    """
    1. compute centroids+surface for all leaves
    2. compute distance
    3. filter leaves (OAxW < 0)+ cone
    4. Compute angle
    5. Sort
    6. Compute Qc iteratively
    
    Qr = Qc S exp(-cid . d)(\alpha_0 - \alpha) / (100 \alpha_0)
    """    
    g.add_property('centroid')
    g.add_property('surface')
    geometries = g.property('geometry')
    centroids = g.property('centroid')
    surfaces = g.property('surface')
    dus = g.property('dispersal_units')
    tesselator = pgl.Tesselator()
    bbc = pgl.BBoxComputer(tesselator)
    
    def centroid(vid):
        bbc.process(pgl.Scene(geometries[vid]))
        center = bbc.result.getCenter()
        centroids[vid] = center
    
    def surface(vid):
        surfaces[vid] = pgl.surface(geometries[vid][0])*1000
    
    # Compute centroids and surface
    centroid(source)
    for vid in targets:
        centroid(vid)
        surface(vid)
        
    O = ptSrc = centroids[source]
    vects = {vid:(centroids[vid]-O) for vid in targets if (centroids[vid]-O)*wind >= 0 }
    angles = {vid:degrees(pgl.angle(vect, wind)) for vid, vect in vects.iteritems() if degrees(pgl.angle(vect, wind))<= a0 }
    
    distances = {vid:pgl.norm(vects[vid]) for vid in angles}
    
    # Sort the vid based on the distance
    distances = OrderedDict(sorted(distances.iteritems(), key=lambda x: x[1]))
    
    DUs = dus[source]
    shuffle(DUs)
    n = len(DUs)
    print '######### DU : ', n
    for leafid in distances:
        qc = n * surfaces[leafid]/100 * exp(-cid * distances[leafid]) * (a0 - angles[leafid]) / a0
        if qc < 1:
            break
        if qc > n:
            qc = n
        
        dus[leafid] = DUs[:int(qc)]
        del DUs[:int(qc)]
        print 'Leaf ', leafid, ' receive ', qc, ' spores'
        print len(DUs), ' DUs left in stock'
        print ''
        if len(DUs) < 1:
            break
        
        sleep(0.5)
        update_plot(g)
        
        
# Dispersal ########################################################
# Generate a fully developed vine canopy
vine = Vine()
g = vine.setup_canopy(age=50)

# List all leaves in the scene
labels = g.property('label')
leaves = [k for k,l in labels.iteritems() if l.startswith('lf')]

# Select the source leaf
seed(10)
source = choice(leaves)

# Make it emit a stock of dispersal units
diseases=plugin.discover('alep.disease')
powdery_mildew = diseases['powdery_mildew'].load()
nb_du = 1000
dispersal_units = generate_stock_du(nb_du, disease=powdery_mildew)
g.node(source).dispersal_units = dispersal_units

# Define the wind direction
wind = (1,0,0)

# Call dispersal function
wind_dispersal(g, source, targets=leaves, wind=wind)

# Display vine with affected leaves
# update_plot(g)
