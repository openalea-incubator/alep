""" Estimate the distribution of wind speed within a grapevine canopy """

from alinea.alep.vine import Vine
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer
from alinea.astk.plantgl_utils import *
from alinea.weather.mini_models import wind_speed_on_leaf
from alinea.alep.architecture import get_leaves, set_property_on_each_id
from alinea.alep.alep_color import alep_colormap, green_yellow_red
    
import numpy as np

# Color leaves according to wind speed at their height #############
def update_plot(g, global_wind_speed=5.):
    g = alep_colormap(g, 'wind_speed', cmap='jet', 
                      lognorm=False, zero_to_one=False, vmax=global_wind_speed)
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
    
# Wind on canopy ###################################################
# Generate a fully developed vine canopy
vine = Vine()
g = vine.setup_canopy(age=50)

# Force a value for wind speed
global_wind_speed = 5.

# Estimate the distribution of wind speed within the canopy
leaves = get_leaves(g, label='lf')
# geometries = {k:geom for k, geom_list in g.property('geometry').iteritems() for 
                # geom in geom_list if k in leaves}
geometries = {k:geom for k, geom in g.property('geometry').iteritems() if k in leaves}
lai = get_lai(geometries, vine.domain_area)
h = get_height(geometries)
surf,_=get_area_and_normal(geometries)
heights = dict([(k,np.average(h[k],weights= surf[k])) for k in h])
wind_speeds = {}
for leaf in leaves:
    wind_speeds.update({leaf:wind_speed_on_leaf(wind_speed=global_wind_speed, 
                                                leaf_height=heights[leaf],
                                                canopy_height=max(heights.itervalues()),
                                                lai=2.)})
# lai=get_lai(geometries)                                                
set_property_on_each_id(g, 'wind_speed', wind_speeds, label = 'lf')

# Visualization
update_plot(g, global_wind_speed=global_wind_speed)