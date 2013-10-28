from alinea.alep.wheat import initialize_stand
from alinea.astk.caribu_interface import *
from openalea.plantgl.all import Viewer

g, wheat, domain_area, domain = initialize_stand(age=1500., length=0.1,
                                                 width=0.2, sowing_density=150,
                                                 plant_density=150, inter_row=0.12)

geometries = g.property('geometry')
shapes = [geom2shape(k,v) for k,v in geometries.iteritems()]                                          
cs = CaribuScene(scene=shapes, pattern=domain)
cs.runPeriodise()
shapes = cs.generate_scene()
newgeom = dict([(s.id,s.geometry) for s in shapes])
geometries.update(newgeom)
Viewer.display(shapes)
