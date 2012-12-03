execfile("test_division.py")

from vplants.plantgl.math import Vector3
from vplants.plantgl.scenegraph import Material,Text,Font,Shape
from openalea.pglviewer import SceneView,display
from openalea.tissueview import MeshView
from openalea.tissueshape import centroid

pos = {}
for i,tup in enumerate([(1,1),(-1,1),(-1,-1),(1,-1)]) :
	pos[i] = Vector3(tup,-2)
	pos[4 + i] = Vector3(tup,0)
	pos[8 + i] = Vector3(tup,2)

mv = MeshView(m,pos,1,10,Material( (0,0,0) ) )
mv.redraw()

tmp = MeshView(m,pos,3,20,Material( (0,255,255) ) )
tmp.redraw()
mv.merge(tmp._scene)

font = Font("ariana",8)
mat = Material( (0,0,0) )
for pid,vec in pos.iteritems() :
	geom = Text("%d" % pid,vec,False,font)
	mv.add(Shape(geom,mat) )

mat = Material( (0,0,255) )
for eid in m.wisps(1) :
	geom = Text("%d" % eid,centroid(m,pos,1,eid),False,font)
	mv.add(Shape(geom,mat) )

mat = Material( (0,255,0) )
for fid in m.wisps(2) :
	geom = Text("%d" % fid,centroid(m,pos,2,fid),False,font)
	mv.add(Shape(geom,mat) )

mat = Material( (255,0,0) )
for cid in m.wisps(3) :
	geom = Text("%d" % cid,centroid(m,pos,3,cid),False,font)
	mv.add(Shape(geom,mat) )

display(mv)

