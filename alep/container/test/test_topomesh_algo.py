from pdb import pm
from openalea.container import Topomesh,\
                                is_flip_topo_allowed,flip_edge,\
                                is_collapse_topo_allowed,collapse_edge

m = Topomesh(3)
protected = [1]
##########################
#
#       single triangle
#
##########################
for i in xrange(3) :
    m.add_wisp(0,i)
for i in xrange(3) :
    m.add_wisp(1,i)
    m.link(1,i,i)
    m.link(1,i,(i + 1) % 3)
m.add_wisp(2,3)
for i in xrange(3) :
    m.link(2,3,i)
m.add_wisp(3,0)
m.link(3,0,3)

#flip
for eid in m.wisps(1) :
    assert not is_flip_topo_allowed(m,eid)

#collapse
for eid in m.wisps(1) :
    assert not is_collapse_topo_allowed(m,eid, protected)

##########################
#
#       double triangles
#
##########################
m.add_wisp(0,3)
m.add_wisp(1,3)
m.link(1,3,0)
m.link(1,3,3)
m.add_wisp(1,4)
m.link(1,4,1)
m.link(1,4,3)
m.add_wisp(2,2)
for eid in (0,3,4) :
    m.link(2,2,eid)
m.link(3,0,2)

#flip
for eid in (1,2,3,4) :
    assert not is_flip_topo_allowed(m,eid)

for eid in (0,) :
    assert is_flip_topo_allowed(m,eid)

#collapse
for eid in (1,2,3,4) :
    assert is_collapse_topo_allowed(m,eid, protected)

for eid in (0,) :
    assert not is_collapse_topo_allowed(m,eid, protected)

##########################
#
#       triple triangles
#
##########################
m.add_wisp(1,5)
m.link(1,5,2)
m.link(1,5,3)
m.add_wisp(2,1)
for eid in (2,3,5) :
    m.link(2,1,eid)
m.link(3,0,1)

#flip
for eid in m.wisps(2) :
    assert not is_flip_topo_allowed(m,eid)

#collapse
for eid in (0,2,3) :
    assert is_collapse_topo_allowed(m,eid, protected)

for eid in (1,4,5) :
    assert not is_collapse_topo_allowed(m,eid, protected)

##########################
#
#       tetrahedre
#
##########################
m.add_wisp(2,0)
for eid in (1,4,5) :
    m.link(2,0,eid)
m.link(3,0,0)

#flip
for eid in m.wisps(2) :
    assert not is_flip_topo_allowed(m,eid)

#collapse
for eid in m.wisps(2) :
    assert not is_collapse_topo_allowed(m,eid, protected)
