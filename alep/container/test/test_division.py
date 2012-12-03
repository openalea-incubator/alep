from openalea.container import Topomesh,topo_divide_face,topo_divide_cell

#############################################
#
#		divide face
#
#############################################

#	0----0----1----1----2
#	|                   |
#	5        fid        2
#	|                   |
#	5----4----4----3----3

m = Topomesh(2)

for i in xrange(6) :
	m.add_wisp(0,i)

for i in xrange(6) :
	m.add_wisp(1,i)
	m.link(1,i,i)
	m.link(1,i,(i + 1) % 6)

fid = m.add_wisp(2)
for i in xrange(6) :
	m.link(2,fid,i)

fid1,fid2,eid = topo_divide_face(m,fid,1,4)

#	0----0----1----1----2
#	|         |         |
#	5  fid1  eid  fid2  2
#	|         |         |
#	5----4----4----3----3

assert set(m.borders(1,eid) ) == set([1,4])
assert set(m.regions(1,eid) ) == set([fid1,fid2])
fid1, = m.regions(1,5)
fid2, = m.regions(1,2)
assert set(m.borders(2,fid1) ) == set([0,eid,4,5])
assert set(m.borders(2,fid2) ) == set([1,2,3,eid])

#############################################
#
#		divide cell
#
#############################################

#top
#	0----12---4----16---8
#	|         |         |
#	0    0    4    4    8
#	|         |         |
#	1----13---5----17---9

#front
#	1----13---5----17---9
#	|         |         |
#	1    1    5    5    9
#	|         |         |
#	2----14---6----18---10

#bottom
#	3----15---7----19---11
#	|         |         |
#	2    2    6    6    10
#	|         |         |
#	2----14---6----18---10

#back
#	0----12---4----16---8
#	|         |         |
#	3    3    7    7    11
#	|         |         |
#	3----15---7----19---11

#left
#	0----0----1
#	|         |
#	3    8    1
#	|         |
#	3----2----2

#middle (no face)
#	4----4----5
#	|         |
#	7         5
#	|         |
#	7----6----6

#right
#	8----8----9
#	|         |
#	11   9    9
#	|         |
#	11---10---10

m = Topomesh(3)

#points
for i in xrange(12) :
	m.add_wisp(0,i)

#edges
for i in xrange(20) :
	m.add_wisp(1,i)

for i in xrange(4) :
	m.link(1,i,i)
	m.link(1,i,(i + 1) % 4)
	
	m.link(1,4 + i,4 + i)
	m.link(1,4 + i,4 + (i + 1) % 4)
	
	m.link(1,8 + i,8 + i)
	m.link(1,8 + i,8 + (i + 1) % 4)
	
	m.link(1,12 + i,i)
	m.link(1,12 + i,4 + i)
	
	m.link(1,16 + i,4 + i)
	m.link(1,16 + i,8 + i)

#faces
for i in xrange(10) :
	m.add_wisp(2,i)

for i in xrange(4) :
	m.link(2,8,i)
	m.link(2,9,8 + i)

for eid in (0,4,12,13) :
	for i in xrange(3) :
		m.link(2,i,i + eid)
		m.link(2,4 + i,4 + i + eid)

for eid in (3,12,7,15) :
	m.link(2,3,eid)
	m.link(2,7,4 + eid)

#cell
cid = m.add_wisp(3)
for i in xrange(10) :
	m.link(3,cid,i)

#divide
cid1,cid2,fid = topo_divide_cell(m,cid,(4,5,6,7) )


#middle (no face)
#	4----4----5
#	|         |
#	7   fid   5
#	|         |
#	7----6----6

assert set(m.borders(2,fid) ) == set([4,5,6,7])
assert set(m.regions(2,fid) ) == set([cid1,cid2])
cid1, = m.regions(2,8)
cid2, = m.regions(2,9)
assert set(m.borders(3,cid1) ) == set([0,1,2,3,8,fid])
assert set(m.borders(3,cid2) ) == set([4,5,6,7,9,fid])







