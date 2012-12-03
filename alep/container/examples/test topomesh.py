from openalea.container import Topomesh

m=Topomesh(2)

pts=[m.add_wisp(0) for i in xrange(6)]
edges=[m.add_wisp(1) for i in xrange(7)]
faces=[m.add_wisp(2) for i in xrange(2)]

for i,(pid1,pid2) in enumerate([(0,1),(1,2),(2,3),(3,4),(4,5),(5,1),(1,4)]) :
    m.link(1,edges[i],pts[pid1])
    m.link(1,edges[i],pts[pid2])

for i,eids in enumerate([(0,4,5,6),(1,2,3,6)]) :
    for eid in eids :
        m.link(2,faces[i],edges[eid])

for eid in edges :
    print eid
    print "B",list(m.borders(1,eid)),m.nb_borders(1,eid)
    print "R",list(m.regions(1,eid)),m.nb_regions(1,eid)

for fid in faces :
    print list(m.borders(2,fid)),m.nb_borders(2,fid)
