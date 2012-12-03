from openalea.container import Graph

g=Graph()
sid=g.add_vertex()
tid=g.add_vertex(10)
tid2=g.add_vertex()
tid3=g.add_vertex(None)

print "add vertex",sid
print "add vertex",tid
print "add vertex",tid2

eid=g.add_edge(sid,tid)
eid2=g.add_edge(sid,tid2,34)
eid3=g.add_edge(sid,tid3,None)
print "add edge",eid,g.source(eid),g.target(eid)
print "add edge",eid2,g.source(eid2),g.target(eid2)
print "add edge",eid3,g.source(eid3),g.target(eid3)


print "has_vertex",g.has_vertex(sid),sid in g,1000 in g
print "has_edge",g.has_edge(eid),g.has_edge(1000)

print "in edges",list(g.in_edges(sid)),list(g.in_edges(tid))
print "nb in edges",g.nb_in_edges(sid),g.nb_in_edges(tid)
print "out edges",list(g.out_edges(sid)),list(g.out_edges(tid))
print "nb out edges",g.nb_out_edges(sid),g.nb_out_edges(tid)
print "edges",list(g.edges(sid)),list(g.edges(tid))
print "nb edges",g.nb_edges(sid),g.nb_edges(tid)

print "vertices",len(g),g.nb_vertices(),list(g.vertices())
print "in neighbors",list(g.in_neighbors(sid)),list(g.in_neighbors(tid))
print "nb in neighbors",g.nb_in_neighbors(sid),g.nb_in_neighbors(tid)
print "out neighbors",list(g.out_neighbors(sid)),list(g.out_neighbors(tid))
print "nb out neighbors",g.nb_out_neighbors(sid),g.nb_out_neighbors(tid)
print "neighbors",list(g.neighbors(sid)),list(g.neighbors(tid))
print "nb neighbors",g.nb_neighbors(sid),g.nb_neighbors(tid)

print "remove vertex",g.remove_vertex(tid)

print "edges",list(g.edges())
print "neighbors",list(g.neighbors(sid))

print "remove edge",g.remove_edge(eid2)
print "edges",list(g.edges())
print "neighbors",list(g.neighbors(sid))

print "clear",g.clear(),len(g),g.nb_edges()
sid=g.add_vertex()
tid=g.add_vertex(10)
tid2=g.add_vertex()
eid=g.add_edge(sid,tid)
eid2=g.add_edge(sid,tid2,34)
print "clear edges",g.clear_edges(),len(g),g.nb_edges()
print list(g.neighbors(sid))
