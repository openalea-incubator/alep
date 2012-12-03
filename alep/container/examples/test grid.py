from openalea.container import Grid

g=Grid( (2,3) )
print "dim",g.dim()
print "shape",g.shape()
print "size",g.size(),len(g)
for i in g :
    print "iter",i,"coord",g.coordinates(i),"index",g.index(g.coordinates(i))

g=Grid( (2,3,4) )
print "dim",g.dim()
print "shape",g.shape()
print "size",g.size(),len(g)
for i in g :
    print "iter",i,"coord",g.coordinates(i),"index",g.index(g.coordinates(i))
