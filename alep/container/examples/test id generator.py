from openalea.container import IdGenerator

g=IdGenerator()
print g.get_id(5)
for i in xrange(10) :
    print g.get_id()

print "none",g.get_id(None)

g.release_id(5)
for i in xrange(3) :
    print g.get_id()

g.clear()
for i in xrange(3) :
    print g.get_id(i)

for i in xrange(3) :
    print g.get_id()
