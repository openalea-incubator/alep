from openalea.container import IdMapFloat

m=IdMapFloat()
print m.add(1.,5)
for i in xrange(10) :
    print m.add(i*3.)

print "none",m.add(4.,None)

print "len",len(m)
print "contains",1 in m
for i in m.iteritems() :
    print "items",i
for i in m :
    print "keys",i,m[i]

m[1]=3.14
print "setitem",m[1]
m[100]=3.14
print "setitem",m[100]
print "pop",m.pop(100)
print "popitem",m.popitem()
print "setdefault",m.setdefault(1,7.64)
print "setdefault",m.setdefault(100,7.64)

print "clear",m.clear()
print "clear",list(m.iteritems())

print "copy",dict(m.copy())
