from openalea.container import *

r=Relation()
left=[r.add_left_element() for i in xrange(10)]
right=[r.add_right_element() for i in xrange(5)]
for i in xrange(5) :
    r.add_link(left[i],right[i])
    r.add_link(left[5+i],right[i])

print left
print list(r.left_elements())
print right
print list(r.right_elements())
print list(r.links())
print [(r.left(lid),r.right(lid)) for lid in r.links()]
