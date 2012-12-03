from openalea.container import PropertyMap,FactorMap

p=PropertyMap()
def pr () :
    print list(p.iteritems())

for i in xrange(5) :
    p[i]=float(i)
pr()
p+=1.
pr()
p-=1.
pr()
p*=2.
pr()
p/=2.
pr()

f=FactorMap(p,10.)
for i in f :
    print i,f[i]

print list(f.iteritems())
