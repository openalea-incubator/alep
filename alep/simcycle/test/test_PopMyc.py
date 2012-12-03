from alinea.simcycle.PopMyc import *


def test_init_del():
    c = ClassMyc()
    del(c)
    s = StackCmyc()
    del(s)
    s = StackBpop()
    del(s)
    p = PopMyc()
    del(p)
    
    
def test_field_PopMyc():
     
    p = PopMyc(agemin =0, agemax = 100, largClass = 10, popType = 'Les')
    assert p.type == 0
    assert p.nbClass == 10
    assert p. agemin == 0
    assert p.agemax == 100
    assert p.largClass == 10
    assert p.pC0 == 0
    assert p.TimeBuffer == 0
    q = p.getQ()
    assert sum(q) == 0
    n = p.getN()
    assert sum(n) == 0
    d = p.getdecT()
    assert sum(d) == 0
    p.to_dict() # works only if fields are correctly aligned
    
#def test_repr():
#    p = PopMyc()
#    print(p)