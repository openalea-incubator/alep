from alinea.simcycle.Lesions import *


def test_init_del():
    p = ParCycle()
    l = Lesions(p)
    del(p)
    del(l)

def test_field():
    p = ParCycle()
    d = p.to_dict()
    l = Lesions(p)
    d = l.to_dict()
    