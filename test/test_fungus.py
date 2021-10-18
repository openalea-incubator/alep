""" test template fungus class"""
from alinea.alep.fungus import Fungus, Lesion, DispersalUnit, DotDict


def test_base_api():
    fungus = Fungus(a_parameter='dummy')
    assert hasattr(fungus, '_parameters')
    fp = fungus.parameters()
    assert isinstance(fp, DotDict)
    assert 'name' in fp
    assert 'a_parameter' in fp
    assert fp.a_parameter == 'dummy'
    assert hasattr(fungus, 'Lesion')
    assert fungus.Lesion is Lesion
    assert hasattr(fungus, 'DispersalUnit')
    assert fungus.DispersalUnit is DispersalUnit
    assert hasattr(fungus, 'lesion') and callable(fungus.lesion)
    assert hasattr(fungus, 'dispersal_unit') and callable(fungus.dispersal_unit)


def test_fungus_transmission():
    """Test whether fungus generate appropriate lesion and dispersal units"""
    fungus = Fungus()
    lesion = fungus.lesion()
    assert type(lesion.fungus) is type(fungus)
    assert lesion.fungus._parameters ==  fungus._parameters

    du = fungus.dispersal_unit()
    assert type(du.fungus) is type(fungus)
    assert du.fungus._parameters ==  fungus._parameters

    my_fungus = Fungus(name='my_fungus')
    my_lesion = my_fungus.lesion()
    assert type(my_lesion.fungus) is type(my_fungus)
    assert my_lesion.fungus._parameters == my_fungus._parameters
    my_du = my_fungus.dispersal_unit()
    assert type(my_du.fungus) is type(my_fungus)
    assert my_du.fungus._parameters == my_fungus._parameters
    # check fungus has not been affected
    assert type(lesion.fungus) is type(fungus)
    assert lesion.fungus._parameters ==  fungus._parameters
    assert type(du.fungus) is type(fungus)
    assert du.fungus._parameters ==  fungus._parameters


def test_parameters():
    fungus = Fungus(fixed='origin', mutable='origin')
    base_lesion = fungus.lesion()
    lesion = fungus.lesion(mutable='muted')

    fp = fungus.parameters()
    p = lesion.parameters()
    bp = base_lesion.parameters()
    assert p.fixed == fp.fixed == bp.fixed == 'origin'
    assert p.mutable == 'muted'
    assert fp.mutable == bp.mutable == 'origin'

    du = fungus.dispersal_unit(mutable='du_muted')
    fp = fungus.parameters()
    p = lesion.parameters()
    dup = du.parameters()
    assert p.fixed == fp.fixed == dup.fixed == 'origin'
    assert p.mutable == 'muted'
    assert dup.mutable == 'du_muted'
    assert fp.mutable == 'origin'

    
def test_emission():
    fungus = Fungus(latency=100)
    lesion = fungus.lesion(latency=200)

    du = lesion.emission(emission_demand=3, latency=300)
    assert isinstance(du, lesion.fungus.DispersalUnit)
    assert du.nb_dispersal_units == 3
    p = du.parameters()
    assert p.latency == 300
    pp = lesion.parameters()
    assert pp.latency == 200
    ppp = fungus.parameters()
    assert ppp.latency == 100

def test_is_like():
    fungus = Fungus(latency=100)
    les = fungus.lesion()
    du = les.emission()

    les2 = fungus.lesion(latency=100)
    du2 = les2.emission(emission_demand=10)
    assert les.is_like(les2)
    assert du.is_like(du2)

    les2 = fungus.lesion(latency=200)
    du2 = les2.emission()
    assert not les.is_like(les2)
    assert not du.is_like(du2)

    fungus = Fungus()
    les = fungus.lesion()
    fungus2 = Fungus(name='my_fungus')
    les2 = fungus2.lesion()
    assert not les.is_like(les2)

    class MyLesion(Lesion): pass
    fungus = Fungus()
    les = fungus.lesion()
    my_fungus = Fungus(lesion=MyLesion)
    les2 = my_fungus.lesion()
    assert not les.is_like(les2)


    