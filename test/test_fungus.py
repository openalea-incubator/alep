""" test template fungus class"""
from alinea.alep.fungus import Fungus, Lesion, DispersalUnit, DotDict


def test_base_api():
    fungus = Fungus(a_parameter='dummy')
    assert hasattr(fungus, 'length_unit')
    assert hasattr(fungus, 'name')
    assert hasattr(fungus, 'parameters')
    assert isinstance(fungus.parameters, DotDict)
    assert 'a_parameter' in fungus.parameters
    assert fungus.parameters.a_parameter == 'dummy'
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
    assert lesion.fungus.name == fungus.name
    du = fungus.dispersal_unit()
    assert du.fungus.name == fungus.name


def test_parameters():
    fungus = Fungus(fixed='origin', mutable='origin')
    base_lesion = fungus.lesion()
    lesion = fungus.lesion(mutable='muted')

    fp = fungus.parameters
    p = lesion.parameters()
    bp = base_lesion.parameters()
    assert p.fixed == fp.fixed == bp.fixed == 'origin'
    assert p.mutable == 'muted'
    assert fp.mutable == bp.mutable == 'origin'

    du = fungus.dispersal_unit(mutable='du_muted')
    fp = fungus.parameters
    p = lesion.parameters()
    dup = du.parameters()
    assert p.fixed == fp.fixed == dup.fixed == 'origin'
    assert p.mutable == 'muted'
    assert dup.mutable == 'du_muted'
    assert fp.mutable == 'origin'

    
def test__emission():
    fungus = Fungus(latency=100)
    lesion = fungus.lesion(latency=200)

    du = lesion.emission(emission_demand=3, latency=300)
    assert isinstance(du, lesion.fungus.DispersalUnit)
    assert du.nb_dispersal_units == 3
    p = du.parameters()
    assert p.latency == 300
    pp = lesion.parameters()
    assert pp.latency == 200
    ppp = fungus.parameters
    assert ppp.latency == 100



    