""" test template fungus class"""
from alinea.alep.fungus import Fungus, Lesion, DispersalUnit

def test_instantiation():
    fungus = Fungus()
    assert hasattr(fungus, 'length_unit')
    assert hasattr(fungus, 'Lesion_class')
    assert fungus.Lesion_class is Lesion
    assert hasattr(fungus, 'DispersalUnit_class')
    assert fungus.DispersalUnit_class is DispersalUnit
    assert hasattr(fungus, 'parameter_names')
    attr = vars(fungus)
    assert all([name in attr for name in fungus.parameter_names])
    
    assert hasattr(fungus, 'lesion') and callable(fungus.lesion)
    assert hasattr(fungus, 'dispersal_unit') and callable(fungus.dispersal_unit)
    
def test_lesion_du_instantiation():
    """Test whether fungus generate appropriate lesion and dispersal units"""
    fungus = Fungus()
    lesion = fungus.lesion()
    assert lesion.fungus.name == fungus.name
    du = fungus.dispersal_unit()
    assert du.fungus.name == fungus.name
    
    #todo : test mutation of parameters when fungus -> lesion -> du -> lesion
    
def test_lesion_emission():
    fungus = Fungus()
    lesion = fungus.lesion()
    du = lesion.emission(nb_DU=3)
    assert isinstance(du, lesion.fungus.DispersalUnit_class)
    assert du.nb_dispersal_units == 3
    