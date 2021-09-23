from alinea.alep.dispersal import Emission
from alinea.alep.fungus import Fungus


def test_emission():
    fungus = Fungus()
    lesion = fungus.lesion()
    lesions = {'leaf': [lesion]}
    emission_model = Emission()
    
    assert not lesion.is_sporulating
    dus = emission_model.get_dispersal_units(lesions)
    assert len(dus) == 0
    
    lesion.is_sporulating = True
    emissions = emission_model.get_dispersal_units(lesions)
    assert 'leaf' in emissions
    dus = emissions['leaf']
    assert len(dus) == 1
    du = dus[0]
    assert isinstance(du, lesion.fungus.DispersalUnit_class)
    assert du.nb_dispersal_units == 1
    
    