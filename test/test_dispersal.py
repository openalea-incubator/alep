from alinea.alep.dispersal import Dispersal, count_dus
from alinea.alep.fungus import Fungus


def test_generic_methods():
    fungus = Fungus()
    lesion = fungus.lesion()
    assert lesion.is_sporulating()

    lesions = {'leaf': [lesion]}
    disp = Dispersal()

    spo = disp.get_sporulating_lesions(lesions)
    assert 'leaf' in spo
    les = spo['leaf'][0]
    assert isinstance(les, lesion.fungus.Lesion_class)

    dem = disp.emission_demands(spo)
    assert 'leaf' in dem
    n = dem['leaf'][0]
    assert n == 1

    dus = disp.get_dispersal_units(spo, dem)
    assert 'leaf' in dus
    du = dus['leaf'][0]
    assert isinstance(du, lesion.fungus.DispersalUnit_class)
    assert du.nb_dispersal_units == 1

    src = count_dus(dus)
    assert 'leaf' in src
    assert src['leaf'] == 1

    tmap = disp.transport_map(src)
    assert 'leaf' in tmap
    assert 'leaf' in tmap['leaf']
    assert tmap['leaf']['leaf'] == 1

    dep = disp.deposits(tmap, dus)
    assert 'leaf' in dep
    du = dep['leaf'][0]
    assert isinstance(du, lesion.fungus.DispersalUnit_class)
    assert du.nb_dispersal_units == 1

    deposits, loss = disp.disperse(lesions)
    assert loss == 0

