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
    assert isinstance(les, lesion.fungus.Lesion)

    dem = disp.emission(spo)
    assert 'leaf' in dem
    n = dem['leaf'][0]
    assert n == 1

    dus = disp.get_dispersal_units(spo, dem)
    assert 'leaf' in dus
    du = dus['leaf'][0]
    assert isinstance(du, lesion.fungus.DispersalUnit)
    assert du.nb_dispersal_units == 1

    src = count_dus(dus)
    assert 'leaf' in src
    assert src['leaf'] == 1

    tmap = disp.transport(src)
    assert 'leaf' in tmap
    assert 'leaf' in tmap['leaf']
    assert tmap['leaf']['leaf'] == 1

    dep = disp.deposits(tmap, dus)
    assert 'leaf' in dep
    du = dep['leaf'][0]
    assert isinstance(du, lesion.fungus.DispersalUnit)
    assert du.nb_dispersal_units == 1

    deposits, loss = disp.disperse(lesions)
    assert loss == 0


def test_emission():
    fungus = Fungus()
    disp = Dispersal()

    lesions = {'leaf': [fungus.lesion() for _ in range(10)]}
    spo = disp.get_sporulating_lesions(lesions)
    dem = disp.emission(spo)
    dus = disp.get_dispersal_units(spo, dem)
    assert len(dus['leaf']) == 1
    assert dus['leaf'][0].nb_dispersal_units == 10

    lesions = {'leaf': [fungus.lesion(), fungus.lesion(mutation='muted')]}
    spo = disp.get_sporulating_lesions(lesions)
    dem = disp.emission(spo)
    dus = disp.get_dispersal_units(spo, dem)
    assert len(dus['leaf']) == 2
    for du in dus['leaf']:
        assert du.nb_dispersal_units == 1

def test_deposits():
    fungus = Fungus()
    disp = Dispersal()
    lesions = {'leaf_1': [fungus.lesion()],
               'leaf_2': [fungus.lesion()]}
    spo = disp.get_sporulating_lesions(lesions)
    dem = disp.emission(spo)
    dus = disp.get_dispersal_units(spo, dem)
    transport = {'leaf_1': {'leaf_1' : 1, 'leaf_2':1}}
    dep = disp.deposits(transport, dus)
    assert len(dep) == 1
    assert len(dep['leaf_1']) == 1
    assert dep['leaf_1'][0].nb_dispersal_units == 2

    