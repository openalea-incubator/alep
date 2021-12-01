import matplotlib.pyplot as plt

from alinea.alep.diseases import get_disease


class StubLeaf(object):
    def __init__(self):
        self.wetness = True
        self.temp = 18.
        self.age = 1.
        self.healthy_surface = 10.
        self.rain_intensity = 0.
        self.relative_humidity = 85.
        self.senesced_length=0
        self.lesions = []

    def properties(self):
        return {'temperature_sequence': [self.temp] * 24, 'relative_humidity_sequence': [self.relative_humidity]*24 }
        

def setup():
    leaf = StubLeaf()
    dt = 1.
    nb_steps = 1000
    return leaf, dt, nb_steps


def test_powdery():
    leaf, dt, nb_steps = setup()
    fungus = get_disease('powdery_mildew')
    s = fungus.lesion()
    leaf.lesions.append(s)

    for i in range(nb_steps):
        if i%6 == 0:
            leaf.rain_intensity = 5
        else:
            leaf.rain_intensity = 0
        s.update(dt, leaf)

    return s, leaf

def test_septoria():
    leaf, dt, nb_steps = setup()

    fungus = get_disease('septoria')
    s = fungus.lesion()
    leaf.lesions.append(s)

    for i in range(nb_steps):
        s.update(dt, leaf)
        #s.control_growth(s.growth_demand)

    # Emptying of the lesions
    leaf.rain_intensity = 5
    s.update(dt, leaf)
    leaf.rain_intensity = 0
    s.update(dt, leaf)
    leaf.rain_intensity = 5
    s.update(dt, leaf)
    leaf.rain_intensity = 0
    s.update(dt, leaf)
    leaf.rain_intensity = 5
    s.update(dt, leaf)

    return s

def test_septoria_rain():
    leaf, dt, nb_steps = setup()

    fungus = get_disease('septoria')
    s = fungus.lesion()
    leaf.lesions.append(s)

    for i in range(nb_steps):
        if i%6 == 0:
            leaf.rain_intensity = 5
        else:
            leaf.rain_intensity = 0
        s.update(dt, leaf)

    return s

def test_death_by_drought():
    leaf, dt, nb_steps = setup()

    fungus = get_disease('septoria')
    s = fungus.lesion()
    leaf.lesions.append(s)

    for i in range(nb_steps):
        if i == 4:
            leaf.wetness = False
        else:
            leaf.wetness = True
        s.update(dt, leaf)
    if len(leaf.lesions) > 0:
        assert s.status == 6, s.status

def test_status_changes():
    leaf, dt, nb_steps = setup()

    fungus = get_disease('septoria')
    s = fungus.lesion()
    leaf.lesions.append(s)

    for i in range(nb_steps):
        s.update(dt, leaf)
        # Test of the emergence of the first ring after 10 h of continuous wetness duration
        if i == 10 :
            assert s.status == 1, s.status

        # Test of the begin of incubation 220DD after the emergence of the lesion
        if s.age_tt >= s.fungus.degree_days_to_chlorosis:
            assert s.status > 2, s.status

        # Test of the begin of incubation 330DD after the emergence of the lesion
        if s.age_tt >= s.fungus.degree_days_to_chlorosis + s.fungus.degree_days_to_sporulation:
            assert s.status > 3, s.status

        # Test of the depletion of the lesion after 3 rain events
        if s.status == 4:
            if i%2 == 0:
                leaf.rain_intensity = 5
            else:
                leaf.rain_intensity = 0

        assert s.status > 4, s.status

def test_Smax():
    leaf, dt, nb_steps = setup()

    fungus = get_disease('septoria')
    s = fungus.lesion()
    leaf.lesions.append(s)

    for i in range(nb_steps):
        s.update(dt, leaf)

    assert s.surface <= s.fungus.Smax

def test_septoria_hist():
    leaf, dt, nb_steps = setup()

    fungus = get_disease('septoria')
    s = fungus.lesion()
    leaf.lesions.append(s)

    hist_lesion_age = []
    hist_lesion_age_dday = []
    hist_lesion_status = []
    hist_lesion_surface = []
    hist_ring_status = []
    hist_all_rings_status = []

    for i in range(nb_steps):
        s.update(dt, leaf)
        hist_lesion_age.append(s.age_physio)
        hist_lesion_age_dday.append(s.age_tt)
        hist_lesion_status.append(s.status)
        hist_lesion_surface.append(s.surface)
        # [hist_all_rings_status.append(hist_ring_status.append(r.status)) for r in s.rings]

    return (hist_lesion_age, hist_lesion_age_dday,
            hist_lesion_status, hist_lesion_surface, hist_all_rings_status)
            
