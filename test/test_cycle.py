from nose import with_setup
import matplotlib.pyplot as plt

from alinea.alep import septoria
from alinea.alep import powdery_mildew

class StubLeaf(object):
    def __init__(self):
        self.wetness = True
        self.temp = 18.
        self.age = 1.
        self.healthy_surface = 10.
        self.rain_intensity = 0.
        self.relative_humidity = 85.
        self.lesions = []
        
class TestClass:
    def setUp(self):
        self.leaf = StubLeaf()
        self.dt = 1.
        self.nb_steps = 1000

    def tearDown(self):
        print("TearDown")

    def test_powdery(self):
        leaf = self.leaf
        
        s = cycle.Lesion(fungus = powdery_mildew())
        leaf.lesions.append(s)
        
        for i in range(self.nb_steps):
            if i%6 == 0:
                leaf.rain_intensity = 5
            else:
                leaf.rain_intensity = 0
            s.update(self.dt, leaf)
        
        return s, leaf

    def test_septoria(self):
        dt = self.dt
        leaf = self.leaf
        
        s = cycle.Lesion(fungus = septoria())
        leaf.lesions.append(s)
       
        for i in range(self.nb_steps):
            s.update(self.dt, leaf)
        
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
        
    def test_septoria_rain(self):
        leaf = self.leaf
        
        s = cycle.Lesion(fungus = septoria())
        leaf.lesions.append(s)
        
        for i in range(self.nb_steps):
            if i%6 == 0:
                leaf.rain_intensity = 5
            else:
                leaf.rain_intensity = 0
            s.update(self.dt, leaf)
        
        return s
        
    def test_death_by_drought(self):
        leaf = self.leaf
        
        s = cycle.Lesion(fungus = septoria())
        leaf.lesions.append(s)
        
        for i in range(self.nb_steps):
            if i == 4:
                leaf.wetness = False
            else:
                leaf.wetness = True
            s.update(self.dt, leaf)
        if len(leaf.lesions) > 0:
            assert s.status == 6, s.status
    
    def test_status_changes(self):
        leaf = self.leaf

        s = cycle.Lesion(fungus = septoria())
        leaf.lesions.append(s)
        
        for i in range(self.nb_steps):
            s.update(self.dt, leaf)
            # Test of the emergence of the first ring after 10 h of continuous wetness duration
            if i == 10 :
                assert s.status == 1, s.status
                
            # Test of the begin of incubation 220DD after the emergence of the lesion
            if s.age_dday >= s.fungus.degree_days_to_chlorosis:
                assert s.status > 2, s.status
                
            # Test of the begin of incubation 330DD after the emergence of the lesion
            if s.age_dday >= s.fungus.degree_days_to_chlorosis + s.fungus.degree_days_to_sporulation:
                assert s.status > 3, s.status
            
            # Test of the depletion of the lesion after 3 rain events
            if s.status == 4:
                if i%2 == 0:
                    leaf.rain_intensity = 5
                else:
                    leaf.rain_intensity = 0
            if s.rings[0].cumul_rain_event >= 3:
                assert s.status > 4, s.status

    def test_Smax(self):
        leaf = self.leaf
        
        s = cycle.Lesion(fungus = septoria())
        leaf.lesions.append(s)
        
        for i in range(self.nb_steps):
            s.update(self.dt, leaf)
        
        assert sum([r.surface for r in s.rings]) <= s.fungus.Smax

    def test_septoria_hist(self):
        leaf = self.leaf

        s = cycle.Lesion(fungus = septoria())
        leaf.lesions.append(s)
        
        hist_lesion_age = []
        hist_lesion_age_dday = []
        hist_lesion_status = []
        hist_lesion_surface = []
        hist_ring_status = []
        hist_all_rings_status = []
        
        for i in range(self.nb_steps):
            s.update(self.dt, leaf)
            hist_lesion_age.append(s.age)
            hist_lesion_age_dday.append(s.age_dday)
            hist_lesion_status.append(s.status)
            hist_lesion_surface.append(sum([r.surface for r in s.rings]))
            # [hist_all_rings_status.append(hist_ring_status.append(r.status)) for r in s.rings]
            
        return (hist_lesion_age, hist_lesion_age_dday, 
                hist_lesion_status, hist_lesion_surface, hist_all_rings_status)
            
    
# create a du
# simulate its growth
# check the number of spores, rings, age, ...

if __name__ == '__main__':
    test = TestClass()
    test.setUp()
    s = test.test_septoria()
    print('TEST 1 : Run Septoria')
    print(' ------------------------------------------------------------------------------- ')
    print('1st ring status must be 5 = EMPTY. It is equal to ', s.rings[0].status, end=' ')  
    
    test.setUp()
    s = test.test_septoria_rain()
    print(2*'\n')
    print('TEST 2 : Run Septoria with rain for infection')
    print(' ------------------------------------------------------------------------------- ')
    print('1st ring status must be 5 = EMPTY or 6 = DEAD. It is equal to', s.rings[0].status)
    
    test.setUp()
    test.test_death_by_drought()
    
    test.setUp()
    test.test_status_changes()
    
    test.setUp()
    test.test_Smax()
    
    test.setUp()
    age, age_dday, status, surface, ring_status = test.test_septoria_hist()
    print(5*'\n')
    plt.plot(age_dday, surface)
    plt.xlabel('Lesion thermal age in degree-days')
    plt.ylabel('Lesion surface in cm2')
    plt.show()