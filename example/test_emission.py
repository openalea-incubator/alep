import numpy as np
import pandas as pd
from alinea.alep.wheat import adel_one_leaf
from alinea.alep.protocol import disperse
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport
from alinea.alep.fungal_objects import DispersalUnit, Lesion
from alinea.alep.septoria import plugin_septoria

class DummyLesion(Lesion):
    def __init__(self, nb_spores=None, position=None):
        super(DummyLesion, self).__init__(nb_spores=nb_spores, position=position)
        class params():
            def __init__(self, name="dummy"):
                self.name = name
                self.spo_rate = 10000.
        self.fungus = params()
        self.surface_spo = 10000.
        self.surface_empty = 0.
        self.hist_emission = [0.]
        self.hist_empty = [0.]
        self.stock_spores = self.fungus.spo_rate*self.surface_spo
        self.hist_stock = [self.stock_spores]
        self.hist_spo = [self.surface_spo]
        
    def is_sporulating(self):
        return True
    
    def update_empty_surface(self, nb_spores_emitted):
        new_surface_empty = (nb_spores_emitted/self.stock_spores) * self.surface_spo
        self.surface_empty += new_surface_empty
        self.surface_spo = max(0., self.surface_spo-new_surface_empty)
        self.hist_emission.append(nb_spores_emitted)
        self.hist_empty.append(self.surface_empty)
        self.hist_spo.append(self.surface_spo)
        
    def reduce_stock(self, nb_spores_emitted):
        self.stock_spores -= nb_spores_emitted
        self.hist_stock.append(self.stock_spores)
        
class DummyEmission():        
    def get_dispersal_units(self, g, fungus_name="dummy", label='LeafElement'):
        DU={}
        lesions = {k:[l for l in les if l.fungus.name is fungus_name and l.is_sporulating()] 
                    for k, les in g.property('lesions').iteritems()} 
        density_dus = 2000.
        for vid, l in lesions.iteritems():
            for lesion in l:
                emissions = lesion.emission(density_dus, rain_exposition=0.5)
                print 'Emissions: %f' % emissions
                try:
                    DU[vid] += []
                except:
                    DU[vid] = []
        return DU

class DummyTransport():
    def disperse(self, g, dispersal_units, time_control = None):
        return {}
        
g = adel_one_leaf()
septoria = plugin_septoria()
lesion = septoria.lesion()
lesion.surfaces_spo = np.array([1., 0., 0.])
lesion.status = lesion.fungus.SPORULATING
lesion.status_edge = lesion.fungus.SPORULATING
print lesion.is_sporulating()
g.node(10).lesions = [lesion]
emitter = PopDropsEmission()
transporter = PopDropsTransport()
nb_steps = 25
rain = {'rain':[1.,1.,1.,1.]}
weather_data = pd.DataFrame(rain)
# for i in range(4):
for i in range(nb_steps):
    les = g.property('lesions')
    lesion = les[10][0]
    print '-----------------------------------'
    print 'Dispersion %d' % i
    print 'Surfaces spo: ' +  str(lesion.surfaces_spo)
    print 'Surface spo total: %f' % lesion.surface_spo
    print 'Surface empty: %f' % lesion.surface_empty
    # if i<3:
    if i<nb_steps-1:
        disperse(g, emitter, transporter, "septoria", label='LeafElement', weather_data=weather_data)