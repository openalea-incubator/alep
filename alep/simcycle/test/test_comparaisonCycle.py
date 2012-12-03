from alinea.simcycle.Lesions import *

from alinea.simcycle import cycle
from alinea.simcycle.cycle import septoria



    #nh = 24
    #T,PPFD, Rh = [20] * nh,[0]*24, [100]*24

environment = {}
environment['relative_humidity'] = 100
environment['PPFD'] = 0
    
class StubLeaf(object):
    def __init__(self):
        self.wetness = True
        self.temp = 20.
        self.healthy_surface = 10.
        self.rain_intensity = 0.
        self.lesions = []

 

#init lesion cycle
leaf = StubLeaf()
septoria = septoria()
septoria.loss_rate = 0
s = cycle.Lesion(fungus = septoria)
leaf.lesions.append(s)

#init lesion simcycle avec 10 lesions
leafsc = StubLeaf()
leafsc.healthy_surface *= 10
septo = ParCycle()
#force identical parameter values for equivalent parameters
septo.txPerteUdin = 0
septo.Slmax = septoria.Smax
septo.Slmin = septoria.Smin
septo.TbasePath = septoria.basis_for_dday
septo.DInc =  septoria.degree_days_to_chlorosis
septo.DChlo = septoria.degree_days_to_sporulation
septo.SIncMin = septoria.epsilon
septo.txCroi =  septoria.growth_rate
septo.DminGerm =  septoria.wd_min
lsc = Lesions(septo)
lsc.addUdin(N=10,Q=10)
leafsc.lesions.append(lsc)

#simul

steps = 240

for i in range(steps):
    for s in leaf.lesions:
        sbef = s.surface
        s.update(1, 1./ 24 * (leaf.temp - septoria.basis_for_dday), leaf, environment)
        leaf.healthy_surface -= (s.surface - sbef)
    #for s in leafsc.lesions:
    leafsc.healthy_surface = lsc.devLes(leafsc.healthy_surface, 1, [leafsc.temp], [environment['PPFD']], [environment['relative_humidity']])