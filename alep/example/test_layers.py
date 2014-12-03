from alinea.adel.data_samples import adel_two_metamers_stand
from alinea.alep.fungal_objects import Fungus
from alinea.alep.architecture import get_leaves
from alinea.alep.protocol import disperse
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport
from alinea.septo3d.dispersion.alep_interfaces import Septo3DTransport
from alinea.astk.Weather import sample_weather
from alinea.astk.TimeControl import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.ion()
import numpy
import pandas
from itertools import product


from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer


class DummyEmission():
    def __init__(self, domain):
        self.domain = domain
        
    def get_dispersal_units(self, g, fungus_name="dummy", label='LeafElement', weather_data=None):
        DU={}
        lesions = {k:[l for l in les if l.fungus.name is fungus_name] for k, les in g.property('lesions').iteritems()} 
        for vid, l in lesions.iteritems():
            for lesion in l:
                emissions = lesion.emission(1e7)
                try:
                    DU[vid] += emissions
                except:
                    DU[vid] = emissions
        return DU

        
def test_transport(interleaf=10., density = 350., layer_thickness=0.01):
    seq, weather = sample_weather()
    every_rain = rain_filter(seq, weather)  
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    evalvalue = rain_timing.next()
    
    g, domain_area, domain, convunit = adel_two_metamers_stand(density = density, 
                                        interleaf = interleaf, leaf_length = 20, leaf_width = 1, Einc = 0)
    leaves = get_leaves(g)
    lf = g.node(leaves[0])
    fungus = Fungus(name='dummy')
    lf.lesions = [fungus.lesion()]
    emitter = DummyEmission(domain=domain)
    # transporter = PopDropsTransport(domain = domain, domain_area = domain_area, dh = layer_thickness, convUnit = convunit)
    transporter = Septo3DTransport(domain = domain, domain_area = domain_area, dh = layer_thickness, convUnit = convunit, wash = False)
    g = disperse(g, emitter, transporter, fungus_name="dummy", label='LeafElement', weather_data=evalvalue.value)

    dus = g.property('dispersal_units')
    if 19 in dus:
        return len(g.node(19).dispersal_units)
    else:
        return 0.
    
def plot_test():
    interleaf = numpy.arange(1, 51, 1)
    layer_thickness = numpy.arange(0.1, 5.1, 0.1)
    df = pandas.DataFrame(index = interleaf, columns = layer_thickness)
    for int_lf in df.index:
        for dh in df.columns:
            print int_lf, dh
            df.ix[int_lf, dh] = test_transport(interleaf=float(int_lf), layer_thickness=dh*0.01)
            
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    X, Y = numpy.meshgrid(df.columns, df.index)
    ax.plot_surface(X,Y, df.values/max(df.max()), rstride=1, cstride=1, color='b', alpha=0.5)
    ax.set_xlabel('layer thickness (cm)')
    ax.set_ylabel('inter leaf (cm)')
    ax.set_zlabel('normalised number of deposits on leaf 2')
    ax.set_zlim([0,1])
    
# Faire des feuilles de 20 cm
# les écarter plus
# faire test pour dh jusqu'à 5
# faire test pour distance jusqu'à 50 cm