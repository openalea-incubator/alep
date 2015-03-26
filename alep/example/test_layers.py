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
    def __init__(self, domain, to_emit = 1e7):
        self.domain = domain
        self.to_emit = to_emit
        
    def get_dispersal_units(self, g, fungus_name="dummy", label='LeafElement', weather_data=None):
        DU={}
        lesions = {k:[l for l in les if l.fungus.name is fungus_name] for k, les in g.property('lesions').iteritems()} 
        for vid, l in lesions.iteritems():
            for lesion in l:
                emissions = lesion.emission(self.to_emit)
                try:
                    DU[vid] += emissions
                except:
                    DU[vid] = emissions
        return DU

        
def test_transport(interleaf=10., density = 350., layer_thickness=0.01, leaf_sectors = 1, by_sector = False):
    seq, weather = sample_weather()
    every_rain = rain_filter(seq, weather)  
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    evalvalue = rain_timing.next()
    
    g, domain_area, domain, convunit = adel_two_metamers_stand(leaf_sectors = leaf_sectors, density = density, 
                                        interleaf = interleaf, leaf_length = 20, leaf_width = 1, Einc = 0)
                                        
    #domain_area *= 10
                                        
    labels = g.property('label')
    bids = [v for v,l in labels.iteritems() if l.startswith('blade')]
    leaves = [[vid for vid in g.components(blade) if labels[vid].startswith('LeafElement')] for blade in bids]
                                        
    fungus = Fungus(name='dummy')
    for lf in leaves[0]:
        g.node(lf).lesions = [fungus.lesion()]
    emitter = DummyEmission(domain=domain, to_emit = 1e7 / float(leaf_sectors))
    transporter = PopDropsTransport(domain = domain, domain_area = domain_area, dh = layer_thickness, convUnit = convunit)
    # transporter = Septo3DTransport(domain = domain, domain_area = domain_area, dh = layer_thickness, convUnit = convunit, wash = False, show = False)
    g = disperse(g, emitter, transporter, fungus_name="dummy", label='LeafElement', weather_data=evalvalue.value)

    dus = g.property('dispersal_units')
    if by_sector == False:
        dus_on_source = sum([dus[lf] for lf in leaves[0] if lf in dus], []) 
        dus_on_target = sum([dus[lf] for lf in leaves[1] if lf in dus], [])
        return len(dus_on_source), len(dus_on_target)
    else:
        nb_dus_on_source = {lf:len(dus[lf]) if lf in dus else 0. for lf in leaves[0]}
        nb_dus_on_target = {lf:len(dus[lf]) if lf in dus else 0. for lf in leaves[1]}
        return nb_dus_on_source, nb_dus_on_target

def plot_interleaf(density = 350., layer_thickness=0.01, leaf_sectors = 1, 
                    linestyle = '-', ax = None, return_ax = True):
    dus_s = []
    dus_t = []
    nb_steps = 30
    for i in range(nb_steps):
        dus_s.append(test_transport(interleaf = i, density = density, 
                                    layer_thickness = layer_thickness, 
                                    leaf_sectors = leaf_sectors)[0])
        dus_t.append(test_transport(interleaf = i, density = density, 
                                    layer_thickness = layer_thickness, 
                                    leaf_sectors = leaf_sectors)[1])
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(range(nb_steps), dus_s, 'b', linestyle = linestyle)
    ax.plot(range(nb_steps), dus_t, 'r', linestyle = linestyle)
    ax.legend(['source leaf', 'target leaf'], loc = 'best')
    ax.set_ylabel('Number of DUs intercepted', fontsize = 14)
    ax.set_ylabel('Distance between leaves (cm)', fontsize = 14)
    if return_ax == True:
        return ax
        
def plot_layer_thickness_x_interleaf():
    interleaf = numpy.arange(1, 51, 1)
    layer_thickness = numpy.arange(0.1, 5.1, 0.1)
    df = pandas.DataFrame(index = interleaf, columns = layer_thickness)
    for int_lf in df.index:
        for dh in df.columns:
            print int_lf, dh
            _, df.ix[int_lf, dh] = test_transport(interleaf=float(int_lf), layer_thickness=dh*0.01)
            
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    X, Y = numpy.meshgrid(df.columns, df.index)
    ax.plot_surface(X,Y, df.values/max(df.max()), rstride=1, cstride=1, color='b', alpha=0.5)
    ax.set_xlabel('layer thickness (cm)')
    ax.set_ylabel('inter leaf (cm)')
    ax.set_zlabel('normalised number of deposits on leaf 2')
    ax.set_zlim([0,1])
    
def plot_sectors(by_sector = False):
    if by_sector == False:
        df = pandas.DataFrame(columns = ['source', 'target'])
        for nb_sect in range(1,11):
            print nb_sect
            df.loc[nb_sect, :] = test_transport(interleaf = 5., layer_thickness=0.01, leaf_sectors = nb_sect, by_sector = False)

        fig, ax = plt.subplots(1,1)
        df.plot(ax = ax)
        ax.set_xlabel('number of sectors')
        ax.set_ylabel('number of deposits on leaf')
        
    else:
        sectors = range(1, 11)
        df_source = pandas.DataFrame(index = sectors, columns = sectors)
        df_target = pandas.DataFrame(index = sectors, columns = sectors)
        for nb_sect in sectors:
            print nb_sect
            nb_dus_on_source, nb_dus_on_target = test_transport(interleaf = 5., layer_thickness=0.01, leaf_sectors = nb_sect, by_sector = True)
            df_source.loc[nb_sect, :nb_sect] = nb_dus_on_source.values()
            df_target.loc[nb_sect, :nb_sect] = nb_dus_on_target.values()

        fig, ax = plt.subplots(1,2)
        df_source.plot(ax = ax[0], marker='*')
        ax[0].set_xlabel('number of sectors', fontsize = 18)
        ax[0].set_ylabel('number of deposits on leaf', fontsize = 18)
        ax[0].set_title('Source leaf', fontsize = 18)
        
        df_target.plot(ax = ax[1], marker='*')
        ax[1].set_xlabel('number of sectors', fontsize = 18)
        ax[1].set_ylabel('number of deposits on leaf', fontsize = 18)
        ax[1].set_title('Target leaf', fontsize = 18)
    
    
# Faire des feuilles de 20 cm
# les écarter plus
# faire test pour dh jusqu'à 5
# faire test pour distance jusqu'à 50 cm