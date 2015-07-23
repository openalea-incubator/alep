from alinea.adel.data_samples import adel_two_metamers_stand
from alinea.alep.fungal_objects import Fungus
from alinea.alep.protocol import disperse
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport
from alinea.septo3d.dispersion.alep_interfaces import Septo3DTransport
from alinea.astk.Weather import sample_weather
from alinea.astk.TimeControl import *
import matplotlib.pyplot as plt
plt.ion()
import numpy
import pandas
from math import ceil

from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer

# Examples with two metamers ##################################################
class DummyEmission():
    def __init__(self, domain, to_emit = 1e7):
        self.domain = domain
        self.to_emit = to_emit
        
    def get_dispersal_units(self, g, fungus_name="dummy", label='LeafElement', weather_data=None):
        DU={}
        for vid, l in lesions.iteritems():
            for lesion in l:
                emissions = lesion.emission(self.to_emit)
                import pdb
                pdb.set_trace()
                try:
                    DU[vid] += emissions
                except:
                    DU[vid] = emissions
        return DU

def sample_weather_with_rain():
    seq, weather = sample_weather()
    every_rain = rain_filter(seq, weather)  
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    return rain_timing.next().value

def get_leaf_ids(g):
    labels = g.property('label')
    bids = [v for v,l in labels.iteritems() if l.startswith('blade')]
    return [[vid for vid in g.components(blade) if labels[vid].startswith('LeafElement')] for blade in bids]

def emission_source(g, leaves, domain_area, density_emitted=1e5):
    DU = {}    
    source = leaves[0]
    lf = source[int(ceil(len(source)/2))]
    DU[lf] = density_emitted*domain_area
    return DU
    
def count_du_source_target(deposits, leaves, by_sector=False):
    if by_sector == False:
        nb_dus_on_source = sum([deposits[lf][0].nb_dispersal_units for lf in leaves[0] if lf in deposits and len(deposits[lf])>0]) 
        nb_dus_on_target = sum([deposits[lf][0].nb_dispersal_units for lf in leaves[1] if lf in deposits and len(deposits[lf])>0])
    else:
        nb_dus_on_source = {lf:deposits[lf][0].nb_dispersal_units if lf in dus else 0. 
                            for lf in leaves[0]}
        nb_dus_on_target = {lf:deposits[lf][0].nb_dispersal_units if lf in dus else 0. 
                            for lf in leaves[1]}
    return nb_dus_on_source, nb_dus_on_target
    
def test_transport_rain_two_metamers(interleaf=10.,
                                     density = 350.,
                                     layer_thickness=0.01, 
                                     leaf_sectors = 2,
                                     density_emitted=1e5,
                                     by_sector = False):
    g, domain_area, domain, convunit = adel_two_metamers_stand(leaf_sectors=leaf_sectors, 
                                                               density=density,
                                                               interleaf=interleaf,
                                                               leaf_length=20,
                                                               leaf_width=1,
                                                               Einc=0)
    weather = sample_weather_with_rain()
    leaves = get_leaf_ids(g)
    DU = emission_source(g, leaves, 
                         domain_area=domain_area,
                         density_emitted=density_emitted)
    transporter = PopDropsTransport(domain=domain, domain_area=domain_area,
                                    dh=layer_thickness, convUnit=convunit,
                                    group_dus=True)
    if sum(DU.values())>0:
        deposits = transporter.disperse(g, DU, weather)
    return count_du_source_target(deposits, leaves, by_sector=False)

def plot_interleaf_two_metamers(density = 350., layer_thickness=0.01,
                                leaf_sectors = 1, linestyle = '-', 
                                ax = None, return_ax = True):
    dus_s = []
    dus_t = []
    interleaves = numpy.arange(0,31,5)
    for i in interleaves:
        s, t = test_transport_rain_two_metamers(interleaf = i, density = density, 
                                                layer_thickness = layer_thickness, 
                                                leaf_sectors = leaf_sectors)
        dus_s.append(s)
        dus_t.append(t)
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(interleaves, dus_s, 'b', linestyle = linestyle)
    ax.plot(interleaves, dus_t, 'r', linestyle = linestyle)
    ax.legend(['source leaf', 'target leaf'], loc = 'best')
    ax.set_ylabel('Number of DUs intercepted', fontsize = 14)
    ax.set_ylabel('Distance between leaves (cm)', fontsize = 14)
    if return_ax == True:
        return ax
        
def plot_layer_thickness_x_interleaf_two_metamers():
    interleaf = numpy.arange(1, 31, 1)
    layer_thickness = numpy.arange(0.01, 0.31, 0.01)
    df_source = pandas.DataFrame(index = interleaf, columns = layer_thickness)
    df_target = pandas.DataFrame(index = interleaf, columns = layer_thickness)
    for int_lf in df_source.index:
        for dh in df_source.columns:
            print int_lf, dh
            (df_source.loc[int_lf, dh],
             df_target.loc[int_lf, dh]) = test_transport_rain_two_metamers(interleaf=float(int_lf), 
                                                                             layer_thickness=dh)
            
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    X, Y = numpy.meshgrid(df_source.columns, df_source.index)
    ax.plot_surface(X,Y, df_source.values, 
                    rstride=1, cstride=1, color='b', alpha=0.5)
    ax.set_xlabel('Layer thickness (m)')
    ax.set_ylabel('inter leaf (cm)')
    ax.set_zlabel('normalised number of deposits on source leaf')
    # ax.set_zlim([0,1])
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    X, Y = numpy.meshgrid(df_target.index, df_target.columns)
    ax.plot_surface(X,Y, df_target.values, 
                    rstride=1, cstride=1, color='b', alpha=0.5)
    ax.set_xlabel('Layer thickness (m)')
    ax.set_ylabel('Inter leaf (cm)')
    ax.set_zlabel('normalised number of deposits on target leaf')
    # ax.set_zlim([0,1])
    
def plot_nb_sectors_x_interleaf_two_metamers():
    interleaf = numpy.arange(1, 31, 1)
    nb_sectors = numpy.arange(1, 11, 1)
    df_source = pandas.DataFrame(index = interleaf, columns = nb_sectors)
    df_target = pandas.DataFrame(index = interleaf, columns = nb_sectors)
    for int_lf in df_source.index:
        for nb_sect in df_source.columns:
            print int_lf, nb_sect
            (df_source.loc[int_lf, nb_sect],
             df_target.loc[int_lf, nb_sect]) = test_transport_rain_two_metamers(interleaf=float(int_lf), 
                                                                                leaf_sectors=nb_sect)
            
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    X, Y = numpy.meshgrid(df_source.columns, df_source.index)
    ax.plot_surface(X,Y, df_source.values/max(df_source.max()), 
                    rstride=1, cstride=1, color='b', alpha=0.5)
    ax.set_xlabel('Number of sectors by leaf')
    ax.set_ylabel('inter leaf (cm)')
    ax.set_zlabel('normalised number of deposits on source leaf')
    ax.set_zlim([0,1])
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    X, Y = numpy.meshgrid(df_target.columns, df_target.index)
    ax.plot_surface(X,Y, df_target.values/max(df_target.max()), 
                    rstride=1, cstride=1, color='b', alpha=0.5)
    ax.set_xlabel('Number of sectors by leaf')
    ax.set_ylabel('inter leaf (cm)')
    ax.set_zlabel('normalised number of deposits on target leaf')
    ax.set_zlim([0,1])
    
def plot_sectors_two_metamers(by_sector = False):
    if by_sector == False:
        df = pandas.DataFrame(columns = ['source', 'target'])
        for nb_sect in range(1,11):
            print nb_sect
            df.loc[nb_sect, :] = test_transport_rain_two_metamers(interleaf=5.,
                                                                  layer_thickness=0.01, 
                                                                  leaf_sectors=nb_sect,
                                                                  by_sector = False)

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
            (nb_dus_on_source, 
             nb_dus_on_target) = test_transport_rain_two_metamers(interleaf=5.,
                                                                  layer_thickness=0.01,
                                                                  leaf_sectors=nb_sect,
                                                                  by_sector=True)
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

# Examples in stand ###########################################################
from alinea.echap.architectural_reconstructions import soisson_reconstruction

def test_transport_wind_canopy(nplants=30, density=250.,
                               inter_row=0.15, nsect=3, 
                               age_canopy=1400., position_source = 3./5):
    # Generate canopy                                   
    adel = soisson_reconstruction(nplants=nplants, density=density,
                                  inter_row=inter_row, nsect=nsect)
    g = adel.setup_canopy(age_canopy)
    
    # Get source leaf
    

