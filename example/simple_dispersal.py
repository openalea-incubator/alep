from alinea.echap.architectural_reconstructions import EchapReconstructions
from alinea.alep.protocol import *
from alinea.alep.architecture import get_leaves
from alinea.alep.disease_outputs import plot_severity_by_leaf
from alinea.alep.septo3d_v2 import plugin_septoria
from alinea.alep.inoculation import InoculationLowerLeaves
from alinea.alep.fungal_objects import Lesion, DispersalUnit
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport
from septo_decomposed import get_weather
from alinea.alep.alep_weather import add_rain_dispersal_events
from datetime import timedelta
import numpy
import matplotlib.pyplot as plt

class DummyLesion(Lesion):
    """ Undetermined lesion with only a method emission."""
    def __init__(self):
        super(DummyLesion, self).__init__()
        class params():
            def __init__(self, name="dummy"):
                self.name = name
        self.fungus = params()
        
    def is_sporulating(self):
        return True
        
class DummyEmission():
    """ Emission model with forced emission of dispersal units. """
    def __init__(self, number = 1e3):
        self.nb = number
        
    def get_dispersal_units(self, g, fungus_name="dummy", label='LeafElement'):
        DU={}
        lesions = {k:[l for l in les if l.fungus.name is fungus_name and l.is_sporulating()] 
                    for k, les in g.property('lesions').iteritems()} 
        for vid, l in lesions.iteritems():
            for lesion in l:
                emissions = [DispersalUnit() for i in range(int(self.nb))]
                try:
                    DU[vid] += emissions
                except:
                    DU[vid] = emissions
        return DU

def tremie12_main_stem(nplants = 1, nsect = 1, age=1500, only_MS = True):
    reconst = EchapReconstructions()
    adel = reconst.get_reconstruction(name='Tremie12', nplants = nplants, nsect = nsect, disc_level = 5, aspect = 'line')
    g = adel.setup_canopy(age=age)
    labels = g.property('label')
    vids = labels.keys()
    vids.reverse()
    if only_MS == True:
        for vid in vids:
            if g.node(vid).complex_at_scale(2).label.startswith('T'):
                g.remove_tree(vid)
    return g, adel

def get_source_leaf(g):
    leaves = get_leaves(g)
    geometries = g.property('geometry')
    leaves = [lf for lf in leaves if lf in geometries]
    return g.node(leaves[2])
    
def dummy_dispersal(nplants=1, nsect = 1, age = 1500., only_MS = True):
    # Create wheat
    g, adel = tremie12_main_stem(nplants = nplants, nsect = nsect, age = age, only_MS = only_MS)    
        
    # Initiate models for emission and dispersal
    emitter = DummyEmission(number = 1e3)
    transporter = PopDropsTransport(domain = adel.domain, domain_area = adel.domain_area)

    # Deposit emitting lesion on source leaf
    source_leaf = get_source_leaf(g)
    source_leaf.lesions = [DummyLesion()]

    DU = emitter.get_dispersal_units(g, "dummy", 'LeafElement')
    deposits = transporter.disperse(g, DU, weather_data)
    
    if only_MS == False:
        MS_ids = [vid for vid in g if vid>0 and not g.node(vid).complex_at_scale(2).label.startswith('T')]
        deposits = {k:v for k,v in deposits.iteritems() if k in MS_ids}
    
    return len(sum(deposits.values(), []))

# Find a rain sequence in weather
weather = get_weather()
seq = weather.data[weather.data.rain>0]
rain_seqs = numpy.where(seq.index[:-1]+timedelta(hours=1)!=seq.index[1:])[0]
start_date = rain_seqs[numpy.diff(rain_seqs).argmax()]
end_date = rain_seqs[numpy.diff(rain_seqs).argmax()+1]
weather_data = seq[start_date:end_date]

def test_transport():
    deposits = []
    # nb_plants = [1, 5, 10, 15, 20, 25, 30]
    nb_plants = [1, 10, 20, 30]
    # nb_sects = [1, 3, 5, 7, 10]
    # ages = [500, 1000, 1500, 2000]
    for npl in nb_plants:
        print '--------------------------'
        print 'number of plant: %d' % npl
        # print 'age: %d' % age
        deposits.append(dummy_dispersal(nplants=npl, nsect=5, age=1500, only_MS = True))
    return deposits
    
def test_emission(nlesions = 100, nplants=1, nsect = 1, age = 1500., only_MS = True):
    nplants = 1
    nsect = 1
    age = 1500.
    only_MS = True
    
    # Create wheat
    g, adel = tremie12_main_stem(nplants = nplants, nsect = nsect, age = age, only_MS = only_MS)    

    # Deposit emitting lesion on source leaf 
    septoria = plugin_septoria()
    source_leaf = get_source_leaf(g)
    lesion = septoria.lesion()
    lesion.surfaces_spo = numpy.array([0.05, 0., 0.])
    lesion.status = 3
    lesions = [lesion for i in range(nlesions)]
    initiation_model = InoculationLowerLeaves()
    initiation_model.allocate(g, lesions, 'LeafElement')

    emitter = PopDropsEmission(domain=adel.domain)
    DU = emitter.get_dispersal_units(g, "septoria", 'LeafElement', weather_data)
    return DU