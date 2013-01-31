""" Test the protocol of communication between adel and alep """

# Imports #########################################################################

import random

from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe

from alinea.alep import cycle2
from alinea.alep.cycle2 import septoria
from alinea.alep.cycle2 import SeptoriaDU
from alinea.alep.cycle2 import powdery_mildew

# from alinea.alep.Septo3DDispersion import SplashDispersal

from alinea.alep.protocol import *

# Construction ####################################################################

def adelR(nplants,dd):
    devT = devCsv('../../adel/test/data/axeTCa0N.csv','../../adel/test/data/dimTCa0N.csv','../../adel/test/data/phenTCa0N.csv','../../adel/test/data/earTCa0N.csv','../../adel/test/data/ssi2sen.csv')
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants)
    cantable = RunAdel(dd,pars)
    return pars,cantable

def leaves_db():
    import cPickle as Pickle
    fn = r'../../adel/adel/data/leaves_simple.db'
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves = fitting.fit_leaves(leaves, 9)
    return leaves

def adel_mtg():
    """ create a very simple adel mtg """
    d = {'plant':[1,1],'axe':[0,1],'numphy':[1,1], 
         'Laz': [0,90], 'Ll' :[3,3], 'Lv' :[3,3] , 'Lsen':[0,0], 'L_shape':[3,3], 'Lw_shape':[.3,.3], 'Linc':[0,0],
         'Einc':[0,45],'El':[1,1],'Ev':[1,1],'Esen':[0,0],'Ed': [0.1,0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=leaves_db())
    g=mtg_interpreter(g)
    return g
    
def adel_mtg2():
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    g=mtg_factory(d,adel_metamer,leaf_db=leaves_db(),stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
    g=mtg_interpreter(g)
    return g
   
def update_climate(g, label = 'LeafElement'):
    """ simulate an environmental program """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.wetness = True
        n.temp = 18.
        n.age = 1.
        n.healthy_surface = 10.
        n.rain_intensity = 0.
        n.relative_humidity = 85.

def update_climate_with_rain(g, label = 'LeafElement'):
    """ simulate an environmental program with rain """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.wetness = True
        n.temp = 18.
        n.age = 1.
        n.healthy_surface = 10.
        n.rain_intensity = 5.
        n.relative_humidity = 85.

# Representation ##################################################################

def plot_lesions(g):
    """ plot the plant with infected elements in red """
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)  

def plot_DU(g):
    """ plot the plant with elements carrying dispersal units in red """
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'dispersal_units' in n.properties():
            n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)

def plot_lesions_after_DU(g):
    """ plot the plant with :
        - elements carrying dispersal units in yellow
        - infected elements in red
        
    """
    green = (0,180,0)
    yellow = (247, 220, 17)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            n.color = red
        else :
            if 'dispersal_units' in n.properties():
                n.color = yellow
            else : 
                n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)

def plot_lesions_in_state(g, state):
    """ plot the plant with elements carrying lesions in given state in red """
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            for les in n.lesions:
                if les.status == state:
                    n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)

# Tests ########################################################################### 

def test_adel_mtg():
    """ Check the proper functioning of 'adel_mtg'.
    
    """
    g = adel_mtg()
    scene = plot3d(g)
    Viewer.display(scene)
    
def test_adel_mtg2():
    """ Check the proper functioning of 'adel_mtg2'.
    
    """
    g = adel_mtg2()
    scene = plot3d(g)
    Viewer.display(scene)

def test_initiate():
    """ Check if 'initiate' from 'protocol.py' deposits dispersal units on the MTG.
    
    """
    g = adel_mtg2()
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    plot_DU(g)
    return g
    
def test_infect():
    """ Check if 'infect' from 'protocol.py' leads to infection by dispersal units on the MTG.

    """
    g = adel_mtg2()
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    dt = 1
    nb_steps = 50
    for i in range(nb_steps):
        update_climate(g)
        infect(g, dt)
    
    plot_lesions_after_DU(g)
    return g
    
def test_update():
    """ Check if 'update' from 'protocol.py' provokes the growth of a lesion instantiated on the MTG.

    """
    g = adel_mtg2()
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(1000)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    dt = 1
    nb_steps = 1000
    for i in range(nb_steps):
        print('time step %d' % i)
        
        if i%100 == 0:
            update_climate_with_rain(g)
        else:
            update_climate(g)
            
        #grow(g)
        infect(g, dt)        
        update(g,dt)
        # lesions = g.property('lesions')
        # for vid,l in lesions.iteritems():
            # for lesion in l:
                # print('statut = %d' % lesion.status)
                # print('nb rings = %d' % len(lesion.rings))
                # print([ring.age_dday for ring in lesion.rings])
        
    plot_lesions_after_DU(g)
    return g
   

class StubDispersal(object):

    def __init__(self):
        pass
        
        
    def disperse(self,scene, DU):
        print DU
        import random
        vids = [geom.id for geom in scene]
        n = len(vids)
        deposits = {}
    
        for vid,dlist in DU.iteritems():
            for d in dlist:
                idx = random.randint(0,n-1)
                v = vids[idx]
                if v not in deposits:
                    deposits[v] = []
                deposits[v].append(d)
        
        print deposits
        return(deposits)
            
            
def test_disperse():
    """ Check if 'disperse' from 'protocol.py' disperse new dispersal units on the MTG.

    """
    g = adel_mtg2()
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(10)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    dispersor = StubDispersal()
    
    dt = 1
    nb_steps = 1000
    for i in range(nb_steps):
        print('time step %d' % i)
        
        if i%100 == 0:
            update_climate_with_rain(g)
            rain = 10
        else:
            update_climate(g)
            rain=0
            
            
        # grow(g)
        infect(g, dt)
        update(g,dt)
        if rain != 0:
            scene = plot3d(g)
            disperse(g, scene, dispersor, "septoria")
        
    plot_lesions_after_DU(g)
    return g
    