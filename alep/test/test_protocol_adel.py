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
from alinea.alep.cycle2 import proba

# from alinea.alep.Septo3DDispersion import SplashDispersal

from alinea.alep.protocol import *

# Plant ###########################################################################

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

# Climate #########################################################################

# def create_climate_list():
    # wetness = [True]

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

def rain_interception(g, rain_intensity = None, label = 'LeafElement'):
    """ simulate an environmental program with rain """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.rain_intensity = rain_intensity

class RandomInoculation:
    """ Example of class created to allocate inoculum on MTG randomly.
    
    """
    
    def disperse(self, g, inoculum, label='LeafElement'):
        """ Select randomly elements of the MTG and allocate them a random part of the inoculum.

        :Parameters:
          - `g` : MTG representing the canopy (and the soil).
          - `inoculum` : source of dispersal units to disperse in the scene.
        """        
        import random
        vids = [n for n in g if g.label(n).startswith(label)]
        n = len(vids)
        for vid in vids:
            g.node(vid).dispersal_units = []
            
        for i in inoculum:
            idx = random.randint(0,n-1)
            v = vids[idx]
            # Deposit a DU from inoculum on node v of the MTG  
            i.deposited()
            g.node(v).dispersal_units.append(i)
            
class StubDispersal(object):

    def __init__(self):
        pass

    def disperse(self, scene, DU):
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
         
        return deposits
        
class StubWashing(object):
    def __init__(self):
        pass
        
    def wash(self, DU, leaf):
        rain_intensity = leaf.rain_intensity
        # rain_duration = leaf.rain_duration
        healthy_surface = leaf.healthy_surface
        rain_duration = 1 # temp
        
        if rain_duration == 0.:
            washing_rate = 0.
        else:
            washing_rate = max(0,min(1, rain_intensity / (healthy_surface + rain_intensity)))
        
        if proba(washing_rate):
            DU.inactive()            
            
# Display #########################################################################

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
    """ plot the plant with elements carrying dispersal units in yellow """
    green = (0,180,0)
    yellow = (247, 220, 17)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'dispersal_units' in n.properties():
            n.color = yellow
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
 
class DisplayLesions(object):
    """ Print the ID of Leaf Elements where new lesions appear. """
    
    def __init__(self):
        self.old_lesions = []
        
    def print_new_lesions(self, g):
        lesions = g.property('lesions')
        
        for vid, l in lesions.iteritems():
            for lesion in l:
                if not lesion in self.old_lesions:
                    self.old_lesions.append(lesion)
                    print('New Lesions on : ' + g.label(vid) + ' %d' % vid)
        
    def print_all_lesions(self, g):
        from pprint import pprint
        lesions = g.property('lesions')
        ldict = {}
        for vid, l in lesions.iteritems():
            for lesion in l:
                if vid not in ldict:
                    ldict[vid] = 0
                ldict[vid] += 1
        print('\n' + 'Number of lesions by leaf element : ' + '\n')             
        pprint(ldict)
        # print('You can find lesions on LeafElements : ' + str(llist).strip('[]'))
    
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
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    dt = 1
    nb_steps = 1000
    for i in range(nb_steps):
        print('time step %d' % i)
        
        update_climate(g)
        if i%100 == 0:
            rain_interception(g, rain_intensity = 5.)          
            
        #grow(g)
        infect(g, dt)        
        update(g,dt)
        # lesions = g.property('lesions')
        # for vid,l in lesions.iteritems():
            # for lesion in l:
                # print('statut = %d' % lesion.status)
                # print('nb rings = %d' % len(lesion.rings))
                # print([ring.age_dday for ring in lesion.rings])
    
    displayer = DisplayLesions()
    displayer.print_all_lesions(g)
    
    plot_lesions_after_DU(g)
    return g

def test_disperse():
    """ Check if 'disperse' from 'protocol.py' disperse new dispersal units on the MTG.

    """
    g = adel_mtg2()
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(2)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    displayer = DisplayLesions()
    
    dt = 1
    nb_steps = 1000
    for i in range(nb_steps):
        print('time step %d' % i)
        
        update_climate(g)
        if i%100 == 0:
            rain_intensity = 5.
            rain_interception(g, rain_intensity)
        else:
            rain_intensity = 0.
            
        # grow(g)
        infect(g, dt)
        update(g,dt)
        if rain_intensity != 0:
            scene = plot3d(g)
            dispersor = StubDispersal()
            disperse(g, scene, dispersor, "Septoria")

        displayer.print_new_lesions(g)
    
    print('-----------------------------------')
    displayer.print_all_lesions(g)
    plot_lesions(g)
    
    return g

def test_washing():
    """ Check if 'washing' from 'protocol.py' washes dispersal units out of the MTG.

    """
    g = adel_mtg2()
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    washor = StubWashing()
    
    dt = 1
    nb_steps = 100
    for i in range(nb_steps):
    
        update_climate(g)
        if i%5 == 0:
            rain_intensity = 5.
            rain_interception(g, rain_intensity)
        else:
            rain_intensity = 0.
    
        if rain_intensity != 0:
            wash(g, washor)
        
            dispersal_units = g.property('dispersal_units')
            nbDU = 0.
            for vid, dlist in dispersal_units.iteritems():
                for d in dlist:
                    if d.active:
                        nbDU += 1
            print('\n')
            print('   _ _ _')
            print('  (pluie)')
            print(' (_ _ _ _)')
            print('   |  |')
            print('    |  | ')
            print('\n')

            print('Sur le MTG il y a %d DU actives en tout' % nbDU)
        # plot_DU(g)
    
    return g