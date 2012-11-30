""" Test the protocol of communication between adel and alep """

from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting

from alinea.alep import cycle
from alinea.alep.cycle import septoria
from alinea.alep.cycle import powdery_mildew

from alinea.alep.protocol import *


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
    
    
def test_adel():
    g = adel_mtg()
    scene = plot3d(g)
    Viewer.display(scene)    

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

    
    
def test_initiate():
    g = adel_mtg()
    les = cycle.LesionFactory(fungus = septoria())
    initiate(g,les)
    plot_lesions(g)
    return g
    
def test_update():
    g = adel_mtg()
    update_climate(g)
    les = cycle.LesionFactory(fungus = septoria())
    initiate(g,les)
    dt = 1
    nb_steps = 1000
    for i in range(nb_steps):
        #grow(g)
        update_climate(g)   
        update(g,dt)
    return g
    