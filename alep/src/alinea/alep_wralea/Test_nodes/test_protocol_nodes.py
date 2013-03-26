""" Test the protocol of communication between adel and alep """

# Imports #########################################################################

import random
import numpy
import pandas
from pylab import *

from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe

from alinea.alep import cycle2
from alinea.alep import septoria
from alinea.alep.septoria import *
from alinea.alep import powdery_mildew

from alinea.alep.dispersal import RandomDispersal
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.inoculation import RandomInoculation

from alinea.alep.protocol import *

from datetime import datetime
import time
from math import ceil, sqrt

# Plant ###########################################################################

def adelR(nplants,dd):
    devT = devCsv('../../adel/example/data/axeTCa0N.csv','../../adel/example/data/dimTCa0N.csv','../../adel/example/data/phenTCa0N.csv','../../adel/example/data/earTCa0N.csv','../../adel/example/data/ssi2sen.csv')
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
    leaves,discard = fitting.fit_leaves(leaves, 9)
    return leaves
   
def leaves_db_flow(fn):
    import cPickle as Pickle
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves,discard = fitting.fit_leaves(leaves, 9)
    return leaves

def adel_mtg():
    """ create a very simple adel mtg """
    d = {'plant':[1,1],'axe_id':['MS','T1'],'ms_insertion':[0,1],'numphy':[1,1], 
         'Laz': [0,90], 'Ll' :[3,3], 'Lv' :[3,3] , 'Lsen':[0,0], 'L_shape':[3,3], 'Lw_shape':[.3,.3], 'Linc':[0,0],
         'Einc':[0,45],'El':[1,1],'Ev':[1,1],'Esen':[0,0],'Ed': [0.1,0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=leaves_db(), leaf_sectors=1)
    g=mtg_interpreter(g)
    return g
    
def adel_mtg2(nb_sect=1):
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaves_db(),stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
    g=mtg_interpreter(g)
    return g
    
def adel_mtg3(nb_sect=1, leaf_db=None, d=None, p=None):
    """ create a less simple adel mtg """
    if p: # nb_plants
        size = int(ceil(sqrt(p)))
        stand = numpy.array([(i, j) for i in range(size) for j in range(size)])
        numpy.random.shuffle(stand)
        stand = [((int(i)-10*size/2., int(j)-10*size/2., 0),random.randint(0,90)) for i, j in 10*stand[:p]]
    else:
        stand = [((0,0,0),0),((10,0,0),90), ((0,10,0), 0)]
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaf_db,stand=stand)
    g=mtg_interpreter(g)
    return g

# Climate #########################################################################

def read_weather(filename):
    """ Use pandas to save data from a .txt file into a dataframe.
    
    """
    raw_data = pandas.read_table(filename)
    
    # Addition of wetness variable :
    wet = dict(wetness=[])
    for i_line in range(len(raw_data.An)):
        if raw_data.Pluie[i_line] > 0. or (raw_data.PAR[i_line] < 644. and raw_data.HR[i_line] > 85.):
            wet['wetness'].append(True)
        else:
            wet['wetness'].append(True)
    wetness = pandas.DataFrame(wet)
    
    weather_data = raw_data.join(wetness)
    
    # Creation of date time indexes :
    first_date = datetime.fromordinal(datetime(weather_data.An[0], 1, 1).toordinal() + weather_data.Jour[0] - 1)
    rng = pandas.date_range(first_date, periods = len(weather_data.An), freq ='H')
    
    weather_data.index = rng
    
    return weather_data

def update_on_leaves(data, date, g, label = 'LeafElement'):
    """ Read weather data for a step of simulation and apply it to each leaf.
    
    """        
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        
        n.temp = data.Tair[date]
        n.rain_intensity = data.Pluie[date]
        n.relative_humidity = data.HR[date]
        n.wetness = data.wetness[date]

    return g
    
def initiate_g(g, label = 'LeafElement'):
    """ Give initial values for plant properties of each LeafElement. 
    
    """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.age = 1.
        n.surface = 10.
        n.healthy_surface = n.surface # TODO : Manage properly
        
    return g

def update_climate(g, label = 'LeafElement'):
    """ simulate an environmental program """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.wetness = True
        n.temp = 18.
        n.rain_intensity = 0.
        n.relative_humidity = 85.
    
    return g

def rain_interception(g, rain_intensity = None, label = 'LeafElement'):
    """ simulate an environmental program with rain """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.rain_intensity = rain_intensity
        
    return g

# Stub models for disease #########################################################

def generate_stock_DU(fungus = septoria(), nb_Spores=random.randint(1,100), nb_DU = 100):
    """ Generate a stock of DU as a list of DU 
    
    """
    SeptoriaDU.fungus = fungus
    stock_DU = [SeptoriaDU(nb_spores=nb_Spores, status='emitted') for i in range(nb_DU)]
    return stock_DU

def inoculator():            
    """ Instantiate the class RandomInoculation().
    
    """
    inoculator = RandomInoculation()
    return inoculator

def dispersor():            
    """ Instantiate the class RandomInoculation().
    
    """
    dispersor = RandomDispersal()
    return dispersor
    
def washor():            
    """ Instantiate the class RandomInoculation().
    
    """
    washor = RapillyWashing()
    return washor

# Display #########################################################################

def temp_plot3D(g):
    """ plot g """
    scene = plot3d(g)
    Viewer.display(scene)


def plot_DU(g):
    """ plot the plant with elements carrying dispersal units in yellow """
    green = (0,180,0)
    yellow = (247, 220, 17)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'dispersal_units' in n.properties() and n.dispersal_units:
            n.color = yellow
        else : 
            n.color = green
    
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

def plot_DU_deposited(g):
    """ plot the plant with elements carrying deposited dispersal units in yellow """
    green = (0,180,0)
    yellow = (247, 220, 17)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        ind_deposited = 0.
        if 'dispersal_units' in n.properties():
            for du in n.dispersal_units:
                if du.status == 'deposited':
                    ind_deposited += 1
            
        if ind_deposited > 0 :
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

def count_lesions_in_state(g, state):
    """ count lesions of the mtg in give state.
    
    """
    count = 0.
    lesions = g.property('lesions')
    return sum(1 for l in lesions.itervalues() for lesion in l if lesion.status == state)
    
def count_DU(g):
    """ count DU of the mtg.
    
    """
    count = 0.
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'dispersal_units' in n.properties():
            for du in n.dispersal_units:
                if du.active:
                    count +=1
    return count

def count_lesions(g):
    """ count DU of the mtg.
    
    """
    count = 0.
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            for l in n.lesions:
                count +=1
    return count
    
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
        
    def print_lesion_surfaces(self, g):
        from pprint import pprint
        lesions = g.property('lesions')
        ldict = {}
        for vid, l in lesions.iteritems():
            for lesion in l:
                if vid not in ldict:
                    ldict[vid] = 0
                ldict[vid] += lesion.surface
        print('\n' + 'Sum of lesion surfaces by leaf element : ' + '\n')             
        pprint(ldict)
    
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

def test_adel_mtg3():
    """ Check the proper functioning of 'adel_mtg2'.
    
    """
    p, d = adelR(3,1000)
    g = adel_mtg3(leaf_db=leaves_db(fn = r'../../adel/adel/data/leaves_simple.db'), d=d, p=p)
    scene = plot3d(g)
    Viewer.display(scene)
    
def test_initiate():
    """ Check if 'initiate' from 'protocol.py' deposits dispersal units on the MTG.
    
    """
    g = adel_mtg2()
    initiate_g(g)
    fungus = septoria()
    SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    plot_DU(g)
    return g
    
def test_infect():
    """ Check if 'infect' from 'protocol.py' leads to infection by dispersal units on the MTG.

    """
    g = adel_mtg2()
    initiate_g(g)
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    dt = 1
    nb_steps = 100
    plot_DU(g)
    for i in range(nb_steps):
        update_climate(g)
        infect(g, dt)
            
    plot_lesions(g)
    return g
       
def test_update():
    """ Check if 'update' from 'protocol.py' provokes the growth of a lesion instantiated on the MTG.

    """
    g = adel_mtg2()
    g = initiate_g(g)
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    dt = 1
    nb_steps = 750
    nb_les_inc = numpy.zeros(nb_steps)
    # nb_les_chlo = numpy.array([0. for i in range(nb_steps)])
    # nb_les_nec = numpy.array([0. for i in range(nb_steps)])
    nb_les_spo = numpy.zeros(nb_steps)
    # nb_les_empty = numpy.array([0. for i in range(nb_steps)])
    for i in range(nb_steps):

        ts = time.clock()
        #print('time step %d' % i)
        
        update_climate(g)
        # if i%100 == 0:
            # global_rain_intensity = 4.
            # rain_interception(g, rain_intensity = global_rain_intensity*0.75)  
        # else:
            # global_rain_intensity = 0.
            
        #grow(g)
        infect(g, dt)        
        update(g,dt)
        
        # Count of lesions :
        nb_les_inc[i] = count_lesions_in_state(g, state = 0)
        # nb_les_chlo[i] = count_lesions_in_state(g, state = 1)
        # nb_les_nec[i] = count_lesions_in_state(g, state = 2)
        nb_les_spo[i] = count_lesions_in_state(g, state = 3)
        # nb_les_empty[i] = count_lesions_in_state(g, state = 4)

        te = time.clock()
        #print "time ", i, " : ", te-ts
    
    # Display results
    plot(nb_les_inc)
    plot(nb_les_spo)
    ylabel('Nombre de lesions dans cet etat sur le MTG')
    xlabel('Pas de temps de simulation')
    ylim([0, 120])
       
    # displayer = DisplayLesions()
    # displayer.print_all_lesions(g)

    # plot_lesions(g)
    return g

def test_disperse():
    """ Check if 'disperse' from 'protocol.py' disperse new dispersal units on the MTG.

    """
    g = adel_mtg2()
    initiate_g(g)
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=random.randint(1,100), status='emitted') for i in range(10)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    dt = 1
    nb_steps = 750
    nb_les = numpy.array([0 for i in range(nb_steps)])
    for i in range(nb_steps):
        print('time step %d' % i)
        
        update_climate(g)
        if i%100 == 0:
            global_rain_intensity = 4.
            rain_interception(g, rain_intensity = global_rain_intensity*0.75)  
        else:
            global_rain_intensity = 0.
        
        # grow(g)
        infect(g, dt)
        update(g,dt)
        
        if global_rain_intensity != 0.:
            scene = plot3d(g)
            disperse(g, scene, dispersor(), "Septoria")
               
        # Count of lesions :
        nb_les[i] = count_lesions(g)
        
        # Display results
        plot_lesions(g)
        
        # displayer = DisplayLesions()
        # displayer.print_new_lesions(g)
    
    plot(nb_les)
    ylabel('Nombre de lesions sur le MTG')
    xlabel('Pas de temps de simulation')
    show()
    
    return g

def test_washing():
    """ Check if 'washing' from 'protocol.py' washes dispersal units out of the MTG.

    """
    g = adel_mtg()
    initiate_g(g)
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    washor = RapillyWashing()
    
    dt = 1
    nb_steps = 100
    nb_DU = numpy.array([0. for i in range(nb_steps)])
    for i in range(nb_steps):
        g = update_climate(g)
        
        # Add artificial rain:
        if i>2 and i%5 == 0 or (i-1)%5 == 0:
                global_rain_intensity = 4.
                rain_interception(g, rain_intensity = global_rain_intensity*0.75)      
        else:
            global_rain_intensity = 0.
        
        # Compute washing 
        # (needs to be done even if no rain to update variables in the washing model)
        wash(g, washor, global_rain_intensity, DU_status='deposited')
           
        # Count of DU :
        nb_DU[i] = count_DU(g)

        # Display results
        if global_rain_intensity != 0. :
            print('\n')
            print('   _ _ _')
            print('  (pluie)')
            print(' (_ _ _ _)')
            print('   |  |')
            print('    |  | ')
            print('\n')
            print('Sur le MTG il y a %d DU actives en tout' % nb_DU[i])

    # Display results
    plot(nb_DU)
    ylim([0, 120])
    ylabel('Nombre de DU sur le MTG')
    xlabel('Pas de temps de simulation')
    show()
    
    return g

def test_growth_control():
    g = adel_mtg()
    initiate_g(g)
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=random.randint(1,100), status='emitted') for i in range(1000)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    controler = NoPriorityGrowthControl()

    dt = 1
    nb_steps = 500
    sum_surface = numpy.array([0. for i in range(nb_steps)])
    for i in range(nb_steps):
        print('time step %d' % i)
        
        g = update_climate(g)
        if i%100 == 0:
            global_rain_intensity = 4.
            rain_interception(g, rain_intensity = global_rain_intensity*0.75)  
        else:
            global_rain_intensity = 0.
            
        #grow(g)
        infect(g, dt)        
        update(g,dt)
        growth_control(g, controler)
        
        vids = [v for v in g if g.label(v).startswith("LeafElement")]
        count_surf = 0.
        for v in vids:
            leaf = g.node(v)
            count_surf += leaf.healthy_surface
            # if 'lesions' in leaf.properties():
                # count_surf += sum([l.surface for l in leaf.lesions])
                # print('leaf element %d - ' % v + 'Healthy surface : %f'  % leaf.healthy_surface) 
        
        sum_surface[i] = count_surf
        
    # Display results:
    plot(sum_surface)
    ylabel('Surface saine totale sur le MTG')
    xlabel('Pas de temps de simulation')
    show()
        
    return g
        
def test_simul_with_weather():
    """ Check if 'update' from 'protocol.py' provokes the growth of a lesion instantiated on the MTG.

    """
    # Read weather data : 
    filename = 'meteo01.txt'
    weather_data = read_weather(filename)
    
    # Generate a MTG with required properties :
    g = adel_mtg2()
    initiate_g(g)
    
    # Deposit first dispersal units on the MTG :
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

  
    # Prepare the simulation loop
    dt = 1
    start_date = datetime(2000, 10, 1)
    end_date = datetime(2000, 12, 31)
    # end_date = datetime(2001, 7, 31)
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        
        # Update climate on leaf elements : 
        update_on_leaves(weather_data, date, g)     
        
        # Run a disease cycle
        #grow(g)
        infect(g, dt)        
        update(g,dt)
        growth_control(g, controler)
               
        global_rain_intensity = weather_data.Pluie[date]
        if global_rain_intensity != 0. :
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria") 
        wash(g, washor, global_rain_intensity)
    
        # displayer = DisplayLesions()
        # displayer.print_all_lesions(g)
    
    # Display results:
    plot_lesions(g)
    return g
    
def test_all():
    
    # Generate a MTG with required properties :
    g = adel_mtg2()
    initiate_g(g)
    
    # Deposit first dispersal units on the MTG :
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=random.randint(1,100), status='emitted') for i in range(10)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()
   
    # Prepare the simulation loop
    dt = 10
    nb_steps = 750

    for i in range(0,nb_steps,dt):
        update_climate(g)
        if i%100 == 0:
            global_rain_intensity = 4.
            rain_interception(g, rain_intensity = global_rain_intensity*0.75)  
        else:
            global_rain_intensity = 0.
        
        # grow(g)
        infect(g, dt)
        update(g,dt)
        growth_control(g, controler)
        
        if global_rain_intensity != 0.:
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria")
        
        wash(g, washor, global_rain_intensity, DU_status='deposited')
        
    return g
