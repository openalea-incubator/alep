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
from alinea.alep.cycle2 import septoria
from alinea.alep.cycle2 import SeptoriaDU
from alinea.alep.cycle2 import powdery_mildew
from alinea.alep.cycle2 import proba

from alinea.alep.protocol import *

from datetime import datetime

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
    leaves = fitting.fit_leaves(leaves, 9)
    return leaves
   
def leaves_db_flow(fn):
    import cPickle as Pickle
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
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaf_db,stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
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
    stock_DU = [SeptoriaDU(fungus = fungus, nbSpores=nb_Spores, status='emitted') for i in range(nb_DU)]
    return stock_DU

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
            leaf = g.node(v)
            if not 'dispersal_units' in leaf.properties():
                leaf.dispersal_units = []
            # Deposit a DU from inoculum on node v of the MTG
            i.deposited()
            leaf.dispersal_units.append(i)

def inoculator():            
    """ Instantiate the class RandomInoculation().
    
    """
    # Temp : Do better. How to instantiate a class in Dataflow ?
    inoculator = RandomInoculation()
    return inoculator

def dispersor():            
    """ Instantiate the class RandomInoculation().
    
    """
    # Temp : Do better. How to instantiate a class in Dataflow ?
    dispersor = StubDispersal()
    return dispersor
    
def washor():            
    """ Instantiate the class RandomInoculation().
    
    """
    # Temp : Do better. How to instantiate a class in Dataflow ?
    washor = StubWashing()
    return washor
    
class StubDispersal(object):

    def __init__(self):
        pass

    def disperse(self, scene, DU):
        vids = [geom.id for geom in scene]
        n = len(vids)
        deposits = {}
    
        for vid,dlist in DU.iteritems():
            for d in dlist:
                if proba(0.5):
                    idx = random.randint(0,n-1)
                    v = vids[idx]
                    if v not in deposits:
                        deposits[v] = []
                    deposits[v].append(d)
         
        return deposits
        
class StubWashing(object):

    def __init__(self):
        self.rain_duration = 0.
        self.rain_intensity = {}
        
    def compute_washing_rate(self, g, global_rain_intensity, label='LeafElement'):
        """ For each LeafElement of the MTG, compute the washing rate.
        
        """
        if global_rain_intensity > 0.:
            self.rain_duration += 1
            
            vids = [n for n in g if g.label(n).startswith(label)]
            for v in vids : 
                leaf = g.node(v)
                if not v in self.rain_intensity:
                    self.rain_intensity[v] = []
                self.rain_intensity[v].append(leaf.rain_intensity)
                leaf.washing_rate = 0.
        else:
            vids = [n for n in g if g.label(n).startswith(label)]
            for v in vids : 
                leaf = g.node(v)
                if v in self.rain_intensity:
                    mean_rain_intensity = numpy.mean(self.rain_intensity[v])
                    leaf.washing_rate = max(0,min(1, mean_rain_intensity / (leaf.healthy_surface + mean_rain_intensity)*self.rain_duration))
                    self.rain_intensity.pop(v)
                else:
                    leaf.washing_rate = 0.
            self.rain_duration = 0.
        
    
    # def wash(self, dispersal_unit, washing_rate):
        # """ On the given LeafElement, disable the DU as a function of the washing_rate.
        
        # """
        # if proba(washing_rate):
            # dispersal_unit.disable()      

class StubGrowthControl(object):

    def __init__(self):
        pass
        
    def control(self, g, label='LeafElement'):
        """ Limit lesion growth to the healthy surface on leaves.
        
        """
        vids = [v for v in g if g.label(v).startswith(label)]
        for v in vids:
            leaf_element = g.node(v)
            if 'lesions' in leaf_element.properties():
                leaf = [le for le in leaf_element.complex().components() if le.label.startswith(label)] # Gather leaf elements of the leaf
                # Does it work for powdery mildew too ?
                leaf_lesions = []
                for le in leaf:
                    if 'lesions' in le.properties():
                        leaf_lesions.extend(le.lesions)
                
                leaf_surface = sum([lf.surface for lf in leaf])
                leaf_healthy_surface = sum([lf.healthy_surface for lf in leaf])
                
                total_demand = sum(l.growth_demand for l in leaf_lesions)
                
                if total_demand > 0. and total_demand > leaf_healthy_surface:
                    for l in leaf_lesions:
                        growth_offer = leaf_healthy_surface * l.growth_demand / total_demand
                        l.growth_control(reduce_up_to = growth_offer)                 
        
                # Update of 'healthy_surface' by leaf element:
                for le in leaf:
                    if 'lesions' in le.properties():
                        lesions_surface = sum([l.surface for l in le.lesions])
                        hs = round(le.surface - lesions_surface, 5)
                        le.healthy_surface = hs
            
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
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            for les in n.lesions:
                if les.status == state:
                    count +=1
    return count
    
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
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    plot_DU(g)
    return g
    
def test_infect():
    """ Check if 'infect' from 'protocol.py' leads to infection by dispersal units on the MTG.

    """
    g = adel_mtg2()
    initiate_g(g)
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
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
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    dt = 1
    nb_steps = 750
    nb_les_inc = numpy.array([0. for i in range(nb_steps)])
    # nb_les_chlo = numpy.array([0. for i in range(nb_steps)])
    # nb_les_nec = numpy.array([0. for i in range(nb_steps)])
    nb_les_spo = numpy.array([0. for i in range(nb_steps)])
    # nb_les_empty = numpy.array([0. for i in range(nb_steps)])
    for i in range(nb_steps):
        print('time step %d' % i)
        
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
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(10)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    dispersor = StubDispersal()
    
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
            disperse(g, scene, dispersor, "Septoria")

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
    g = adel_mtg2()
    initiate_g(g)
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    washor = StubWashing()
    
    dt = 1
    nb_steps = 100
    nb_DU = numpy.array([0. for i in range(nb_steps)])
    for i in range(nb_steps):
    
        g = update_climate(g)
        if i>2 and i%5 == 0 or (i-1)%5 == 0:
                global_rain_intensity = 4.
                rain_interception(g, rain_intensity = global_rain_intensity*0.75)      
        else:
            global_rain_intensity = 0.
        
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
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(1000)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    controler = StubGrowthControl()

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
            # count_surf += leaf.healthy_surface
            if 'lesions' in leaf.properties():
                count_surf += sum([l.surface for l in leaf.lesions])
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
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Call the models that will be used during the simulation :
    controler = StubGrowthControl()
    washor = StubWashing()
    dispersor = StubDispersal()
   
    # Prepare the simulation loop
    dt = 1
    start_date = datetime(2000, 10, 1)
    end_date = datetime(2001, 7, 31)
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        
        # Update climate on leaf elements : 
        update_on_leaves(weather_data, date, g)     
        
        # Run a disease cycle
        #grow(g)
        infect(g, dt)        
        update(g,dt)
        growth_control(g, controler)
        
        scene = plot3d(g)
        disperse(g, scene, dispersor, "Septoria")
        
        global_rain_intensity = weather_data.Pluie[date]
        wash(g, washor, global_rain_intensity)
    
        # displayer = DisplayLesions()
        # displayer.print_all_lesions(g)
    
    # Display results:
    plot_lesions(g)
    return g
    
