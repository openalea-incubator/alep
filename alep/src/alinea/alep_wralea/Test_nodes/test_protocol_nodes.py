#TODO : DO NOT DUPLICATE CODE
""" Test the protocol of communication between adel and alep """

# Imports #########################################################################

import random
import numpy
import csv

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

class ReadWeather(object):
    """ Class that reads weather data.
    
    """
    def __init__(self):
        """ Initialize lists of weather data.
        
        """
        self.date = []
        self.hhmm = []
        self.PAR = []
        self.temp = []
        self.RH = []
        self.wind_speed = []
        self.rain = []
        self.wetness = []
    
    def read_weather_data(self):
        """ Open a txt file and save required weather data line after line.
        
        Reading is rough and not flexible for this example.
        Filename could be an input argument if all weather files have the same structure.
        """
        from datetime import datetime
        filename = 'meteo01.txt'
        with open(filename) as f:
            reader = csv.reader(f, delimiter='\t')
            reader.next()
            ind_line = 0
            for line in reader:
                year = int(line[0])
                julian_day = int(line[1])
                hour = ind_line%24
                t = datetime.fromordinal(datetime(year, 1, 1).toordinal() + julian_day - 1)
                t = t.replace(hour = hour)
                self.date.append(t)
                
                self.hhmm.append(float(line[3]))
                self.PAR.append(float(line[4]))
                self.temp.append(float(line[5]))
                self.RH.append(float(line[6]))
                self.wind_speed.append(float(line[7]))
                self.rain.append(float(line[8]))
                if float(line[8]) > 0. or (float(line[4]) < 644 and float(line[6]) > 85):
                    self.wetness.append(True)
                else:
                    self.wetness.append(False)
                
                ind_line +=1
        
        return self

def update_on_leaves(data, time_step, g, label = 'LeafElement'):
        """ Read weather data for a step of simulation and apply it to each leaf.
        
        """        
        vids = [n for n in g if g.label(n).startswith(label)]
        for v in vids : 
            n = g.node(v)
            
            n.temp = data.temp[time_step]
            n.rain_intensity = data.rain[time_step]
            n.relative_humidity = data.RH[time_step]
            n.wetness = data.wetness[time_step]

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

def rain_interception(g, rain_intensity = None, label = 'LeafElement'):
    """ simulate an environmental program with rain """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.rain_intensity = rain_intensity

# Stub models for disease #########################################################

def generate_stock_DU(fungus = septoria(), nb_Spores=random.randint(1,100), nb_DU = 100):
    """ Generate a stock of DU as a list of DU 
    
    """
    stock_DU = [SeptoriaDU(fungus = fungus, nbSpores=nb_Spores, nature='emitted') for i in range(nb_DU)]
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
    
    def wash(self, dispersal_unit, washing_rate):
        """ On the given LeafElement, inactive the DU as a function of the washing_rate.
        
        """
        if proba(washing_rate):
            dispersal_unit.inactive()      
            
# Display #########################################################################

def temp_plot3D(g):
    """ plot g """
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

def plot_DU_deposited(g):
    """ plot the plant with elements carrying deposited dispersal units in yellow """
    green = (0,180,0)
    yellow = (247, 220, 17)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        ind_deposited = 0.
        if 'dispersal_units' in n.properties():
            for du in n.dispersal_units:
                if du.nature == 'deposited':
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
    g = initiate_g(g)
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    plot_DU(g)
    return g
    
def test_infect():
    """ Check if 'infect' from 'protocol.py' leads to infection by dispersal units on the MTG.

    """
    g = adel_mtg2()
    g = initiate_g(g)
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
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
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    dt = 1
    nb_steps = 1000
    for i in range(nb_steps):
        print('time step %d' % i)
        
        update_climate(g)
        if i%100 == 0:
            global_rain_intensity = 4.
            rain_interception(g, rain_intensity = global_rain_intensity*0.75)  
        else:
            global_rain_intensity = 0.
            
        #grow(g)
        infect(g, dt)        
        update(g,dt)
    
    displayer = DisplayLesions()
    displayer.print_all_lesions(g)

    plot_lesions(g)
    return g

def test_disperse():
    """ Check if 'disperse' from 'protocol.py' disperse new dispersal units on the MTG.

    """
    g = adel_mtg2()
    g = initiate_g(g)
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
            global_rain_intensity = 4.
            rain_interception(g, rain_intensity = global_rain_intensity*0.75)  
        else:
            global_rain_intensity = 0.
        
        # grow(g)
        infect(g, dt)
        update(g,dt)
        if global_rain_intensity != 0.:
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
    g = initiate_g(g)
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    washor = StubWashing()
    
    dt = 1
    nb_steps = 100
    for i in range(nb_steps):
    
        update_climate(g)
        if i>2 and i%5 == 0 or (i-1)%5 == 0:
                global_rain_intensity = 4.
                rain_interception(g, rain_intensity = global_rain_intensity*0.75)      
        else:
            global_rain_intensity = 0.
        
        wash(g, washor, global_rain_intensity)
        
        if global_rain_intensity != 0:
                   
            dispersal_units = g.property('dispersal_units')
            nbDU = 0.
            for vid, dlist in dispersal_units.iteritems():
                for d in dlist:
                    if d.active:
                        nbDU += 1
        
        if i>2 and (i)%5 == 0:
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

def test_growth_control():
    g = adel_mtg()
    g = initiate_g(g)
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(1000)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
   
    dt = 1
    nb_steps = 1000
    for i in range(nb_steps):
        print('time step %d' % i)
        
        update_climate(g)
        if i%100 == 0:
            global_rain_intensity = 4.
            rain_interception(g, rain_intensity = global_rain_intensity*0.75)  
        else:
            global_rain_intensity = 0.
            
        #grow(g)
        infect(g, dt)        
        update(g,dt)
        growth_control(g)
        
        vids = [v for v in g if g.label(v).startswith("LeafElement")]
        for v in vids:
            leaf = g.node(v)
            if 'lesions' in leaf.properties():
                print('leaf element %d - ' % v + 'Healthy surface : %f'  % leaf.healthy_surface) 
        
    return g
        
def test_simul_with_weather():
    """ Check if 'update' from 'protocol.py' provokes the growth of a lesion instantiated on the MTG.

    """
    from datetime import datetime
    from dateutil import rrule
    
    g = adel_mtg2()
    g = initiate_g(g)
    weather = ReadWeather()
    weather_data = weather.read_weather_data()
    stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    washor = StubWashing()
   
    dt = 1
    start_date = datetime(2000, 10, 1)
    end_date = datetime(2001, 7, 31)
    for date in rrule.rrule(rrule.HOURLY, dtstart = start_date, until = end_date):
        print(date)
        
        for t in range(len(weather_data.date)):
            if weather_data.date[t] == date:
                time_step = t
        update_on_leaves(weather_data, time_step, g)     
        
        #grow(g)
        infect(g, dt)        
        update(g,dt)
        growth_control(g)
        
        scene = plot3d(g)
        dispersor = StubDispersal()
        disperse(g, scene, dispersor, "Septoria")
        
        wash(g, washor, weather_data)
    
        displayer = DisplayLesions()
        displayer.print_all_lesions(g)

    plot_lesions(g)
    return g
    
