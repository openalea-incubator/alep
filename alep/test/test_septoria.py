""" Test basic responses of the wheat septoria model """

# Imports #########################################################################
import random as rd

from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe

from alinea.alep.protocol import *
from alinea.alep.septoria import *
from alinea.alep.inoculation import RandomInoculation

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
    
def adel_one_leaf():
    """ create a very simple adel mtg """
    d = {'plant':[1],'axe_id':['MS'],'ms_insertion':[0],'numphy':[1], 
         'Laz': [0], 'Ll' :[3], 'Lv' :[3] , 'Lsen':[0], 'L_shape':[3], 'Lw_shape':[.3], 'Linc':[0],
         'Einc':[0],'El':[1],'Ev':[1],'Esen':[0],'Ed': [0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=leaves_db(), leaf_sectors=1)
    g=mtg_interpreter(g)
    return g
    
def adel_mtg2(nb_sect=1):
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaves_db(),stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
    g=mtg_interpreter(g)
    return g
    
# Climate #########################################################################
def set_initial_properties_g(g, surface_leaf_element=5., label = 'LeafElement'):
    """ Give initial values for plant properties of each LeafElement. 
    
    Parameters
    ----------
    surface: float
        Initial surface of each leaf element
    label: str
        Label of the part of the MTG concerned by the calculation
    
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy
    """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.surface = surface_leaf_element
        n.healthy_surface = surface_leaf_element # TODO : Manage properly
        
    return g
    
def update_climate_all(g, wetness=True,
                          temp = 18.,
                          rain_intensity=0.,
                          relative_humidity=85.,
                          label = 'LeafElement'):
    """ Simulate an environmental program.
    
    All leaf elements have the same values for all variables.
    
    Parameters
    ----------
    wetness: bool
        True if the leaf element is wet, False otherwise
    temp: float
        Temperature of the leaf element (degrees celsius)
    rain_intensity : float
        Rain intensity on the leaf element (mm/h)
    relative_humidity : float
        Relative humidity on the leaf element (percent)
    label: str
        Label of the part of the MTG concerned by the calculation ('LeafElement')
     
    Returns
    -------
    None
        Update directly the MTG
    """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.wetness = wetness
        n.temp = temp
        n.rain_intensity = rain_intensity
        n.relative_humidity = relative_humidity
    
    return g

# Tests ###########################################################################
def test_initiate():
    """ Check if 'initiate' from 'protocol.py' deposits dispersal units on the MTG.
    
    Generate a wheat MTG and distribute dispersal units randomly on leaf elements.
    Check that all the stock of DU is distributed in properties of leaf elements.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    set_initial_properties_g(g, surface_leaf_element=5.)
    
    # Generate a stock of septoria dispersal units
    fungus = septoria()
    SeptoriaDU.fungus = fungus
    nb_dus_in_stock = 100
    stock = [SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted') for i in range(nb_dus_in_stock)]
    
    # Call the protocol of initiation with a model distributing the DUs randomly
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Check that all the DUs that were in the stock are now instantiated 
    # on leaf elements :
    dispersal_units = g.property('dispersal_units')
    nb_dus_on_leaves = sum(len(dus) for dus in dispersal_units.itervalues())
    assert nb_dus_on_leaves == nb_dus_in_stock
    
def test_infect():
    """ Check if 'infect' from 'protocol.py' leads to infection by dispersal units 
        on the MTG.
    
    Generate a wheat MTG and distribute dispersal units randomly on leaf elements.
    Run a short loop to compute infection. Check that all the stock of DU caused
    appearance of a lesion on leaf elements.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    set_initial_properties_g(g, surface_leaf_element=5.)
    
    # Generate a stock of septoria dispersal units
    fungus = septoria()
    SeptoriaDU.fungus = fungus
    nb_dus_in_stock = 100
    stock = [SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted') for i in range(nb_dus_in_stock)]
    
    # Call the protocol of initiation with a model distributing the DUs randomly
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
        
    # Loop of simulation
    dt = 1
    nb_steps = 11
    for i in range(0,nb_steps,dt):
        # Offer good conditions for at least 10h:
        update_climate_all(g, wetness=True, temp=20.)
        infect(g, dt)
    
    # Check that all the DUs that were in the stock are now lesions on leaf elements :
    lesions = g.property('lesions')
    nb_lesions_on_leaves = sum(len(l) for l in lesions.itervalues())
    assert nb_lesions_on_leaves == nb_dus_in_stock
    
def test_update():
    """ Check if 'update' from 'protocol.py' provokes the growth of a lesion
        instantiated on the MTG.
    
    Generate a wheat MTG and deposit 1 dispersal unit on a leaf element.
    Run a loop to compute infection and update. Check that all the stages
    of a lesion of septoria have been reached eventually.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    set_initial_properties_g(g, surface_leaf_element=5.)
    
    # Generate a stock of septoria dispersal units
    fungus = septoria()
    SeptoriaDU.fungus = fungus
    nb_dus_in_stock = 1
    stock = [SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted') for i in range(nb_dus_in_stock)]
    
    # Call the protocol of initiation with a model distributing the DUs randomly
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
        
    # Loop of simulation
    dt = 1
    nb_steps = 500
    di = 0. # delay before infection
    for i in range(0,nb_steps,dt):
        # After infection, the lesion 'age_dday' will be added 1 DD by time step
        # Note : base temperature for septoria = -2 degrees celsius
        update_climate_all(g, wetness=True, temp=22.)
        infect(g, dt)
        update(g, dt)

        # Check that the lesion is in the right status
        lesion = g.property('lesions')
        if lesion and di==0.:
            di = i 
            assert sum(len(l) for l in lesion.itervalues()) == 1
            l = lesion.values()[0][0]
            print(l.status)
            if (220.+di) <= i < (330.+di):
                assert l.status == 1
            elif (330.+di) <= i < (350.+di):
                assert l.status == 2
            elif i >= (350.+di):
                assert l.status == 3
                
def test_disperse():
    """ Check if 'disperse' from 'protocol.py' disperse new 
        dispersal units on the MTG.
    
    Generate a wheat MTG and distribute dispersal units randomly on leaf elements.
    Run a loop to compute infection and update. Create artificial rains in the loop
    on a regular basis. Check that the number of lesions on the MTG increases.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    set_initial_properties_g(g, surface_leaf_element=5.)
    
    # Generate a stock of septoria dispersal units
    fungus = septoria()
    SeptoriaDU.fungus = fungus
    nb_dus_in_stock = 1
    stock = [SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted') for i in range(nb_dus_in_stock)]
    
    # Call the protocol of initiation with a model distributing the DUs randomly
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
        
    # Loop of simulation
    dt = 1
    nb_steps = 700
    for i in range(0,nb_steps,dt):
        update_climate_all(g, wetness=True, temp=22.)
        infect(g, dt)
        update(g, dt)
        
        # Force rain occurences
        if i>400 & i%100 == 0:
            global_rain_intensity = 4.
            update_climate_all(g, rain_intensity = global_rain_intensity*0.75)  
            
            