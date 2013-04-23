""" Example of use of the model of continuous lesions of septoria with real weather data.

"""
# Imports #########################################################################
import random as rd
import numpy as np
import pandas
from pylab import *
import matplotlib.pyplot as plt

from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe

from alinea.alep.protocol import *
from alinea.alep.septoria import *
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.dispersal import RandomDispersal
from alinea.alep.washing import RapillyWashing

from datetime import datetime
from datetime import timedelta
import time

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
         'Einc':[0],'El':[0],'Ev':[0],'Esen':[0],'Ed': [0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=leaves_db(), leaf_sectors=1)
    g=mtg_interpreter(g)
    return g
    
def adel_mtg2(nb_sect=1):
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaves_db(),stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
    g=mtg_interpreter(g)
    return g

def set_initial_properties_g(g, 
                             surface_leaf_element=5.,
                             position_senescence=None,
                             label = 'LeafElement'):
    """ Give initial values for plant properties of each LeafElement. 
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    surface: float
        Initial surface of each leaf element
    position_senescence: float
        Position of senescence on blade axis
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
        n.position_senescence = None
        
    return g
    
# Climate #########################################################################
def read_weather(filename, call_update="step_by_step"):
    """ Use pandas to save data from a .txt file into a dataframe.
    
    call_update="step_by_step" or "dispersal"
    
    """
    raw_data = pandas.read_table(filename)
    
    # Addition of wetness variable
    wet = dict(wetness=[])
    for i_line in range(len(raw_data.An)):
        wet['wetness'].append(compute_wetness(raw_data.Pluie[i_line], 
                              raw_data.HR[i_line],
                              raw_data.PAR[i_line]))
    wetness = pandas.DataFrame(wet)
    weather_data = raw_data.join(wetness)
    
    # Compute dispersal events and keep only rain max at dispersal occurence
    dispersal = []
    max_rain = 0.
    ind_max = 0.
    rain = zeros(len(raw_data.Pluie))
    if call_update == "dispersal":
        call = zeros(len(raw_data.Pluie))
    elif call_update == "step_by_step":
        call = ones(len(raw_data.Pluie))
        
    for i_line in range(len(raw_data.An)):
        dispersal.append(compute_dispersal_event(raw_data.Pluie[i_line], raw_data.HR[i_line]))
        if dispersal[i_line] == True:
            if raw_data.Pluie[i_line] > max_rain:
                max_rain = raw_data.Pluie[i_line]
                ind_max = i_line
        elif ind_max > 0.:
            rain[ind_max] = max_rain
            if call_update == "dispersal":
                call[ind_max] = True
            max_rain = 0.
            ind_max = 0.
    
    del weather_data['Pluie']
    Pluie = pandas.DataFrame(dict(Pluie=rain))
    weather_data = weather_data.join(Pluie)
    
    dispersal_event = pandas.DataFrame(dict(dispersal_event=dispersal))
    weather_data = weather_data.join(dispersal_event)

    call_update = pandas.DataFrame(dict(call_update=call))
    weather_data = weather_data.join(call_update)
    
    # Creation of date time indexes
    first_date = datetime.fromordinal(datetime(weather_data.An[0], 1, 1).toordinal() + weather_data.Jour[0] - 1)
    rng = pandas.date_range(first_date, periods = len(weather_data.An), freq ='H')
    
    weather_data.index = rng

    return weather_data

def compute_wetness(rain=0., relative_humidity=0., PAR=0.):
    """ Compute leaf wetness as in Rapilly et Jolivet as a function
        of rain or relative humidity and PAR.
    
    Parameters
    ----------
    rain: float
        Rain (in mm)
    relative_humidity: float
        Relative humidity of the air around the leaf (in %)
    PAR: float
        Photosynthetically active radiation around the leaf (in nm)
    
    Returns
    -------
    wet: True or False
        True if the leaf is wet
    """
    if rain>0. or (relative_humidity >= 85. and PAR < 644.):
        return True
    else:
        return False

def compute_dispersal_event(rain=0., relative_humidity=0.):
    """ Compute a dispersal event as a function of rain or relative humidity.
    
    Parameters
    ----------
    rain: float
        Rain (in mm)
    relative_humidity: float
        Relative humidity of the air around the leaf (in %)
    
    Returns
    -------
    dispersal_event: True or False
        True if dispersal event is occuring 
    """
    if rain>0.5 and relative_humidity >= 85.:
        return True
    else:
        return False
        
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

def update_temp_list_on_g(g, weather_data, start_date, 
                          date=None, label = 'LeafElement'):
    """ Update the property list of temperatures on the leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    weather_data: Pandas DataFrame
        Data frame with climatic data
    date: datetime
        Date in datetime format
        
    Returns
    ----------
    g: MTG
        Updated MTG representing the canopy
    """
    # Hack:
    weather_data.call_update[start_date]=True
    delta = timedelta(hours=1)
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids :
        n = g.node(v)
        t = date
        n.temp_list = [weather_data.Tair[t]]
        t-=delta
        if t!=(start_date-delta):
            while weather_data.call_update[t] == False:
                if t == (start_date-delta):
                    break
                else:
                    n.temp_list.append(weather_data.Tair[t])
                    t-=delta
    

## TEMP
def update_climate_all(g, wetness=True,
                          temp = 22.,
                          rain_intensity=0.,
                          relative_humidity=85.,
                          label = 'LeafElement'):
    """ Simulate an environmental program.
    
    All leaf elements have the same values for all variables.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
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
    g: MTG
        Updated MTG representing the canopy
    """
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.wetness = wetness
        n.temp = temp
        n.rain_intensity = rain_intensity
        n.relative_humidity = relative_humidity
    
    return g
    
# Measure variables ###############################################################
def compute_total_severity(g):
    """ Compute disease severity on the whole plant.
    
    Severity is the ratio between disease surface and total surface of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    severity: float
        Ratio between disease surface and total surface of leaves (in %)
    """
    # Initiation
    surfaces = g.property('surface')
    healthy_surfaces = g.property('healthy_surface')
    
    total_leaf_surface = sum(surfaces.values())
    total_disease_surface = total_leaf_surface - sum(healthy_surfaces.values())
    
    # # Select all the leaves
    # bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
    # for blade in bids:
        # leaf = [vid for vid in g.components(blade) if labels[vid].startswith('LeafElement')]
        # leaf_surface = sum(surfaces[lf] for lf in leaf)
        # # Compute total surface of leaves
        # total_leaf_surface += leaf_surface
        # # Compute disease surface on leaves
        # leaf_healthy_surface = sum(healthy_surfaces[lf] for lf in leaf)
        # disease_surface = max(0., leaf_surface - leaf_healthy_surface)
        # total_disease_surface += disease_surface
        
    # Compute ratio, i.e. severity
    severity = 100 * total_disease_surface / total_leaf_surface if total_leaf_surface > 0. else 0.
    
    return severity

def compute_total_necrosis(g):
    """ Compute necrosis percentage on the whole plant.
    
    Necrosis percentage ratio between necrotic (and sporulating) disease surface and total surface of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    necrosis_percentage: float
        Ratio between necrotic (and sporulating) disease surface and total surface of leaves (in %)
    """
    # Leaf
    surfaces = g.property('surface')
    total_leaf_surface = sum(surfaces.values())
    
    # Disease
    lesions = g.property('lesions')
    if lesions:
        lesions = [l for les in lesions.values() for l in les if l.status>=l.fungus.NECROTIC]
        total_necrotic_surface = sum(l.surface for l in lesions)
    else:
        total_necrotic_surface = 0.
    
    
    # total_disease_surface = total_leaf_surface - sum(healthy_surfaces.values())
    
    # Select all the leaves
    # bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
    # for blade in bids:
        # leaf = [vid for vid in g.components(blade) if labels[vid].startswith('LeafElement')]
        # leaf_surface = sum(surfaces[lf] for lf in leaf)
        # Compute disease surface on leaves
        # leaf_healthy_surface = sum(healthy_surfaces[lf] for lf in leaf)
        # disease_surface = max(0., leaf_surface - leaf_healthy_surface)
        # total_disease_surface += disease_surface
    
    # Compute ratio, i.e. severity
    necrosis_percentage = 100 * total_necrotic_surface / total_leaf_surface if total_leaf_surface > 0. else 0.
    
    return necrosis_percentage

def compute_state_one_leaf(g, status="SPORULATING"):
    """ Compute surface of lesions in chosen state on the MTG.
    
    """
    lesions = g.property('lesions')
    if lesions:
        if status=="INCUBATING":
            lesions = [l for les in lesions.values() for l in les if l.status==l.fungus.INCUBATING]
        if status=="CHLOROTIC":
            lesions = [l for les in lesions.values() for l in les if l.status==l.fungus.CHLOROTIC]
        if status=="NECROTIC":
            lesions = [l for les in lesions.values() for l in les if l.status==l.fungus.NECROTIC]
        if status=="SPORULATING":
            lesions = [l for les in lesions.values() for l in les if l.status==l.fungus.SPORULATING]
        surface = sum(l.surface for l in lesions)
    else:
        surface = 0.
    
    return surface
    
    
def count_DU(g, status='deposited'):
    """ count DU of the mtg.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    status: str
        Status of the dus to count ('deposited', 'emitted' or 'all')
        
    Returns
    -------
        Number of dus on the MTG
    
    """
    dispersal_units = g.property('dispersal_units')
    if status=='all':
        return sum(len(du) for du in dispersal_units.itervalues())
    else:
        return len([d for du_list in dispersal_units.itervalues() for d in du_list if d.status==status])
        
def count_lesions(g):
    """ count lesions of the mtg.
    
    """
    lesions = g.property('lesions')
    return sum(len(l) for l in lesions.itervalues())
    
# Septoria ########################################################################
def test_all():
    # Generate a MTG with required properties :
    g = adel_mtg2()
    # g = adel_one_leaf()
    set_initial_properties_g(g)
    
    # Deposit first dispersal units on the MTG :
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted') for i in range(1)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

    # Prepare the simulation loop
    dt = 1
    nb_steps = 3500
    severity = np.zeros(nb_steps)
    nb_DUs = np.zeros(nb_steps)
    nb_lesions = np.zeros(nb_steps)
    for i in range(0,nb_steps,dt):
        print('time step %d' % i)
        dt_pluie = 1
        # Update climate and force rain occurences       
        if i>0 and i%dt_pluie == 0:
            global_rain_intensity = 4.
        else:
            global_rain_intensity = 0.
        update_climate_all(g, wetness=True, temp=22., rain_intensity = global_rain_intensity*0.75)
        
        # grow(g)
        infect(g, dt)
        update(g, dt_pluie, growth_control_model=controler)

        if global_rain_intensity != 0.:
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria")
        
        wash(g, washor, global_rain_intensity, DU_status='deposited')
        
        # Measure severity
        severity[i]=compute_total_severity(g)
        
        nb_DUs[i] = count_DU(g, 'deposited')
        nb_lesions[i] = count_lesions(g)
    
    # Display results
    fig = plt.figure()
    
    ax1 = fig.add_subplot(3,1,1)
    plot(range(0,nb_steps,dt), severity)
    ylim([0, 105])
    ylabel('Total disease severity')
    
    ax2 = fig.add_subplot(3,1,2)
    plt.bar(range(0,nb_steps,dt), nb_DUs, color='r')
    ylabel('Number of DUs')

    ax3 = fig.add_subplot(3,1,3)
    plot(range(0,nb_steps,dt), nb_lesions, color='r')
    ylabel('Number of lesions')
    xlabel('Simulation time step')
    
    fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g

def simulate_severity():
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.txt'
    weather_data = read_weather(filename, "dispersal")
    
    # Generate a MTG with required properties :
    # g = adel_one_leaf()
    g = adel_mtg2()
    set_initial_properties_g(g)
    
    # Deposit first dispersal units on the MTG :
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

    # Prepare the simulation loop
    dt = 1
    start_date = datetime(2000, 10, 1)
    end_date = datetime(2001, 7, 31)
    # end_date = datetime(2000, 11, 30)
    severity = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    nb_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    nb_lesions = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try:
            ind+=1.
        except:
            ind = 0.
        # Update climate on leaf elements : 
        update_on_leaves(weather_data, date, g)     
        
        # Run a disease cycle
        #grow(g)
        infect(g, dt)
        if weather_data.call_update[date]:
            # print(date)
            update_temp_list_on_g(g, weather_data, start_date, date)
            update(g, dt, growth_control_model=controler)
        # control_growth(g, controler)
               
        global_rain_intensity = weather_data.Pluie[date]
        if global_rain_intensity != 0. :
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria") 
        wash(g, washor, global_rain_intensity)
    
        # Measure severity
        severity[ind]=compute_total_severity(g)
        
        nb_DUs[ind] = count_DU(g, 'deposited')
        nb_lesions[ind] = count_lesions(g)
    
    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(3,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), severity)
    ylim([0, 105])
    ylabel('Total disease severity')
    
    ax2 = fig.add_subplot(3,1,2)
    plt.bar(range(nb_steps), nb_DUs, color='r')
    ylabel('Number of DUs')
    xlim([0, nb_steps])

    ax3 = fig.add_subplot(3,1,3)
    plot(range(nb_steps), nb_lesions, color='r')
    ylabel('Number of lesions')
    xlim([0, nb_steps])
    xlabel('Simulation time step')
    
    fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g
 
def simulate_necrosis():
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.txt'
    weather_data = read_weather(filename, "dispersal")
    
    # Generate a MTG with required properties :
    # g = adel_one_leaf()
    g = adel_mtg2()
    set_initial_properties_g(g)
    
    # Deposit first dispersal units on the MTG :
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

    # Prepare the simulation loop
    dt = 1
    start_date = datetime(2000, 10, 1)
    end_date = datetime(2001, 7, 31)
    # end_date = datetime(2000, 11, 30)
    necrosis = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try:
            ind+=1.
        except:
            ind = 0.
        # Update climate on leaf elements : 
        update_on_leaves(weather_data, date, g)     
        
        # Run a disease cycle
        #grow(g)
        infect(g, dt)
        if weather_data.call_update[date]:
            # print(date)
            update_temp_list_on_g(g, weather_data, start_date, date)
            update(g, dt, growth_control_model=controler)
        # control_growth(g, controler)
               
        global_rain_intensity = weather_data.Pluie[date]
        if global_rain_intensity != 0. :
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria") 
        wash(g, washor, global_rain_intensity)
    
        # Measure severity
        necrosis[ind]=compute_total_necrosis(g)
    
    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    # fig = plt.figure()
    
    # ax1 = fig.add_subplot(3,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), necrosis)
    ylim([0, 105])
    ylabel('Total necrosis percentage')
    
    # ax2 = fig.add_subplot(3,1,2)
    # plt.bar(range(nb_steps), nb_DUs, color='r')
    # ylabel('Number of DUs')
    # xlim([0, nb_steps])

    # ax3 = fig.add_subplot(3,1,3)
    # plot(range(nb_steps), nb_lesions, color='r')
    # ylabel('Number of lesions')
    # xlim([0, nb_steps])
    # xlabel('Simulation time step')
    
    # fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g
    
def simulate_states_one_leaf():
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.txt'
    weather_data = read_weather(filename, "step_by_step")
    
    # Generate a MTG with required properties :
    g = adel_one_leaf()
    set_initial_properties_g(g)
    
    # Deposit first dispersal units on the MTG :
    fungus = septoria(); SeptoriaDU.fungus = fungus
    stock = [SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted') for i in range(100)]
    inoculator = RandomInoculation()
    initiate(g, stock, inoculator)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

    # Prepare the simulation loop
    dt = 1
    start_date = datetime(2000, 10, 1)
    end_date = datetime(2000, 11, 30)
    incubating = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    chlorotic = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    necrotic = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    sporulating = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try:
            ind+=1.
        except:
            ind = 0.
        # Update climate on leaf elements : 
        update_on_leaves(weather_data, date, g)     
        
        # Run a disease cycle
        #grow(g)
        infect(g, dt)
        if weather_data.call_update[date]:
            # print(date)
            update_temp_list_on_g(g, weather_data, start_date, date)
            update(g, dt, growth_control_model=controler)
        # control_growth(g, controler)
               
        global_rain_intensity = weather_data.Pluie[date]
        if global_rain_intensity != 0. :
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria") 
        wash(g, washor, global_rain_intensity)
    
        # Measure severity
        incubating[ind] = compute_state_one_leaf(g, status="INCUBATING")
        chlorotic[ind] = compute_state_one_leaf(g, status="CHLOROTIC")
        necrotic[ind] = compute_state_one_leaf(g, status="NECROTIC")
        sporulating[ind]= compute_state_one_leaf(g, status="SPORULATING")
    
    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    y_max = max([max(incubating), max(chlorotic), max(necrotic), max(sporulating)])+1
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(4,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), incubating)
    ylim([0, y_max])
    ylabel('Incubating surface')
    
    ax2 = fig.add_subplot(4,1,2)
    plot(pandas.date_range(start_date, end_date, freq='H'), chlorotic)
    ylabel('Chlorotic surface')
    ylim([0, y_max])
    # xlim([0, nb_steps])

    ax3 = fig.add_subplot(4,1,3)
    plot(pandas.date_range(start_date, end_date, freq='H'), necrotic)
    ylabel('Necrotic surface')
    ylim([0, y_max])
    # xlim([0, nb_steps])

    ax4 = fig.add_subplot(4,1,4)
    plot(pandas.date_range(start_date, end_date, freq='H'), sporulating)
    ylabel('Sporulating surface')
    ylim([0, y_max])
    # xlim([0, nb_steps])
    xlabel('Simulation time step')
    
    fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g