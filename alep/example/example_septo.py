""" Example of use of the model of continuous lesions of septoria with real weather data.

"""
# Imports #########################################################################
import random as rd
import numpy as np
import pandas
from pylab import *
import matplotlib.pyplot as plt
import string

from alinea.alep.wheat_examples import adel_mtg, adel_mtg2, adel_mtg3, adel_one_leaf
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

from alinea.alep.protocol import *
from alinea.alep.septoria import *
from alinea.alep.microclimate_leaf import *
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.du_position_checker import BiotrophDUProbaModel
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.dispersal import RandomDispersal
from alinea.alep.washing import RapillyWashing

from datetime import datetime, timedelta
import time

# Plant ###########################################################################
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
def read_weather(filename, update_events="step_by_step"):
    """ Use pandas to save data from a .txt file into a dataframe.
    
    update_events="step_by_step" or "dispersal"
    
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
    rain_counter = 0.
    rain = zeros(len(raw_data.Pluie))
    rain_duration = zeros(len(raw_data.Pluie))
    if update_events == "dispersal":
        call = zeros(len(raw_data.Pluie))
    elif update_events == "step_by_step":
        call = ones(len(raw_data.Pluie))
    else:
        raise Exception("Wrong str for update_events")
        
    for i_line in range(len(raw_data.An)):
        dispersal.append(compute_dispersal_event(raw_data.Pluie[i_line], raw_data.HR[i_line]))
        if dispersal[i_line] == True:
            rain_counter += 1.
            if raw_data.Pluie[i_line] > max_rain:
                max_rain = raw_data.Pluie[i_line]
                ind_max = i_line
        elif ind_max > 0.:
            rain[ind_max] = max_rain
            rain_duration[ind_max] = rain_counter
            if update_events == "dispersal":
                call[ind_max] = True
            max_rain = 0.
            ind_max = 0.
            rain_counter = 0.
    
    del weather_data['Pluie']
    Pluie = pandas.DataFrame(dict(Pluie=rain))
    weather_data = weather_data.join(Pluie)
    
    rain_duration = pandas.DataFrame(dict(rain_duration=rain_duration))
    weather_data = weather_data.join(rain_duration)
    
    dispersal_event = pandas.DataFrame(dict(dispersal_event=dispersal))
    weather_data = weather_data.join(dispersal_event)

    update_events = pandas.DataFrame(dict(update_events=call))
    weather_data = weather_data.join(update_events)
    
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
    # # # Use CaribuMicroclimModel
    # # mean_globalclimate = {'PAR':data.PAR[date],
                          # # 'Pluie':data.Pluie[date],
                          # # 'Tair':data.Tair[date],
                          # # 'HR':data.HR[date],
                          # # 'Vent':data.Vent[date]}
                          
    # # microclim_model = CaribuMicroclimModel()
    # # microclimate = microclim_model.microclim(mean_globalclimate, scene=plot3d(g))
    
    # # # Distribution on leaf elements
    # # for vid, variables in microclimate.iteritems():
        # # leaf = g.node(vid)
        # # leaf.PAR = variables['PAR']
        # # leaf.rain_intensity = variables['rain']
        # # leaf.rain_duration = data.rain_duration[date]
        # # leaf.temp = variables['Tleaf']
        # # leaf.relative_humidity = variables['humidity']
        # # leaf.wetness = compute_wetness(relative_humidity=variables['humidity'],
                                       # # PAR=variables['PAR'])
        
    # Equals weather data on the entire plant for the rest
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.temp = data.Tair[date]
        n.PAR = data.PAR[date]
        n.rain_intensity = data.Pluie[date]
        n.rain_duration = data.rain_duration[date]
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
    # Hack for first line:
    weather_data.update_events[start_date]=True
    delta = timedelta(hours=1)
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids :
        n = g.node(v)
        t = date
        n.temp_list = [weather_data.Tair[t]]
        t-=delta
        if t!=(start_date-delta):
            while weather_data.update_events[t] == False:
                if t == (start_date-delta):
                    break
                else:
                    n.temp_list.append(weather_data.Tair[t])
                    t-=delta
    
def update_climate_all(g, wetness=True,
                          temp = 22.,
                          rain_intensity=0.,
                          rain_duration=0.,
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
    rain_intensity: float
        Rain intensity on the leaf element (mm/h)
    rain_duration: float
        Rain duration (in hours)
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
        n.rain_duration = rain_duration
        n.relative_humidity = relative_humidity
    
    return g

def compute_hourly_delta_ddays(temp=0., basis_for_dday=-2):
    """ Compute delta degree days in the hour.
    
    Parameters
    ----------
    dt: int
        Time step of the simulation (in hours)
    temp: float
        Temperature from weather data
    basis_for_dday: float
        Basis temperature for degree days calculation
    """
   
    return max(0,(temp - basis_for_dday)/24.)
    
# Fungus ##########################################################################
def distribute_dispersal_units(g, nb_dus=1, model="SeptoriaExchangingRings"):
    """ Distribute new dispersal units on g. 
    
    Call the method 'initiate' from the protocol with dispersal units.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    nb_dus: int
        Number of dispersal units to put on g
        
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy
    """
    fungus = septoria(model=model)
    SeptoriaDU.fungus = fungus
    dispersal_units = ([SeptoriaDU(nb_spores=rd.randint(1,100), status='emitted')
                        for i in range(nb_dus)])

    inoculator = RandomInoculation()
    initiate(g, dispersal_units, inoculator)
    
    return g
    
def distribute_lesions(g, nb_lesions=1, model="SeptoriaExchangingRings"):
    """ Distribute new lesions on g. 
    
    Call the method 'initiate' from the protocol with lesions.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    nb_lesions: int
        Number of lesions to put on g
    model: str
        Type of model of septoria lesion
        
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy
    """
    fungus = septoria(model=model)
    models = ({"SeptoriaExchangingRings":SeptoriaExchangingRings,
                    "SeptoriaWithRings":SeptoriaWithRings, 
                    "ContinuousSeptoria":ContinuousSeptoria})
    if model in models:
        models[model].fungus = fungus
        lesions = [models[model](nb_spores=rd.randint(1,100)) for i in range(nb_lesions)]

    inoculator = RandomInoculation()
    initiate(g, lesions, inoculator)
    
    return g

def execute_one_step(g, weather_data, start_date, date, dt,
                    position_checker=BiotrophDUProbaModel(),
                    controler=NoPriorityGrowthControl(),
                    dispersor=RandomDispersal(), 
                    washor=RapillyWashing()):
    """ Execute one time step of complete simulation for septoria.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    weather_data: pandas DataFrame
        Data frame of meteorological data (see 'read weather')
    start_date: format datetime
        Start date of the simulation
    date: format datetime
        Date of the simulation
    dt: int
        Time step of the simulation (in hours)
    position_checker: model
        Model that disable the DU if it falls on an existing lesion or senescent tissue
        Requires a method 'check position' (see doc)
    controler: model
        Model with rules of competition between the lesions
    dispersor: model
        Model used to position each DU in stock on g
    washor: model
        Model used to wash the DUs out of the leaf
       
    Returns
    -------
        Number of dus on the MTG
    """
    # Update climate on leaf elements : 
    update_on_leaves(weather_data, date, g)     
    
    # Run a disease cycle
    # grow(g)
    infect(g, dt, position_checker)
    if weather_data.update_events[date]:
        update_temp_list_on_g(g, weather_data, start_date, date)
        update(g, dt, growth_control_model=controler)
           
    global_rain_intensity = weather_data.Pluie[date]
    if global_rain_intensity != 0. :
        scene = plot3d(g)
        seed(1)
        disperse(g, scene, dispersor, "Septoria")
        seed(1)
        wash(g, washor, global_rain_intensity)
    
    return g

# Measure variables ###############################################################
def int2str(integer):
    """ Convert an integer to a string.
    """
    return "%d" % integer

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

    # Compute ratio, i.e. severity
    necrosis_percentage = 100 * total_necrotic_surface / total_leaf_surface if total_leaf_surface > 0. else 0.
    
    return necrosis_percentage

class LeafInspector:
    def __init__(self, g, leaf_number=1):
        """ Find the ids of the leaf elements on the chosen blade of the main stem.
        
        First leaf is the upper one. This method returns the id of the 'blade'
        carrying 'LeafElements'.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        leaf_number: int
            Number of the chosen leaf
        """
        labels = g.property('label')
        surfaces = g.property('surface')
        mets = [n for n in g if g.label(n).startswith('metamer') and g.order(n)==0]
        bids = [v for v,l in labels.iteritems() if l.startswith('blade')]
        blade_id = bids[len(mets)-leaf_number]
        # Find leaf elements on the blade
        self.leaf_elements = [id for id in g.components(blade_id) if labels[id].startswith('LeafElement')]
        # Initialize ratios (surfaces in state compared to green surface on leaf)
        self.ratio_inc = []
        self.ratio_chlo = []
        self.ratio_nec = []
        self.ratio_spo = []
        # Initialize total severity
        self.severity = []
        # Initialize necrosis percentage
        self.necrosis_percentage = []
    
    def compute_states(self, g):
        """ Compute surface of lesions in chosen state on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        surfaces = g.property('surface')
        leaf_elements = self.leaf_elements
        
        # Compute total surface on the blade
        total_surface = sum(surfaces[id] for id in leaf_elements)
        
        # Compute disease surface
        lesions = g.property('lesions')
        les=[]
        for id in leaf_elements:
            if id in lesions.keys():
                les += lesions[id]

        if les:
            for l in les:
                l.compute_all_surfaces()
            
            surface_inc = sum(l.surface_inc for l in les)
            surface_chlo = sum(l.surface_chlo for l in les)
            surface_nec = sum(l.surface_nec for l in les)
            surface_spo = sum(l.surface_spo for l in les)
        else:
            surface_inc = 0.
            surface_chlo = 0.
            surface_nec = 0.
            surface_spo = 0.

        self.ratio_inc.append(100 * surface_inc / total_surface if total_surface>0. else 0.)
        self.ratio_chlo.append(100 * surface_chlo / total_surface if total_surface>0. else 0.)
        self.ratio_nec.append(100 * surface_nec / total_surface if total_surface>0. else 0.)
        self.ratio_spo.append(100 * surface_spo / total_surface if total_surface>0. else 0.)
    
    def compute_severity(self, g):
        """ Compute severity on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        surfaces = g.property('surface')
        leaf_elements = self.leaf_elements
        
        # Compute green surface on the blade
        total_surface = sum(surfaces[id] for id in leaf_elements)
        # Compute disease surface
        lesions = g.property('lesions')
        disease_surface = sum([l.surface for id in leaf_elements for l in lesions[id]])

        self.severity.append(100 * disease_surface / total_surface if total_surface>0. else 0.)
    
    def count_DU(self, g):
        """ count DU of the leaf.
   
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
           
        Returns
        -------
            Number of dus on the leaf
        """
        dispersal_units = g.property('dispersal_units')
        leaf_elements = self.leaf_elements
        return sum(1 for id in leaf_elements for d in dispersal_units[id] if d.is_active)
        
    def count_DU_on_green(self, g, nb_unwashed):
        """ Count DU of the leaf that are on green tissue.
        
        Same calculation as in BiotrophProbaModel. 
        Might not be the actual number of DUs on green tissue.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        nb_unwashed: int
            Number of DUs staying after washing
           
        Returns
        -------
            Number of dus on green tissue on the leaf
        """
        leaf_elements = self.leaf_elements
        healthy = self.compute_healthy_surface(g)
        surface = sum(g.node(vid).surface for vid in leaf_elements)
        ratio = healthy / surface if (surface>0. and healthy>0.) else 0.
        
        return round(ratio * nb_unwashed)
    
    def compute_healthy_surface(self, g):
        """ Compute healthy surface on the leaf.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
           
        Returns
        -------
            Leaf healthy surface
        """
        leaf_elements = self.leaf_elements
        return sum(g.node(vid).healthy_surface for vid in leaf_elements)
    
def count_DU(g):
    """ count DU of the mtg.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
       
    Returns
    -------
        Number of dus on the MTG
    
    """
    dispersal_units = g.property('dispersal_units')
    return sum(len(du) for du in dispersal_units.itervalues())
  
def count_lesions(g):
    """ count lesions of the mtg.
    
    """
    lesions = g.property('lesions')
    return sum(len(l) for l in lesions.itervalues())
    
# Septoria ########################################################################
def test_microclimate():
    """
    """
    # Read weather data : 
    filename = 'meteo01.txt'
    weather_data = read_weather(filename, "step_by_step")
    # Generate a MTG with required properties :
    g = adel_one_leaf()
    # g = adel_mtg2()
    # Loop of simulation
    start_date = datetime(2000, 10, 1)
    end_date = datetime(2001, 7, 31)
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        # Update climate on leaf elements : 
        update_on_leaves(weather_data, date, g)
    
def test_all(model="SeptoriaExchangingRings"):
    # Generate a MTG with required properties :
    g = adel_mtg2()
    # g = adel_one_leaf()
    set_initial_properties_g(g)
    
    # Generate a stock of septoria dispersal units and distribute it on g
    nb_dus_in_stock = 1
    distribute_dispersal_units(g, nb_dus=nb_dus_in_stock, model=model)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

    # Prepare the simulation loop
    dt = 1
    nb_steps = 350
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
        update_climate_all(g, wetness=True, temp=22., 
                            rain_intensity = global_rain_intensity*0.75, 
                            rain_duration = 2.)
        
        # grow(g)
        infect(g, dt)
        update(g, dt_pluie, growth_control_model=controler)

        if global_rain_intensity != 0.:
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria")
            wash(g, washor, global_rain_intensity, DU_status='deposited')
        
        # Measure severity
        severity[i]=compute_total_severity(g)
        
        nb_DUs[i] = count_DU(g)
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

def compare_models():
    """
    """
    start_date = datetime(2000, 10, 1)
    end_date=datetime(2000, 11, 30)
    nb_lesions = 100
    # Run the simulation with each model
    severity1, nb_DUs1, nb_lesions1 = simulate_severity(start_date,end_date, nb_lesions,
                                                        model="SeptoriaExchangingRings")
    severity2, nb_DUs2, nb_lesions2 = simulate_severity(start_date,end_date, nb_lesions,
                                                        model="SeptoriaWithRings")
    severity3, nb_DUs3, nb_lesions3 = simulate_severity(start_date,end_date, nb_lesions,
                                                        model="ContinuousSeptoria")
    
    # Display results
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(3,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), severity1, color='b')
    plot(pandas.date_range(start_date, end_date, freq='H'), severity2, color='r')
    plot(pandas.date_range(start_date, end_date, freq='H'), severity3, color='g')
    ylim([0, 105])
    ylabel('Total disease severity')
    
    ax2 = fig.add_subplot(3,1,2)
    plt.bar(range(nb_steps), nb_DUs1, color='b')
    plt.bar(range(nb_steps), nb_DUs2, color='r')
    plt.bar(range(nb_steps), nb_DUs3, color='g')
    ylabel('Number of DUs')
    xlim([0, nb_steps])

    ax3 = fig.add_subplot(3,1,3)
    plot(range(nb_steps), nb_lesions1, color='b')
    plot(range(nb_steps), nb_lesions2, color='r')
    plot(range(nb_steps), nb_lesions3, color='g')
    ylabel('Number of lesions')
    xlim([0, nb_steps])
    xlabel('Simulation time step')
    
    fig.subplots_adjust(hspace=1)
    
    plt.show()

def compare_time_steps(model="ContinuousSeptoria"):
    """
    """
    start_date = datetime(2000, 10, 1)
    # end_date=datetime(2000, 11, 30)
    end_date=datetime(2001, 07, 31)
    nb_lesions = 1
    # Run the simulation with each model
    severity1, nb_DUs1, nb_lesions1 = simulate_severity(start_date, end_date, nb_lesions,
                                                        model=model,
                                                        update_events="step_by_step")                                                        
    severity2, nb_DUs2, nb_lesions2 = simulate_severity(start_date, end_date, nb_lesions,
                                                        model=model,
                                                        update_events = "dispersal")
    
    # Display results
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(4,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), severity1, color='b')
    plot(pandas.date_range(start_date, end_date, freq='H'), severity2, color='r')
    ylim([0, 105])
    ylabel('Total disease severity')
    
    ax2 = fig.add_subplot(4,1,2)
    plt.bar(range(nb_steps), nb_DUs1, color='b')
    ylabel('Number of DUs')
    xlim([0, nb_steps])
    ax3 = fig.add_subplot(4,1,3)
    plt.bar(range(nb_steps), nb_DUs2, color='r')
    ylabel('Number of DUs')
    xlim([0, nb_steps])

    ax4 = fig.add_subplot(4,1,4)
    plot(range(nb_steps), nb_lesions1, color='b')
    plot(range(nb_steps), nb_lesions2, color='r')
    ylabel('Number of lesions')
    xlim([0, nb_steps])
    xlabel('Simulation time step')
    
    fig.subplots_adjust(hspace=1)
    
    plt.show()
    
def simulate_severity(start_date=datetime(2000, 10, 1),
                      end_date=datetime(2000, 11, 30),
                      nb_lesions = 100,
                      model="SeptoriaExchangingRings",
                      update_events = "step_by_step"):
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.txt'
    weather_data = read_weather(filename, update_events)
    
    # Generate a MTG with required properties :
    g = adel_one_leaf()
    # g = adel_mtg2()
    set_initial_properties_g(g)
    
    # Generate a stock of septoria lesions and distribute them on g
    nb_lesions_in_stock = nb_lesions
    # distribute_dispersal_units(g, nb_dus=nb_lesions_in_stock, model=model)
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, model=model)

    # Prepare the simulation loop
    dt = 1
    severity = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    nb_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    nb_lesions = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try: ind+=1.
        except: ind = 0.

        g = execute_one_step(g, weather_data, start_date, date, dt)
        
        # Measure severity
        severity[ind]=compute_total_severity(g)
        
        nb_DUs[ind] = count_DU(g)
        nb_lesions[ind] = count_lesions(g)
    
    # Display results:
    # nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    # fig = plt.figure()
    
    # ax1 = fig.add_subplot(3,1,1)
    # plot(pandas.date_range(start_date, end_date, freq='H'), severity)
    # ylim([0, 105])
    # ylabel('Total disease severity')
    
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
    
    # plt.show()
    
    return severity, nb_DUs, nb_lesions

def imitate_septo3D(start_date=datetime(2000, 10, 1),
                      end_date=datetime(2000, 11, 30),
                      nb_lesions = 100,
                      model="SeptoriaExchangingRings",
                      update_events = "step_by_step"):
    """ Run a simulation with the model of continuous lesions of septoria.
    
    Lesions are updated only on dispersal events. Weather data is managed
    by an adapter that pass it to the fungal model in the good format.

    Parameters
    ----------
        None
    """
    # Read weather data : 
    filename = 'meteo01.txt'
    weather_data = read_weather(filename, update_events)
    
    # Generate a MTG with required properties :
    g = adel_one_leaf()
    # g = adel_mtg2()
    set_initial_properties_g(g, surface_leaf_element=10.)
    
    # Generate a stock of septoria lesions and distribute them on g
    nb_lesions_in_stock = nb_lesions
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, model=model)
    
    # Call the models that will be used during the simulation :
    position_checker=BiotrophDUProbaModel()
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()
    leaf_inspector1 = LeafInspector(g, leaf_number=1)
    # leaf_inspector2 = LeafInspector(g, leaf_number=2)
    # leaf_inspector3 = LeafInspector(g, leaf_number=3)
    # leaf_inspector4 = LeafInspector(g, leaf_number=4)

    # Prepare the simulation loop
    dt = 1
    degree_days = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    healthy_surface = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    total_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    unwashed_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    on_green_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    incubating_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
       
    for date in pandas.date_range(start_date, end_date, freq='H'):
        print(date)
        try: 
            ind+=1.
            # Estimate delta degree days
            degree_days[ind] = degree_days[ind-1] + compute_hourly_delta_ddays(temp=weather_data.Tair[date])
        except: 
            ind = 0.
            # Estimate delta degree days
            degree_days[ind] = compute_hourly_delta_ddays(temp=weather_data.Tair[date])
        
        # Update climate on leaf elements
        update_on_leaves(weather_data, date, g)     
        
        # Run a disease cycle
        infect(g, dt, position_checker)
        if weather_data.update_events[date]:
            update_temp_list_on_g(g, weather_data, start_date, date)
            update(g, dt, growth_control_model=controler)

        # Healthy surface
        healthy_surface[ind] = leaf_inspector1.compute_healthy_surface(g)
        
        # Surfaces in state
        leaf_inspector1.compute_states(g)
        # leaf_inspector2.compute_states(g)
        # leaf_inspector3.compute_states(g)
        # leaf_inspector4.compute_states(g)

        global_rain_intensity = weather_data.Pluie[date]
        if global_rain_intensity != 0. :
            initial_DUs = leaf_inspector1.count_DU(g)
            scene = plot3d(g)
            seed(1)
            disperse(g, scene, dispersor, "Septoria")
            seed(1)
            total_DUs[ind] = leaf_inspector1.count_DU(g) - initial_DUs
            wash(g, washor, global_rain_intensity)
            unwashed_DUs[ind] = leaf_inspector1.count_DU(g) - initial_DUs
            on_green_DUs[ind] = leaf_inspector1.count_DU_on_green(g, unwashed_DUs[ind])
            # print('pluie')

    # raise Exception('')
    
    # Draw the outputs:
    # 1. Build the frame
    fig = plt.figure()
    nb_graphs = 10
    length = 5
    width = nb_graphs/length
    all_letters = string.ascii_uppercase
    for i in range(nb_graphs):
        graph_handle = 'ax'+int2str(i+1)
        # if i < 4:
            # globals()[graph_handle] = fig.add_subplot(length,width,i+1,
                                                      # xticklabels=[], yscale='log')
        # elif i < 8:
        if i < 8:
            globals()[graph_handle] = fig.add_subplot(length,width,i+1, xticklabels=[])
        else:
            globals()[graph_handle] = fig.add_subplot(length,width,i+1)
        graph_letter = all_letters[i]
        globals()[graph_letter] = plt.text(0.05, 0.9, graph_letter, 
                                        fontsize=18, ha='center',
                                        va='center', transform=globals()[graph_handle].transAxes)

    # Titles
    droplets = plt.text(-0.275, 2.1, 'Number of Infectious Droplets', 
                        fontsize=18, rotation='vertical', transform=ax3.transAxes)
    leaf_area = plt.text(-0.275, 1.4, 'Leaf Area', 
                        fontsize=18, rotation='vertical', transform=ax7.transAxes)
    coverage = plt.text(-0.275, 0.75, '% Recovered', 
                        fontsize=18, rotation='vertical', transform=ax9.transAxes)
        
    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    # ax1.bar(pandas.date_range(start_date, end_date, freq='H'), total_DUs,  color='b')
    # ax1.set_yscale('log')
    ax1.bar(degree_days, total_DUs,  color='b')
    # ax1.plot(pandas.date_range(start_date, end_date, freq='H'), washing_rate)   
    # ax1.plot(degree_days, washing_rate)   
    ax1.set_ylabel('Total intercepted')
    ax1.set_xlim([0, max(degree_days)])
    ax1.set_ylim([0, max(total_DUs)+20])
    
    ax2.bar(degree_days, unwashed_DUs,  color='b')
    ax2.set_xlim([0, max(degree_days)])
    ax2.set_ylim([0, max(total_DUs)+20])
    ax2.set_ylabel('Total unwashed')
    
    ax3.bar(degree_days, on_green_DUs,  color='b')
    ax3.set_xlim([0, max(degree_days)])
    ax3.set_ylim([0, max(total_DUs)+20])
    ax3.set_ylabel('Total on Green')
    
    ax33 = ax3.twinx()
    ax33.plot(degree_days, healthy_surface,  color='b', linestyle='--')
    ax33.set_xticklabels([])
    ax33.set_xlim([0, max(degree_days)])
    
    ax5.plot(degree_days, leaf_inspector1.ratio_inc, color='b')
    ax5.set_ylim([0, 105])
    ax5.set_xlim([0, max(degree_days)])
    
    ax6.plot(degree_days, leaf_inspector1.ratio_chlo, color='b')
    ax6.set_ylim([0, 105])
    ax6.set_xlim([0, max(degree_days)])
    
    ax7.plot(degree_days, leaf_inspector1.ratio_nec, color='b')
    ax7.set_ylim([0, 105])
    ax7.set_xlim([0, max(degree_days)])
    
    ax8.plot(degree_days, leaf_inspector1.ratio_spo, color='b')
    ax8.set_ylim([0, 105])
    ax8.set_xlim([0, max(degree_days)])
    
    # raise Exception('')
    
    # ax2 = fig.add_subplot(3,1,2)
    # plot(pandas.date_range(start_date, end_date, freq='H'), weather_data.Pluie[start_date:end_date])
    # ylabel('Rain intensity')
    
    # ax3 = fig.add_subplot(3,1,3)
    # plot(pandas.date_range(start_date, end_date, freq='H'), weather_data.rain_duration[start_date:end_date])
    # ylabel('Rain duration')
    
    # ax2 = fig.add_subplot(1,2,2)
    # plt.bar(pandas.date_range(start_date, end_date, freq='H'), unwashed_DUs, color='b')
    # ylabel('Total unwashed')
    # ax2.set_yscale('log')
    # xlim([0, nb_steps])

    # ax3 = fig.add_subplot(3,1,3)
    # plot(range(nb_steps), nb_lesions, color='r')
    # ylabel('Number of lesions')
    # xlim([0, nb_steps])
    # xlabel('Simulation time step')
    
    # fig.subplots_adjust(hspace=1)
    
    plt.show(False)

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
    
    # Generate a stock of septoria lesions and distribute them on g
    nb_lesions_in_stock = 100
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, model=model)
    
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
        if weather_data.update_events[date]:
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
    
def simulate_states_one_leaf(model="SeptoriaExchangingRings"):
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
    
    # Generate a stock of septoria lesions and distribute them on g
    nb_lesions_in_stock = 10
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, model=model)
    
    # Call the models that will be used during the simulation :
    position_checker = BiotrophDUProbaModel()
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()
    
    # Call the variable saver
    outputs = LeafInspector(g, leaf_number=1)
    
    # Prepare the simulation loop
    dt = 1
    start_date = datetime(2000, 10, 1)
    end_date = datetime(2001, 2, 1)
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
        infect(g, dt) # , position_checker)
        if weather_data.update_events[date]:
            # print(date)
            update_temp_list_on_g(g, weather_data, start_date, date)
            update(g, dt, growth_control_model=controler)
                
        global_rain_intensity = weather_data.Pluie[date]
        if global_rain_intensity != 0. :
            scene = plot3d(g)
            disperse(g, scene, dispersor, "Septoria") 
            wash(g, washor, global_rain_intensity)
    
        # Save variables
        outputs.compute_severity(g)
        # outputs.compute_states(g)

    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    # y_max = max([max(incubating), max(chlorotic), max(necrotic), max(sporulating)])+1
    
    # fig = plt.figure()
    
    # ax1 = fig.add_subplot(4,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), outputs.severity)
    # plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_inc)
    # ylim([0, y_max])
    # ylabel('Incubating surface')
    
    # ax2 = fig.add_subplot(4,1,2)
    # plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_chlo, color='g')
    # ylabel('Chlorotic surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])

    # ax3 = fig.add_subplot(4,1,3)
    # plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_nec, color='r')
    # ylabel('Necrotic surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])

    # ax4 = fig.add_subplot(4,1,4)
    # plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_spo, color='k')
    # plot(pandas.date_range(start_date, end_date, freq='H'), created_surface, color='b', linestyle='--')
    # plot(pandas.date_range(start_date, end_date, freq='H'), sum_surface, color='k', linestyle='--')
    # ylabel('Sporulating surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])
    # xlabel('Simulation time step')
    
    # fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g