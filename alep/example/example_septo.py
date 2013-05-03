""" Example of use of the model of continuous lesions of septoria with real weather data.

"""
# Imports #########################################################################
import random as rd
import numpy as np
import pandas
from pylab import *
import matplotlib.pyplot as plt

from alinea.alep.wheat_examples import adel_mtg, adel_mtg2, adel_mtg3, adel_one_leaf
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

from alinea.caribu.CaribuScene import CaribuScene
import alinea.caribu.sky_tools.turtle as turtle

from alinea.alep.protocol import *
from alinea.alep.septoria import *
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
    # Use CaribuMicroclimModel for rain
    # scene = plot3d(g)
    # rain = data.Pluie[date]
    # microclim_model = CaribuMicroclimModel(energy=data.PAR[date])
    # microclimate = {}
    # microclim = microclim_model.microclim(microclimate, rain, scene)
    
    # Equals weather data on the entire plant for the rest
    vids = [n for n in g if g.label(n).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.temp = data.Tair[date]
        # n.rain_intensity = microclim[v]['rain'] * rain
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
    
class CaribuMicroclimModel(object):
    """ Adaptor for Caribu model compliying echap microclimate_model protocol 
    """
    def __init__(self, sectors='46', energy=1):
        self.sectors = sectors
        self.energy = energy
    def microclim(self, microclimate, rain, scene):
        local_meteo = microclimate_leaf(self.sectors, self.energy, microclimate, rain, scene)
        return local_meteo

def microclimate_leaf(sectors, energy, microclimate, rain, scene):
    energy, emission, direction, elevation, azimuth = turtle.turtle(sectors, energy) 
    sources = zip(energy, direction)
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)    
    c_scene.addSources(sources)
    output = c_scene.runCaribu(infinity=False)
    if rain >0:
        rain_leaf = c_scene.output_by_id(output, idmap)['Einc']
    EiInf = c_scene.output_by_id(output, idmap)['EiInf']
    EiSup = c_scene.output_by_id(output, idmap)['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                if rain == 0:
                    microclimate[Infid] = {'radiation': e + a, 'rain': 0} 
                else:
                    microclimate[Infid] = {'radiation': e + a, 'rain': rain_leaf[Infid]} 
    return microclimate
        
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

def compute_states_one_leaf(g, blade_id=None):
    """ Compute surface of lesions in chosen state on a blade of the MTG.
    
    """
    # Find leaf elements on the blade
    surfaces = g.property('surface')
    labels = g.property('label')
    leaf_elements = [id for id in g.components(blade_id) if labels[id].startswith('LeafElement')]
    
    # Compute green surface on the blade
    green_surface = sum(surfaces[id] for id in leaf_elements)
    
    # Compute disease surface
    lesions = g.property('lesions')
    lesions = [l for l in lesions[id] for id in leaf_elements]
    if lesions:
        for l in lesions:
            l.compute_all_surfaces()
        
        surface_inc = sum(l.surface_inc for l in lesions)
        surface_chlo = sum(l.surface_chlo for l in lesions)
        surface_nec = sum(l.surface_inc for l in lesions)
        surface_spo = sum(l.surface_inc for l in lesions)
    else:
        surface_inc = 0.
        surface_chlo = 0.
        surface_nec = 0.
        surface_spo = 0.

    ratio_inc = surface_inc / green_surface if green_surface>0. else 0.
    ratio_chlo = surface_chlo / green_surface if green_surface>0. else 0.
    ratio_nec = surface_nec / green_surface if green_surface>0. else 0.
    ratio_spo = surface_spo / green_surface if green_surface>0. else 0.
    
    return ratio_inc, ratio_chlo, ratio_nec, ratio_spo

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
    
    def compute_states(self, g):
        """ Compute surface of lesions in chosen state on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        surfaces = g.property('surface')
        leaf_elements = self.leaf_elements
        
        # Compute green surface on the blade
        green_surface = sum(surfaces[id] for id in leaf_elements)
        
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

        self.ratio_inc.append(surface_inc / green_surface if green_surface>0. else 0.)
        self.ratio_chlo.append(surface_chlo / green_surface if green_surface>0. else 0.)
        self.ratio_nec.append(surface_nec / green_surface if green_surface>0. else 0.)
        self.ratio_spo.append(surface_spo / green_surface if green_surface>0. else 0.)
    
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
    end_date=datetime(2000, 11, 30)
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
    
    ax1 = fig.add_subplot(3,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), severity1, color='b')
    plot(pandas.date_range(start_date, end_date, freq='H'), severity2, color='r')
    ylim([0, 105])
    ylabel('Total disease severity')
    
    ax2 = fig.add_subplot(3,1,2)
    plt.bar(range(nb_steps), nb_DUs1, color='b')
    plt.bar(range(nb_steps), nb_DUs2, color='r')
    ylabel('Number of DUs')
    xlim([0, nb_steps])

    ax3 = fig.add_subplot(3,1,3)
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
    
    # Generate a stock of septoria dispersal units and distribute it on g
    nb_lesions_in_stock = nb_lesions
    # distribute_dispersal_units(g, nb_dus=nb_dus_in_stock, model=model)
    distribute_lesions(g, nb_lesions=nb_lesions_in_stock, model=model)
    
    # Call the models that will be used during the simulation :
    controler = NoPriorityGrowthControl()
    washor = RapillyWashing()
    dispersor = RandomDispersal()

    # Prepare the simulation loop
    dt = 1
    # start_date = datetime(2000, 10, 1)
    # end_date = datetime(2001, 7, 31)
    # end_date = datetime(2000, 11, 30)
    severity = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    nb_DUs = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    nb_lesions = np.zeros(len(pandas.date_range(start_date, end_date, freq='H')))
    for date in pandas.date_range(start_date, end_date, freq='H'):
        # print(date)
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
    
    # Generate a stock of septoria dispersal units and distribute it on g
    nb_dus_in_stock = 100
    distribute_dispersal_units(g, nb_dus=nb_dus_in_stock)
    
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
    
    # Generate a stock of septoria dispersal units and distribute it on g
    nb_dus_in_stock = 10
    distribute_dispersal_units(g, nb_dus=nb_dus_in_stock, model=model)
    
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
        infect(g, dt, position_checker)
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
        outputs.compute_states(g)

    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    # y_max = max([max(incubating), max(chlorotic), max(necrotic), max(sporulating)])+1
    
    # fig = plt.figure()
    
    # ax1 = fig.add_subplot(4,1,1)
    plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_inc)
    # ylim([0, y_max])
    # ylabel('Incubating surface')
    
    # ax2 = fig.add_subplot(4,1,2)
    plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_chlo, color='g')
    # ylabel('Chlorotic surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])

    # ax3 = fig.add_subplot(4,1,3)
    plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_nec, color='r')
    # ylabel('Necrotic surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])

    # ax4 = fig.add_subplot(4,1,4)
    plot(pandas.date_range(start_date, end_date, freq='H'), outputs.ratio_spo, color='k')
    # plot(pandas.date_range(start_date, end_date, freq='H'), created_surface, color='b', linestyle='--')
    # plot(pandas.date_range(start_date, end_date, freq='H'), sum_surface, color='k', linestyle='--')
    # ylabel('Sporulating surface')
    # ylim([0, y_max])
    # xlim([0, nb_steps])
    # xlabel('Simulation time step')
    
    # fig.subplots_adjust(hspace=1)
    
    plt.show()
    
    return g
    
def fake_simulate_states_one_leaf(model="SeptoriaExchangingRings"):
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
    
    # Generate a stock of septoria dispersal units and distribute it on g
    nb_dus_in_stock = 1
    distribute_dispersal_units(g, nb_dus=nb_dus_in_stock, model=model)
    
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
        update_climate_all(g, wetness=True, temp=22.)  
        
        # Run a disease cycle
        infect(g, dt)
        update(g, dt, growth_control_model=controler)
                   
        # Measure severity
        incubating[ind] = compute_state_one_leaf(g, status="INCUBATING")
        chlorotic[ind] = compute_state_one_leaf(g, status="CHLOROTIC")
        necrotic[ind] = compute_state_one_leaf(g, status="NECROTIC")
        sporulating[ind]= compute_state_one_leaf(g, status="SPORULATING")
    
    # Display results:
    nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
    
    plot(range(nb_steps), incubating)
    plot(range(nb_steps), chlorotic, color='g')
    plot(range(nb_steps), necrotic, color='r')
    plot(range(nb_steps), sporulating, color='k')
    # xlim([0, nb_steps])
    
    plt.show()
    
    return g