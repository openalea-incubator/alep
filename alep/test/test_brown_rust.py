""" Tests for brown rust model.
"""
from alinea.alep.brown_rust import *
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.adel.data_samples import adel_two_metamers_stand
from alinea.alep.architecture import get_leaves, set_properties
from openalea.vpltk import plugin
import numpy as np
import pandas as pd

def get_small_g():
    g, domain_area, domain, convunit = adel_two_metamers_stand(leaf_sectors = 1,
                                                               density = 350., 
                                                               interleaf = 10.,
                                                               leaf_length = 20,
                                                               leaf_width = 1, Einc = 0)
    return g
   
def get_g_and_one_leaf():
    g = get_small_g()
    leaves = get_leaves(g)
    return g, g.node(leaves[0])
   
def get_one_leaf():
    g = get_small_g()
    leaves = get_leaves(g)
    return g.node(leaves[0])

def test_du_instanciation():
    brown_rust = BrownRustFungus()
    du = brown_rust.dispersal_unit()
    assert issubclass(du.__class__, DispersalUnit)
    assert isinstance(du, BrownRustDU)
    assert du.is_active
    
def test_lesion_instanciation():
    brown_rust = BrownRustFungus()
    lesion = brown_rust.lesion()
    assert issubclass(lesion.__class__, Lesion)
    assert isinstance(lesion, BrownRustLesion)
    assert lesion.is_active

def test_brown_rust_plugin():
    diseases=plugin.discover('alep.disease')
    brown_rust = diseases['brown_rust'].load()
    assert issubclass(BrownRustFungus, Fungus)
    assert brown_rust == BrownRustFungus

def test_du_mutability():
    brown_rust = BrownRustFungus()
    du1 = brown_rust.dispersal_unit(mutable = True, temp_opt_inf = 10.)
    du2 = brown_rust.dispersal_unit(mutable = True, temp_opt_inf = 20.)
    assert du1.fungus.temp_opt_inf != du2.fungus.temp_opt_inf
    
def test_lesion_mutability():
    brown_rust = BrownRustFungus()
    les1 = brown_rust.lesion(mutable = True, temp_opt_inf = 10.)
    les2 = brown_rust.lesion(mutable = True, temp_opt_inf = 20.)
    assert les1.fungus.temp_opt_inf != les2.fungus.temp_opt_inf
    
def test_group_dus_in_cohort(nb_dus_in_cohort = 20):
    brown_rust = BrownRustFungus()
    brown_rust.parameters(group_dus=True)
    du = brown_rust.dispersal_unit()
    du.set_position([[np.random.random(), np.random.random()] 
                     for i in range(nb_dus_in_cohort)])
    assert du.nb_dispersal_units == nb_dus_in_cohort
    
def test_group_lesions_in_cohort(nb_lesions_in_cohort = 20):
    brown_rust = BrownRustFungus()
    brown_rust.parameters(group_dus=True)
    lesion = brown_rust.lesion()
    lesion.set_position([[np.random.random(), np.random.random()]
                         for i in range(nb_lesions_in_cohort)])
    assert lesion.nb_lesions == nb_lesions_in_cohort

def test_lesion_creation_by_du():
    # Individual dispersal unit
    brown_rust = BrownRustFungus()
    du1 = brown_rust.dispersal_unit(mutable = True, group_dus = False)
    du1.set_position([np.random.random(), np.random.random()])
    les1 = du1.create_lesion()
    assert issubclass(les1.__class__, Lesion)
    assert isinstance(les1, BrownRustLesion)
    assert les1.is_active
    assert du1.position == les1.position
    assert du1.is_active == False
    
    # Cohort
    du2 = brown_rust.dispersal_unit(mutable = True, group_dus = True)
    du2.set_position([[np.random.random(), np.random.random()] for i in range(20)])
    les2 = du2.create_lesion()
    assert issubclass(les2.__class__, Lesion)
    assert isinstance(les2, BrownRustLesion)
    assert les2.is_active
    assert les2.position == les2.position
    assert du2.is_active == False

def test_infection_impossible_on_senescence():
    g = get_small_g()
    set_properties(g, label = 'LeafElement', senesced_length = 1.)
    leaves = get_leaves(g)
    # Individual dispersal unit
    leaf1 = leaves[0]
    leaf = g.node(leaf1)
    brown_rust = BrownRustFungus()
    du1 = brown_rust.dispersal_unit(mutable = True, group_dus = False)
    du1.set_position([np.random.random(), np.random.random()])
    leaf.dispersal_units = [du1]
    assert du1.is_active
    du1.disable_if_senescent(leaf)
    assert du1.is_active == False

    # Cohort
    leaf2 = leaves[1]
    leaf = g.node(leaf2)
    du2 = brown_rust.dispersal_unit(mutable = True, group_dus = True)
    du2.set_position([[np.random.random(), np.random.random()] for i in range(20)])
    g.node(leaf2).dispersal_units = [du2]
    assert du2.is_active
    du2.disable_if_senescent(leaf)
    assert du2.is_active == False
    
def test_infection_temperature():
    
    def infection_test_one_temp(temp = 15.):
        leaf = get_one_leaf()
        leaf.lesions = []
        len_seq = 30
        df_temp = pd.DataFrame([temp for i in range(len_seq)], columns = ['temp'])
        df_wet = pd.DataFrame([True for i in range(len_seq)], columns = ['wet'])
        leaf.temperature_sequence = df_temp['temp']
        leaf.wetness_sequence = df_wet['wet']
        dus = [brown_rust.dispersal_unit(group_dus = False) for i in range(nb_dus)]
        for du in dus:
            du.infect(dt = len_seq, leaf = leaf)
        return len(leaf.lesions)
    
    brown_rust = BrownRustFungus()
    nb_dus = 100
    temp_max = brown_rust.parameters()['temp_max_inf']
    temp_min = brown_rust.parameters()['temp_min_inf']
    temp_opt = brown_rust.parameters()['temp_opt_inf']
    assert infection_test_one_temp(temp = temp_opt) >= 0.85 * nb_dus
    assert infection_test_one_temp(temp = temp_max+1.) == 0.
    assert infection_test_one_temp(temp = temp_min-1.) == 0.
    assert infection_test_one_temp(temp = np.random.randint(temp_min, temp_opt-1)) < nb_dus
    assert infection_test_one_temp(temp = np.random.randint(temp_opt, temp_max-1)) > 0.

def test_infection_temperature_grouped():
    
    def infection_test_one_temp(temp = 15.):
        leaf = get_one_leaf()
        leaf.lesions = []
        len_seq = 30
        df_temp = pd.DataFrame([temp for i in range(len_seq)], columns = ['temp'])
        df_wet = pd.DataFrame([True for i in range(len_seq)], columns = ['wet'])
        leaf.temperature_sequence = df_temp['temp']
        leaf.wetness_sequence = df_wet['wet']
        du = brown_rust.dispersal_unit(group_dus = True)
        du.set_position([[np.random.random(), np.random.random()] for i in range(nb_dus)])
        du.infect(dt = len_seq, leaf = leaf)
        nb_infect = sum([les.nb_lesions for les in leaf.lesions])
        assert nb_infect + du.nb_dispersal_units <= nb_dus
        return nb_infect
    
    brown_rust = BrownRustFungus()
    nb_dus = 100
    temp_max = brown_rust.parameters()['temp_max_inf']
    temp_min = brown_rust.parameters()['temp_min_inf']
    temp_opt = brown_rust.parameters()['temp_opt_inf']
    assert infection_test_one_temp(temp = temp_opt) >= 0.85 * nb_dus
    assert infection_test_one_temp(temp = temp_max+1.) == 0.
    assert infection_test_one_temp(temp = temp_min-1.) == 0.
    assert infection_test_one_temp(temp = np.random.randint(temp_min, temp_opt-1)) < nb_dus
    assert infection_test_one_temp(temp = np.random.randint(temp_opt, temp_max-1)) > 0.
    
def test_infection_wetness():
        
    def infection_test_one_wetness_duration(wet_duration = 20):
        leaf = get_one_leaf()
        leaf.lesions = []
        df_temp = pd.DataFrame([brown_rust.parameters()['temp_opt_inf'] for i in range(len_seq)],
                                columns = ['temp'])
        wd = np.zeros(len_seq)
        wd[:wet_duration] = 1
        np.random.shuffle(wd)
        df_wet = pd.DataFrame(wd, columns = ['wet'])
        leaf.temperature_sequence = df_temp['temp']
        leaf.wetness_sequence = df_wet['wet']
        du = brown_rust.dispersal_unit(group_dus = True)
        du.set_position([[np.random.random(), np.random.random()] for i in range(nb_dus)])
        du.infect(dt = len_seq, leaf = leaf)
        nb_infect = sum([les.nb_lesions for les in leaf.lesions])
        assert nb_infect + du.nb_dispersal_units <= nb_dus
        return nb_infect
    
    brown_rust = BrownRustFungus()
    nb_dus = 100
    len_seq = 30
    wetness_min = brown_rust.parameters()['wetness_min']
    assert infection_test_one_wetness_duration(wet_duration = len_seq) == nb_dus
    assert infection_test_one_wetness_duration(wet_duration = 0) == 0.
    assert infection_test_one_wetness_duration(wet_duration = np.random.randint(wetness_min, len_seq)) <= nb_dus
    
def test_infection_wetness_grouped():
        
    def infection_test_one_wetness_duration(wet_duration = 20):
        leaf = get_one_leaf()
        leaf.lesions = []
        df_temp = pd.DataFrame([brown_rust.parameters()['temp_opt_inf'] for i in range(len_seq)],
                                columns = ['temp'])
        wd = np.zeros(len_seq)
        wd[:wet_duration] = 1
        np.random.shuffle(wd)
        df_wet = pd.DataFrame(wd, columns = ['wet'])
        leaf.temperature_sequence = df_temp['temp']
        leaf.wetness_sequence = df_wet['wet']
        dus = [brown_rust.dispersal_unit(group_dus = False) for i in range(nb_dus)]
        for du in dus:
            du.infect(dt = len_seq, leaf = leaf)
        return len(leaf.lesions)
    
    brown_rust = BrownRustFungus()
    nb_dus = 100
    len_seq = 30
    wetness_min = brown_rust.parameters()['wetness_min']
    assert infection_test_one_wetness_duration(wet_duration = len_seq) == nb_dus
    assert infection_test_one_wetness_duration(wet_duration = 0) == 0.
    assert infection_test_one_wetness_duration(wet_duration = np.random.randint(wetness_min, len_seq)) <= nb_dus
    
def test_loss_delay():
    brown_rust = BrownRustFungus()
    # Individual dispersal unit
    leaf = get_one_leaf()
    loss_delay = int(brown_rust.parameters()['loss_delay'])
    df_temp = pd.DataFrame([brown_rust.parameters()['temp_opt_inf'] for i in range(loss_delay+1)],
                                columns = ['temp'])
    df_wet = pd.DataFrame([False for i in range(loss_delay+1)], columns = ['wet'])
    leaf.temperature_sequence = df_temp['temp']
    leaf.wetness_sequence = df_wet['wet']
    du = brown_rust.dispersal_unit(mutable = True, group_dus = False)
    assert du.is_active
    du.infect(dt = loss_delay+1, leaf = leaf)
    assert du.is_active == False
    
    # Cohort
    leaf = get_one_leaf()
    loss_delay = int(brown_rust.parameters()['loss_delay'])
    df_temp = pd.DataFrame([brown_rust.parameters()['temp_opt_inf'] for i in range(loss_delay+1)],
                                columns = ['temp'])
    df_wet = pd.DataFrame([False for i in range(loss_delay+1)], columns = ['wet'])
    leaf.temperature_sequence = df_temp['temp']
    leaf.wetness_sequence = df_wet['wet']
    du = brown_rust.dispersal_unit(mutable = True, group_dus = True)
    du.set_position([[np.random.random(), np.random.random()] for i in range(20)])
    assert du.is_active
    du.infect(dt = loss_delay+1, leaf = leaf)
    assert du.nb_dispersal_units == 0.
    assert du.is_active == False
    
def test_thermal_time():
    brown_rust = BrownRustFungus()
    # Individual dispersal unit
    leaf = get_one_leaf()
    latency_min = brown_rust.parameters()['latency_min']
    df_temp = pd.DataFrame([brown_rust.parameters()['temp_opt_lat'] for i in range(24)],
                                columns = ['temp'])
    leaf.temperature_sequence = df_temp['temp']
    lesion = brown_rust.lesion(mutable = True, group_dus = False)
    lesion.update(leaf = leaf)
    age1 = lesion.age_tt
    assert age1 == brown_rust.parameters()['temp_opt_lat']
    lesion.update(leaf = leaf)
    age2 = lesion.age_tt
    assert age2 > age1
    leaf.temperature_sequence = [0.]
    lesion.update(leaf = leaf)
    age3 = lesion.age_tt
    assert age3 == age2
    
    # Cohort
    leaf = get_one_leaf()
    leaf.temperature_sequence = df_temp['temp']
    lesion = brown_rust.lesion(mutable = True, group_dus = True)
    lesion.set_position([[np.random.random(), np.random.random()] for i in range(20)])
    lesion.update(leaf = leaf)
    age1 = lesion.age_tt
    assert age1 == brown_rust.parameters()['temp_opt_lat']
    lesion.update(leaf = leaf)
    age2 = lesion.age_tt
    assert age2 > age1
    leaf.temperature_sequence = [0.]
    lesion.update(leaf = leaf)
    age3 = lesion.age_tt
    assert age3 == age2

def test_latent_period():
    brown_rust = BrownRustFungus()
    # Individual dispersal unit
    leaf = get_one_leaf()
    latency_min = brown_rust.parameters()['latency_min']
    df_temp = pd.DataFrame([brown_rust.parameters()['temp_opt_lat'] for i in range(50)],
                                columns = ['temp'])
    leaf.temperature_sequence = df_temp['temp']
    lesion = brown_rust.lesion(mutable = False, group_dus = False)
    for i in range(6):
        lesion.update(leaf = leaf)
    assert lesion.age_physio_lat >= 1.
    assert lesion.is_sporulating
    
    # Cohort
    leaf = get_one_leaf()
    leaf.temperature_sequence = df_temp['temp']
    lesion = brown_rust.lesion(mutable = True, group_dus = True)
    lesion.set_position([[np.random.random(), np.random.random()] for i in range(20)])
    for i in range(6):
        lesion.update(leaf = leaf)
    assert lesion.age_physio_lat >= 1.
    assert lesion.is_sporulating
    
def test_surface():
    brown_rust = BrownRustFungus()
    # Individual dispersal unit
    leaf = get_one_leaf()
    latency_min = brown_rust.parameters()['latency_min']
    df_temp = pd.DataFrame([brown_rust.parameters()['temp_opt_lat']], columns = ['temp'])
    leaf.temperature_sequence = df_temp['temp']
    lesion = brown_rust.lesion(mutable = False, group_dus = False, latency_min = 175.)
    surfs = []
    surfs_lat = []
    surfs_spo = []
    for i in range(1000):
        lesion.update(leaf = leaf)
        lesion.control_growth(growth_offer = lesion.growth_demand)
        surfs.append(lesion.surface)
        surfs_lat.append(lesion.surface_lat)
        surfs_spo.append(lesion.surface_spo)
    
    from alinea.echap.weather_data import read_weather_year
    weather = read_weather_year(2012)
    df = weather.data
    df = df.reset_index()
    df = df[df['degree_days']>0]
    brown_rust = BrownRustFungus()
    # Individual dispersal unit
    leaf = get_one_leaf()
    latency_min = brown_rust.parameters()['latency_min']    
    lesion = brown_rust.lesion(mutable = False, group_dus = False, latency_min = 175.)
    surfs = []
    surfs_lat = []
    surfs_spo = []
    tt = []
    cum_tt = 0.
    count = 0.
    for i, row in df.iterrows():
        if i > count < 3000:
            count += 1
            leaf.temperature_sequence = [row['temperature_air']]
            lesion.update(leaf = leaf)
            lesion.control_growth(growth_offer = lesion.growth_demand)
            surfs.append(lesion.surface)
            surfs_lat.append(lesion.surface_lat)
            surfs_spo.append(lesion.surface_spo)
            cum_tt += lesion.delta_thermal_time_growth(leaf_temperature = [row['temperature_air']])
            tt.append(cum_tt)
    
    
    brown_rust = BrownRustFungus()
    # Individual dispersal unit
    leaf = get_one_leaf()
    latency_min = brown_rust.parameters()['latency_min']
    df_temp = pd.DataFrame([20.], columns = ['temp'])
    leaf.temperature_sequence = df_temp['temp']
    lesion = brown_rust.lesion(mutable = False, group_dus = False, latency_min = 175.)
    surfs = []
    surfs_lat = []
    surfs_spo = []
    tt = []
    cum_tt = 0.
    for i in range(24):
        lesion.update(leaf = leaf)
        lesion.control_growth(growth_offer = lesion.growth_demand)
        cum_tt += lesion.delta_thermal_time_growth(leaf_temperature = df_temp['temp'])
        tt.append(cum_tt)
    
    # Penser a faire update puis control growth avec growth offer non limitante
    # Scenario with updates long time step
    # No negative surface if temperature negative
    # No increase of surface if negative temperature
    # Surface always inferior to Smax
    # Surface tot = surface spo + surface lat
    # Surface spo always < 0.37 * Smax
    # trouver un test : nb lesions * surface ?
    # time step < delai latence : surface spo == 0
    # time step > delai latence : fonctionne ?
    # time step tres long : surface lat = 0
    # time step intermediaire : surface lat > 0 + surface spo > 0
    # balance on one time step
    pass
    
def test_density():
    brown_rust = BrownRustFungus()
    growth_controler = NoPriorityGrowthControl()
    
    def test_single_density(d = 40):
        # Individual dispersal unit
        g, leaf = get_g_and_one_leaf()
        leaf.area = 1.
        leaf.green_area = 1.
        latency_min = brown_rust.parameters()['latency_min']
        df_temp = pd.DataFrame([brown_rust.parameters()['temp_opt_lat']], columns = ['temp'])
        leaf.temperature_sequence = df_temp['temp']
        leaf.lesions = [brown_rust.lesion(mutable = False, group_dus = False,
                        Smax = 0.22, sporulating_capacity = 0.2) for i in range(d)]
        for i in range(1000):
            for les in leaf.lesions:
                les.update(leaf = leaf)
            growth_controler.control(g)
        surfs = [les.surface for les in leaf.lesions]
        surfs_spo = [les.surface_spo for les in leaf.lesions]
        return np.mean(surfs), np.mean(surfs_spo)
    
    surf_mean = []
    surf_spo_mean = []
    for d in np.arange(1, 41, 2):
        print d
        surf, surf_spo = test_single_density(d)
        surf_mean.append(surf)
        surf_spo_mean.append(surf_spo)
        
    import matplotlib.pyplot as plt
    df = pd.read_csv('rust_density.csv', sep = ';')
    fig, ax = plt.subplots()
    ax.plot(np.arange(1, 41, 2), surf_spo_mean)
    ax.plot(df.density, df.surface)
    
def test_logistic():
    # test strictly growing
    # test stable if progress = 0
    pass
    
def test_senescence_response():
    # indiv tout meurt
    # cohort verif nb lesions diminue
    # verif plus de lesion dans partie senescente
    # verif activite
    pass

def test_competition():
    pass
    
def test_emptiness():
    pass
    