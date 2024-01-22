""" Test basic responses of the wheat septoria model """

# Imports #########################################################################
import random as rd

from openalea.core import plugin
from alinea.astk.TimeControl import *

from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *

from alinea.adel.data_samples import adel_one_leaf, adel_one_leaf_element
from alinea.alep.protocol import *
from alinea.alep.septoria import *
from alinea.alep.disease_operation import (distribute_dispersal_units,
                                           distribute_lesions)
from alinea.alep.disease_outputs import count_dispersal_units
from alinea.alep.architecture import set_properties, update_healthy_area
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.growth_control import NoPriorityGrowthControl
#from alinea.alep.senescence import WheatSeptoriaPositionedSenescence
#from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.alep.washing import RapillyWashing
from alinea.alep.diseases import get_disease

# Tests ###########################################################################
def test_initiate(model="septoria"):
    """ Check if 'initiate' from 'protocol.py' deposits dispersal units on the MTG.

    Generate a wheat MTG and distribute dispersal units randomly on leaf elements.
    Check that all the stock of DU is distributed in properties of leaf elements.

    Parameters
    ----------
    model: model
        One version of the model of septoria lesion among:
        'septoria_exchanging_rings', 'septoria_with_rings', 'septoria_continuous'
    """
    # Generate a wheat MTG
    g = adel_one_leaf_element()
    set_properties(g, label = 'LeafElement',
                   area=5., green_area=5., position_senescence=None)

    # Generate a stock of septoria dispersal units and distribute it on g
    nb_initial_dus = 100
    distribute_dispersal_units(g, nb_dus=nb_initial_dus,
                               disease_model=model,
                               initiation_model=RandomInoculation(),
                               label = 'LeafElement')

    # Check that all the DUs that were in the stock are now instantiated
    # on leaf elements :
    dispersal_units = g.property('dispersal_units')
    nb_dus_on_leaves = sum(len(dus) for dus in dispersal_units.values())
    assert nb_dus_on_leaves == nb_initial_dus

def test_infect(model="septoria"):
    """ Check if 'infect' from 'protocol.py' leads to infection by dispersal units
        on the MTG.

    Generate a wheat MTG and distribute dispersal units randomly on leaf elements.
    Run a short loop to compute infection. Check that all the stock of DU caused
    appearance of a lesion on leaf elements.

    Parameters
    ----------
    model: model
        One version of the model of septoria lesion among:
        'septoria_exchanging_rings', 'septoria_with_rings', 'septoria_continuous'
    """
    # Generate a wheat MTG
    g = adel_one_leaf_element()
    set_properties(g, label='LeafElement',
                   area=5., green_area=5., healthy_area=5., position_senescence=None,
                   temperature_sequence = [18] * 24,
                    relative_humidity_sequence = [85] * 24,
                    wetness_sequence = [True] * 24
    )

    # Generate a stock of septoria dispersal units and distribute it on g
    nb_dus = 100
    distribute_dispersal_units(g, nb_dus=nb_dus,
                               disease_model=model,
                               initiation_model=RandomInoculation(),
                               label = 'LeafElement')

    # Loop of simulation
    septo_timing = TimeControl(delay=1, steps = 10)
    timer = TimeControler(disease = septo_timing)
    for t in timer:
        # Offer good conditions for at least 10h
        set_properties(g, label = 'LeafElement', wetness=True, temp=20.)

        # Infect
        infect(g, t['disease'].dt, label = 'LeafElement')

    # Check that all the DUs that were in the stock are now lesions on leaf elements
    lesions = g.property('lesions')
    nb_lesions_on_leaves = sum(len(l) for l in lesions.values())
    assert nb_lesions_on_leaves == nb_dus

def test_update(model="septoria"):
    """ Check if 'update' from 'protocol.py' provokes the growth of a lesion
        instantiated on the MTG.

    Generate a wheat MTG and deposit 1 lesion on a leaf element. Run a loop
    to compute update and growth. Check that all the stages of a lesion of
    septoria have been reached eventually. Check that the lesion size does
    not exceed Smax.

    Parameters
    ----------
    model: model
        One version of the model of septoria lesion among:
        'septoria_exchanging_rings', 'septoria_with_rings', 'septoria_continuous'
    """
    # Generate a wheat MTG
    g = adel_one_leaf_element()
    set_properties(g, label = 'LeafElement',
                   area=5., green_area=5., healthy_area=5., position_senescence=None)

    # Generate a stock of septoria dispersal units and distribute it on g
    distribute_lesions(g, nb_lesions=1, disease_model=model,
                       initiation_model=RandomInoculation(),
                       label = 'LeafElement')

    # Call model of growth control
    controler = NoPriorityGrowthControl()

    # Loop of simulation
    septo_timing = TimeControl(delay=1, steps = 1000)
    timer = TimeControler(disease = septo_timing)
    for t in timer:
        # print(timer.numiter)
        # After infection, the lesion 'age_dday' will be added 1 DD by time step
        # Note : base temperature for septoria = -2 degrees celsius
        set_properties(g, label = 'LeafElement', wetness=True, temp=22., rain_intensity=0.)

        # Update
        update(g, t['disease'].dt, controler, label='LeafElement')

        # Check that the lesion is in the right status and has the right surface
        lesion = g.property('lesions')
        if lesion:
            assert sum(len(l) for l in lesion.values()) == 1
            l = list(lesion.values())[0][0]
            assert l.age_tt == timer.numiter
            f = l.fungus
            print(('lesion surface: %f' % round(l.surface, 6)))
            if l.age_dday < 220.:
                assert l.status == 0
                assert round(l.surface, 6) < round(f.Smin, 6)
            if 220. <= l.age_tt < 330.:
                assert l.status == 1
                assert round(f.Smin, 6) <= round(l.surface, 6)
                assert round(l.surface, 6) < round(f.Smin + f.growth_rate * f.degree_days_to_necrosis, 6)
            elif 330. <= l.age_tt < 350.:
                assert l.status == 2
                assert round(f.Smin + f.growth_rate * f.degree_days_to_necrosis, 6) <= round(l.surface, 6)
                assert round(l.surface, 6) < round(f.Smin + f.growth_rate * (f.degree_days_to_necrosis + f.degree_days_to_sporulation), 6)
            elif l.age_tt >= 350.:
                assert l.status == 3
                assert round(f.Smin + f.growth_rate * f.degree_days_to_sporulation, 6) <= round(l.surface, 6)
                assert round(l.surface, 6) <= round(l.fungus.Smax, 6)

def test_growth_control(model="septoria"):
    """ Check if 'control_growth' from 'protocol.py' limits the lesion growth
        up to available surface on the leaf.

    Generate a wheat MTG and deposit sumultaneously 1000 lesions on a leaf element.
    Run a loop to compute update and growth. Check that the healthy surface of the
    leaf decreases up to 0. Check that lesion growth is stopped after this point.

    Parameters
    ----------
    model: model
        One version of the model of septoria lesion among:
        'septoria_exchanging_rings', 'septoria_with_rings', 'septoria_continuous'
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    initial_leaf_area = 5.
    set_properties(g, label = 'LeafElement',
                   area=initial_leaf_area, green_area=initial_leaf_area,
                   healthy_area=initial_leaf_area, position_senescence=None)
    labels = g.property('label')
    total_initial_area = sum(g.node(vid).area for vid in labels.keys()
                             if labels[vid].startswith('LeafElement'))

    # Generate a stock of septoria dispersal units and distribute it on g
    distribute_lesions(g, nb_lesions=1000, disease_model=model,
                       initiation_model=RandomInoculation(),
                       label = 'LeafElement')

    # Call model of growth control
    controler = NoPriorityGrowthControl()

    # Initiation of healthy area the day before
    healthy_area_before = []

    # Loop of simulation
    septo_timing = TimeControl(delay=1, steps = 100)
    timer = TimeControler(disease = septo_timing)
    for t in timer:
        # print(timer.numiter)
        # After infection, the lesion 'age_dday' will be added 1 DD by time step
        # Note : base temperature for septoria = -2 degrees celsius
        set_properties(g, label = 'LeafElement', wetness=True, temp=22., rain_intensity=0.)

        # Update
        update(g, t['disease'].dt, controler, label='LeafElement')
        update_healthy_area(g, label = 'LeafElement')

        # Find the value of interest on the MTG (total healthy area of the leaf)
        lesions = g.property('lesions')
        healthy_areas = g.property('healthy_area')

        bids = (v for v,l in labels.items() if l.startswith('blade'))
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith('LeafElement')]
            leaf_healthy_area = sum(healthy_areas[lf] for lf in leaf)
            print(leaf_healthy_area)

        # Check that healthy area + lesion surface = initial healthy surface
        if lesions:
            leaf_lesion_area = sum(l.surface for les in lesions.values() for l in les)
            assert round(leaf_healthy_area,6)+round(leaf_lesion_area,6)==round(total_initial_area,6)

        # Check that healthy area decreases after emergence of the first lesion
        if lesions and not healthy_area_before:
            # Emergence of the first lesion
            healthy_area_before.append(leaf_healthy_area)
        elif lesions and healthy_area_before[0] > 0.:
            # Check that healthy area decreases
            assert leaf_healthy_area < healthy_area_before
            # Update of 'healthy_area_before' for next step
            healthy_area_before[0] = leaf_healthy_area
        elif lesions and healthy_area_before[0] == 0.:
            # Check that healthy area stays null
            assert leaf_healthy_area == 0.
            # Update of 'healthy_surface_before' for next step
            healthy_area_before[0] = leaf_healthy_area

def test_disperse(model="septoria"):
    """ Check if 'disperse' from 'protocol.py' disperse new
        dispersal units on the MTG.

    Generate a wheat MTG and distribute lesions randomly on leaf elements.
    Run a loop to compute infection and update. Create artificial rains in the loop
    on a regular basis. Check that the number of lesions on the MTG increases.

    Parameters
    ----------
    model: model
        One version of the model of septoria lesion among:
        'septoria_exchanging_rings', 'septoria_with_rings', 'septoria_continuous'
    """
    # Generate a wheat MTG
    g = adel_one_leaf()
    set_properties(g, label = 'LeafElement',
                   area=5., green_area=5., healthy_area=5.,
                   position_senescence=None)

    # Generate a stock of septoria dispersal units and distribute it on g
    distribute_lesions(g, nb_lesions=1, disease_model=model,
                       initiation_model=RandomInoculation(),
                       label = 'LeafElement')

    # Call models of growth control and dispersal
    controler = NoPriorityGrowthControl()
    emitter = SeptoriaRainEmission()
    transporter = Septo3DSplash()

    # Loop of simulation
    septo_timing = TimeControl(delay=1, steps = 1000)
    timer = TimeControler(disease = septo_timing)
    for t in timer:
        # Update climate and force rain occurences
        if timer.numiter>400 and timer.numiter%100 == 0:
            rain_intensity = 3.
        else:
            rain_intensity = 0.
        set_properties(g, label = 'LeafElement',
                       wetness=True, temp=22.,
                       relative_humidity=85.,
                       rain_intensity=rain_intensity)

        # Update
        update(g, t['disease'].dt, controler, label='LeafElement')

        # Force rain occurences
        if rain_intensity != 0:
            # Count objects on the MTG before dispersal event
            lesions = g.property('lesions')
            dispersal_units = g.property('dispersal_units')
            total_stock_spores_before = sum(l.stock_spores for les in list(lesions.values()) for l in les)
            total_DUs_before = sum(len(du) for du in dispersal_units.values())

            # Dispersal event
            disperse(g, dispersor, "septoria", label='LeafElement')

            # Check that stocks of spores on lesions decrease
            if total_stock_spores_before != 0.:
                total_stock_spores_after = sum(l.stock_spores for les in list(lesions.values()) for l in les)
                # print(total_stock_spores_after)
                assert total_stock_spores_after < total_stock_spores_before

                # Check that new DUs are deposited on the MTG
                total_DUs_after = sum(len(du) for du in dispersal_units.values())
                if total_DUs_after != 0:
                    assert total_DUs_after >= total_DUs_before

def test_senescence(status='CHLOROTIC', model="septoria"):
    """ Check if 'senescence' from 'protocol.py' compute the effects of
    senescence on the lesions of the MTG

    Generate a wheat MTG and distribute 2 DUs with know position on leaf elements.
    Make them grow into lesions until chosen status, then stop updating them.

    Set senescence on first lesion. Check that the activity of this lesion is reduced
    by senescence comparatively to the other lesion, but not the stock of spores.

    Parameters
    ----------
    status: str: 'CHLOROTIC' or 'SPORULATING'
        Status of the lesion when touched by senescence.
    model: model
        One version of the model of septoria lesion among:
        'septoria_exchanging_rings', 'septoria_with_rings', 'septoria_continuous'
    """
    # Generate a wheat MTG
    g = adel_one_leaf_element()
    set_properties(g, label = 'LeafElement',
                   area=5., green_area=5., healthy_area=5., position_senescence=None)

    # Generate a stock of septoria dispersal units and distribute it on g
    distribute_lesions(g, nb_lesions=2, disease_model=model,
                       initiation_model=RandomInoculation(),
                       label = 'LeafElement')

    # create a flat list of lesions and fix their position
    lesions = g.property('lesions')
    les = sum(list(lesions.values()), [])
    les[0].position = [0.3, 0]
    les[1].position = [0.7, 0]
    # Call a models growth control and senescence
    controler = NoPriorityGrowthControl()
    senescence_model = WheatSeptoriaPositionedSenescence(g)

    # Test if lesion is CHLOROTIC when senescence occur
    if status == "CHLOROTIC":
        if model != "septoria_exchanging_rings":
            dt = 300
            set_properties(g, label = 'LeafElement', wetness=True, temp=22., rain_intensity=0.)
            # Simulation to obtain a lesion in chosen status
            update(g, dt = dt, growth_control_model = controler, senescence_model=senescence_model)
        else:
            dt=1
            nb_steps=300
            for i in range(0,nb_steps,dt):
                # Simulation to obtain a lesion in chosen status
                set_properties(g, label = 'LeafElement', wetness=True, temp=22., rain_intensity=0.)
                update(g, dt = dt, growth_control_model = controler, senescence_model=senescence_model)

        # 1. Compare lesions before senescence
        # Compute variables of interest
        l = g.property('lesions')
        l = sum(list(l.values()), [])
        lesion1 = l[0]
        lesion1.compute_all_surfaces()
        lesion2 = l[1]
        lesion2.compute_all_surfaces()
        # Compare the two lesions
        assert lesion1.status == lesion2.status == 1
        assert lesion1.surface_alive == lesion2.surface_alive
        assert lesion1.surface_dead == lesion2.surface_dead == 0.

        # 2. Set senescence, update lesion so they know they
        #    are on a senescent tissue and compare lesions
        set_properties(g, label = 'LeafElement', position_senescence=0.6)
        update(g, dt = 0., growth_control_model = controler, senescence_model=senescence_model)

        # Compute variables of interest
        l = g.property('lesions')
        l = sum(list(l.values()), [])
        lesion1 = l[0]
        lesion1.compute_all_surfaces()
        lesion2 = l[1]
        lesion2.compute_all_surfaces()
        # Compare the two lesions
        assert lesion2.growth_is_active == False
        assert lesion2.is_active == False
        assert lesion1.surface_alive > lesion2.surface_alive
        assert lesion2.surface_alive == 0.
        assert lesion1.surface_dead == 0.
        assert lesion2.surface_dead > 0.
        assert lesion2.surface == lesion2.surface_dead

    # Test if lesion is SPORULATING when senescence occur
    elif status == "SPORULATING":
        if model != "septoria_exchanging_rings":
            dt = 400
            set_properties(g, label = 'LeafElement', wetness=True, temp=22., rain_intensity=0.)
            # Simulation to obtain a lesion in chosen status
            update(g, dt = dt, growth_control_model = controler, senescence_model=senescence_model)
        else:
            dt=1
            nb_steps=400
            for i in range(0,nb_steps,dt):
                # Simulation to obtain a lesion in chosen status
                set_properties(g, label = 'LeafElement', wetness=True, temp=22., rain_intensity=0.)
                update(g, dt = dt, growth_control_model = controler, senescence_model=senescence_model)

        # 1. Compare lesions before senescence
        # Compute variables of interest
        l = g.property('lesions')
        l = sum(list(l.values()), [])
        lesion1 = l[0]
        lesion1.compute_all_surfaces()
        lesion2 = l[1]
        lesion2.compute_all_surfaces()
        # Compare the two lesions
        assert lesion1.status == lesion2.status == 3
        assert lesion1.surface_alive == lesion2.surface_alive
        assert lesion1.surface_dead == lesion2.surface_dead == 0.
        assert lesion1.stock_spores == lesion2.stock_spores > 0.

        # 2. Set senescence, update lesion so they know they
        #    are on a senescent tissue and compare lesions
        set_properties(g, label = 'LeafElement', position_senescence=0.6)
        update(g, dt = 0., growth_control_model=controler, senescence_model=senescence_model)

        # Compute variables of interest
        l = g.property('lesions')
        l = sum(list(l.values()), [])
        lesion1 = l[0]
        lesion1.compute_all_surfaces()
        lesion2 = l[1]
        lesion2.compute_all_surfaces()
        # Compare the two lesions
        assert lesion2.growth_is_active == False
        assert lesion2.is_active == True
        assert lesion1.surface_alive > lesion2.surface_alive
        assert lesion2.surface_alive > 0.
        assert lesion1.surface_dead == 0.
        assert lesion2.surface_dead > 0.
        assert lesion2.surface > lesion2.surface_dead
        assert lesion1.stock_spores == lesion2.stock_spores > 0.

# if __name__ == '__main__':
    # g=test_growth_control()