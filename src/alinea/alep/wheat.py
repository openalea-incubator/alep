#alep/src/alinea/alep/wheat.py recupéré a partir du log ggarin authored and pradal committed on 2 Jun 2016 
#1 parent 446b694 commit 7a19f09a65e17e5271f53cb0d097a0ab8548f523

""" Few functions to call wheat MTGs. """

# Wheat stand ######################################################################
from alinea.adel.astk_interface import initialise_stand, AdelWheat
import alinea.adel.data_samples as adel_data
from alinea.adel.stand.stand import agronomicplot
from alinea.astk.plant_interface import *
from alinea.alep.architecture import update_healthy_area

def initialize_stand(age=0., length=0.1, width=0.2, sowing_density=150, 
                     plant_density=150, inter_row=0.12, nsect=1, seed = None, sample='random'):
    g, wheat, domain_area, domain, convUnit = initialise_stand(age = age, length=length, 
                                                            width=width, 
                                                            sowing_density=sowing_density, 
                                                            plant_density=plant_density,
                                                            inter_row=inter_row, 
                                                            nsect=nsect, seed= seed, sample=sample)
    # Add the property 'healthy_area' on the leaves
    update_healthy_area(g, label = 'LeafElement')
    return g, wheat, domain_area, domain
initialize_stand.__doc__ = initialise_stand.__doc__

def init_stand_height(h_factor=1., age=0., length=0.1, width=0.2, sowing_density=150, 
                      plant_density=150, inter_row=0.12, nsect=1, seed = None, sample='random'):
    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=length, 
                                                            width=width, 
                                                            sowing_density=sowing_density, 
                                                            plant_density=plant_density,
                                                            inter_row=inter_row)
    devT = adel_data.devT()
    devT['dimT']['L_internode']*=h_factor
    devT['dimT']['L_sheath']*=h_factor
    wheat = AdelWheat(nplants=nplants, positions = positions, nsect=nsect, devT=devT, seed=seed, sample=sample)
    g,_ = new_canopy(wheat,age=age)
    # Add the property 'healthy_area' on the leaves
    update_healthy_area(g, label = 'LeafElement')
    return g, wheat, domain_area, domain

# Mock-ups #########################################################################
from alinea.adel.newmtg import *
import alinea.adel.data_samples as adel_data
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import setAdel,RunAdel,genGeoLeaf,genGeoAxe
from math import ceil, sqrt
import random
import numpy

def adelR(nplants,dd):
    devT = adel_data.devT()
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants)
    cantable = RunAdel(dd,pars)
    return pars,cantable

def adel_mtg():
    """ create a very simple adel mtg """
    d = {'plant':[1,1],'axe_id':['MS','T1'],'ms_insertion':[0,1],'numphy':[1,1], 
         'Laz': [0,90], 'Ll' :[3,3], 'Lv' :[3,3] ,'Lr':[0,0], 'Lsen':[0,0], 'L_shape':[3,3], 'Lw_shape':[.3,.3], 'Linc':[0,0],
         'Einc':[0,45],'El':[1,1],'Ev':[1,1],'Esen':[0,0],'Ed': [0.1,0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=adel_data.leaves_db(), leaf_sectors=1)
    g=mtg_interpreter(g)
    return g

def adel_one_leaf():
    """ create a very simple adel mtg """
    d = {'plant':[1],'axe_id':['MS'],'ms_insertion':[0],'numphy':[1], 
         'Laz': [0], 'Ll' :[3], 'Lv' :[3] , 'Lr': [0], 'Lsen':[0], 'L_shape':[3], 'Lw_shape':[.3], 'Linc':[0],
         'Einc':[0],'El':[0],'Ev':[0],'Esen':[0],'Ed': [0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=adel_data.leaves_db(), leaf_sectors=1)
    g=mtg_interpreter(g)
    return g

def adel_one_leaf_element():
    """ create a very simple adel mtg """
    d = {'plant':[1],'axe_id':['MS'],'ms_insertion':[0],'numphy':[1], 
         'Laz': [0], 'Ll' :[3], 'Lv' :[3] , 'Lr': [0], 'Lsen':[0], 'L_shape':[3], 'Lw_shape':[.3], 'Linc':[0],
         'Einc':[0],'El':[0],'Ev':[0],'Esen':[0],'Ed': [0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=adel_data.leaves_db(), leaf_sectors=1)
    g=mtg_interpreter(g)
    g.remove_vertex(13)
    labels = g.property('label')
    labels[13] = 'Removed'
    return g

def adel_mtg2(nb_sect=1):
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=adel_data.leaves_db(),stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
    g=mtg_interpreter(g)
    return g

def adel_sample(nb_sect=1, d=None, p=None):
    """ create a less simple adel mtg """
    if p: # nb_plants
        if d:
            p, d = adelR(p, d)
        else:
            p, d = adelR(p, 1000)
        size = int(ceil(sqrt(len(p))))
        stand = numpy.array([(i, j) for i in range(size) for j in range(size)])
        numpy.random.shuffle(stand)
        stand = [((int(i)-10*size/2., int(j)-10*size/2., 0),random.randint(0,90)) for i, j in 10*stand[:len(p)]]
    else:
        size = 1.
        stand = [((0,0,0),0),((10,0,0),90), ((0,10,0), 0)]
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=adel_data.leaves_db(),stand=stand)
    g=mtg_interpreter(g)
    # approximation of domain area
    domain_area = size*0.05
    return g, domain_area

# Index finders ####################################################################
def find_blade_id(g, leaf_rank = 1, only_visible=True):
    labels = g.property('label')
    if only_visible==True:
        mets = [n for n in g if g.label(n).startswith('metamer') and g.order(n)==0
                and sum(component.visible_length for component in g.node(n).components())>0]
    else:
        mets = [n for n in g if g.label(n).startswith('metamer') and g.order(n)==0]
    bids = [co.index() for n in mets for co in g.node(n).components() if co.label.startswith('blade')]
    blade_id = bids[len(mets)-leaf_rank]
    return blade_id

def find_leaf_ids(g, blade_id):
    labels = g.property('label')
    leaf_elements = [id for id in g.components(blade_id) if labels[id].startswith('LeafElement')]
    return leaf_elements

# Tests ###########################################################################
def test_adel_mtg():
    """ Check the proper functioning of 'adel_mtg'.
    
    """
    g = adel_mtg()
    scene = plot3d(g)
    Viewer.display(scene)
    return g

def test_adel_one_leaf():
    """ Check the proper functioning of 'adel_one_leaf'.
    
    """
    g = adel_one_leaf()
    scene = plot3d(g)
    Viewer.display(scene)
    return g

def test_adel_mtg2():
    """ Check the proper functioning of 'adel_mtg2'.
    
    """
    g = adel_mtg2()
    scene = plot3d(g)
    Viewer.display(scene)
    return g

def test_adel_mtg3():
    """ Check the proper functioning of 'adel_mtg3'.
    
    """
    p, d = adelR(3,1000)
    g = adel_mtg3(d=d, p=p)
    scene = plot3d(g)
    Viewer.display(scene)
    return g