"""
implements a vine plant model based on topvine lsystem
"""

from openalea.lpy import Lsystem,AxialTree, generateScene
from openalea.mtg.io import lpy2mtg, mtg2lpy

import os
vinedir = os.path.dirname(__file__)

def run(lsystem, axiom = '', nbstep = -1, parameters = {}):
    """ Run a lsystem """
    c_iter = lsystem.getLastIterationNb()
    if nbstep < 0:
        nbstep = lsystem.derivationLength - c_iter
    if len(axiom) == 0:
        axiom = lsystem.axiom
    elif type(axiom) == str:
        lsystem.makeCurrent()
        axiom = AxialTree(axiom)
        lsystem.done()
    if len(parameters) > 0:
        lsystem.context().updateNamespace(parameters)
    return lsystem.iterate(axiom,c_iter,nbstep)
    
class Vine(object):
    
    def __init__(self, lpy_filename = vinedir + '/topvine.lpy'):
        self.lpy_filename = lpy_filename
        inter_row = 2.2# m
        inter_plant = 1.1
        self.domain_area = 3 * inter_plant * inter_row
        
    def setup_canopy(self, age=0):
        self.start = age
        self.lsys = Lsystem(self.lpy_filename)
        tree = run(self.lsys, nbstep = int(2 + age))
        g = lpy2mtg(tree,self.lsys)
        return g
        
    def grow(self,g,time_control):
        if len(time_control) > 0:
            axiom = mtg2lpy(g,self.lsys)
            #time_control.check('dt',1)
            # /!\ TEMP /!\ ######################################
            dt = 1  
            tree = run(self.lsys,axiom = axiom, nbstep = dt)
            # /!\ TEMP /!\ ######################################
            # tree = run(self.lsys,axiom = axiom, nbstep = time_control.dt)
            return lpy2mtg(tree,self.lsys)
        else :
            return g
        
    def reset(self):
        if self.start is None:
            self.start = 0
        return setup_canopy(self.start)
        
    def plot(self,g):
        tree = mtg2lpy(g,self.lsys)
        self.lsys.plot(tree)
        return g
        
    def generate_scene(self,g):
        tree = mtg2lpy(g,self.lsys)
        self.lsys.sceneInterpretation(tree)
        return generateScene(tree)
        

    