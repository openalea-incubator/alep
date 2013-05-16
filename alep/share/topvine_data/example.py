from openalea.lpy import Lsystem
from openalea.mtg.io import lpy2mtg, mtg2lpy

from alinea.alep.inoculation import RandomInoculation
from alinea.alep.growth_control import NoPriorityGrowthControl
from alinea.alep.powdery_mildew import *
from alinea.alep.protocol import *


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


lpy_filename = "vinemtg_n_TT.lpy"
lsys = Lsystem(lpy_filename)

#creation de l'axiom
lstring0 = run(lsys,nbstep=8)
g0 = lpy2mtg(lstring0,lsys)
fungus = powdery_mildew()
PowderyMildewDU.fungus=fungus
nb_dus = 100
dispersal_units = ([PowderyMildewDU(nb_spores=1, status="emitted") for i in range(nb_dus)])
inoculator=RandomInoculation()
initiate(g0,dispersal_units,inoculator, label='lf')
axiom = mtg2lpy(g0,lsys)

#autres init
dt=1
controler = NoPriorityGrowthControl()
nsteps = 2

for i in range(nsteps):
    lsys.plot(axiom)
    lstring = run(lsys,axiom=axiom,nbstep=1)
    g = lpy2mtg(lstring,lsys)
    infect(g,dt)
    update(g,dt,controler)
    axiom = mtg2lpy(g,lsys)


