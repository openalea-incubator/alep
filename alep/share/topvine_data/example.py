from openalea.lpy import Lsystem,AxialTree,generateScene
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
    return lsystem.iterate(axiom,c_iter,nbstep), lsystem

def iterate(lsystem, axiom = '', nbstep = -1, parameters = {}, local = {}):
    if parameters and (not isinstance(parameters, dict)) and is_file(parameters):
        with open(parameters) as f:
            py_code = f.read()
            parameters = {}
            exec(py_code, globals(), parameters)

    local_settings = dict(parameters)
    if local:
        local_settings.update(local)

    axial, l = run(lsystem, axiom, nbstep, local_settings)
    params = dict(parameters)
    d = dict()
    l.context().getNamespace(d)
    for k in parameters:
        params[k] = d[k]
    return axial, l, params

lpy_filename = "vinemtg_n_TT.lpy"
l = Lsystem(lpy_filename)

#creation de l'axiom
lstring,l,pars= iterate(l,nbstep=8)
g = lpy2mtg(lstring,l)
fungus = powdery_mildew()
PowderyMildewDU.fungus=fungus
nb_dus = 100
dispersal_units = ([PowderyMildewDU(nb_spores=1, status="emitted") for i in range(nb_dus)])
inoculator=RandomInoculation()
initiate(g,dispersal_units,inoculator, label='lf')
axiom = mtg2lpy(g,l)

#autres init
# dt=1
# controler = NoPriorityGrowthControl()
# nsteps = 2

# for i in range(nsteps):
    # lstring,l,params = iterate(l,axiom=axiom,nbstep=1)
    # g = lpy2mtg(lstring,l)
    # infect(g,dt)
    # update(g,dt,controler)
    # axiom = mtg2lpy(g)


