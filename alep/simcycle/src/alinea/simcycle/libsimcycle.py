#  Wrapping of Septo3d functions for simulating deaseases


import os
from ctypes import *
from ctypes.util import find_library

lib_path = 'd:\\openaleapkg\\simcycle\\build-scons\\lib\\simcycle.dll'
if not os.path.exists(lib_path):
    lib_path = find_library('simcycle')
if lib_path is None:
    if os.name == 'posix':
        lib_path = 'libsimcycle.so'
lib_simcycle = cdll.LoadLibrary(lib_path)

# declaration of cStructs : be carefull !  order and type of parameter matters a lot!

# PopMyc

class cClassMyc(Structure):
    """
    A class of N lesion of age age whose surface is Q
    """
    _fields_ = [('Q', c_float),
                ('age',c_float),
                ('N',c_float)]

class cBarPop(Structure):
    """
    """
    _fields_ = [('Q', c_float),
                ('decT', c_float),
                ('pos', c_int),
                ('N', c_float)]

class cStackCmyc(Structure):
    """
    Priority stack of Cmyc Classes
    """
    _fields_ = [('stack', POINTER(cClassMyc)),
                ('nclass', c_int),
                ('ilim', c_int)]

class cStackBpop(Structure):
    """
    Priority stack of Bpop Classes
    """
    _fields_ = [('stack', POINTER(cBarPop)),
                ('nclass', c_int),
                ('ilim', c_int)]

class cPopMyc(Structure):
    """
    """
    _fields_ = [('type', c_int),
                ('nbClass', c_int),
                ('agemin', c_float),
                ('agemax', c_float),
                ('Q', POINTER(c_float)),
                ('N', POINTER(c_float)),
                ('decT', POINTER(c_float)),
                ('largClass', c_float),
                ('pC0', c_int),
                ('TimeBuffer', c_float)]

# Lesions

nlesrec = lib_simcycle.getnLesRec()

class cParCycle(Structure):
    """
    """
    _fields_ = [('txPerteUdin',c_float),
                ('HrGerm',c_float),
                ('PARgerm',c_float),
                ('DminGerm',c_float),
                ('pinc',c_float),
                ('TbasePath',c_float),
                ('TmaxPath',c_float),
                ('Lcdd',c_float),
                ('FixedLat',c_int),
                ('DInc',c_float),
                ('del2',c_float),
                ('fdhr',c_float),
                ('SIncMin',c_float),
                ('txInc',c_float),
                ('DChlo',c_float),
                ('Slmin',c_float),
                ('Slmax',c_float),
                ('txCroi',c_float),
                ('nbEvSplaMax',c_int),
                ('pNoSpo',POINTER(c_float)),
                ('pEv', c_double),
                ('dstom',c_float),
                ('pStoInf',c_float),
                ('NbTotSpore',c_int),
                ('pSporeRes',c_float),
                ('pSpoEv',POINTER(c_float)),
                ('LaiSpoSol',c_float),
                ('TTinfSol',c_float)]

class cLesions(Structure):
    """
    """
    _fields_ = [('NewUdin',cClassMyc),
                ('PopUdin',cPopMyc),
                ('PopInc',cPopMyc),
                ('PopHlat',cPopMyc),
                ('PopCroi',cPopMyc),
                ('PopMat',cPopMyc),
                ('PopChlo',cPopMyc),
                ('PopSpo',cPopMyc),
                ('nbLesAd',c_int),
                ('SEmpLes',c_float),
                ('SLesInSen',c_float),
                ('nhGerm',c_int),
                ('dTTmax',c_float),
                ('skipInc',c_int),
                ('skipChlo',c_int),
                ('skipCroi',c_int),
                ('nLesRec',c_int),
                ('LesRec',(c_float * nlesrec)),
                ('par',cParCycle)]

# declaration of simple restypes (other than int or structures)


# PopMyc

lib_simcycle.borneInf.restype = c_float
lib_simcycle.decTAge.restype = c_float
lib_simcycle.ageC.restype = c_float
lib_simcycle.bornemax.restype = c_float
lib_simcycle.Qtot.restype = c_float
lib_simcycle.Ntot.restype = c_float
lib_simcycle.Ntot_sc.restype = c_float
lib_simcycle.Qdisp_sc.restype = c_float
lib_simcycle.Qdisp_sbp.restype = c_float
lib_simcycle.NQtot.restype = c_float
lib_simcycle.NtotEntre.restype = c_float
lib_simcycle.QtotEntre.restype = c_float
lib_simcycle.NQtotEntre.restype = c_float
lib_simcycle.strClass.restype = c_char_p
lib_simcycle.RemoveQtoPop.restype = c_float
lib_simcycle.RemoveNtoPop.restype = c_float
lib_simcycle.propEntiere.restype = c_float
lib_simcycle.RemoveQbyNd.restype = c_float
lib_simcycle.killInStack.restype = c_float
lib_simcycle.killInStackBp.restype = c_float
lib_simcycle.expandClass.restype = c_float
lib_simcycle.calcdemC.restype = c_float
lib_simcycle.expandPop.restype = c_float
lib_simcycle.strPop.restype = c_char_p

#Lesions

lib_simcycle.getLesRec.restype = c_float
lib_simcycle.SLesInc.restype = c_float
lib_simcycle.NLesInc.restype = c_float
lib_simcycle.SLesSpo.restype = c_float
lib_simcycle.SLesNec.restype = c_float
lib_simcycle.SLesChlo.restype = c_float
lib_simcycle.NLesChloNec.restype = c_float
lib_simcycle.NLesNec.restype = c_float
lib_simcycle.SLesTot.restype = c_float
lib_simcycle.NLesTot.restype = c_float
lib_simcycle.SLesOnGreen.restype = c_float
lib_simcycle.Lesions_somT.restype = c_float
lib_simcycle.Dlat.restype = c_double
lib_simcycle.TxSat.restype = c_float
lib_simcycle.fertUdin.restype = c_float
lib_simcycle.freeSpaceIncC.restype = c_float
lib_simcycle.dSChlo.restype = c_float
lib_simcycle.freeSpace.restype = c_float
lib_simcycle.freeSpaceNewChlo.restype = c_float
lib_simcycle.getEclinSoil.restype = c_double

