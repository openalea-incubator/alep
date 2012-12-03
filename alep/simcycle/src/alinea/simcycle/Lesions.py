#  Python access to septo3D/Lesions.h structure and functions


from ctypes import *
from alinea.simcycle.PopMyc import *
from alinea.simcycle.libsimcycle import *
lc = lib_simcycle

# Python classes for syntaxic sugar

class ParCycle(cParCycle):
    """
    Parameter list for Septoria as found in Septo3D RefParameter.h
    """

    def __init__(self):
        super(ParCycle,self).__init__()
        lc.getParCycle_void(byref(self))

    def __del__(self):
        lc.freeParCycle(byref(self))

    def to_dict(self):
        d = dict([(f,getattr(self,f)) for f,_ in self._fields_])
        d['pNoSpo'] = [self.pNoSpo[i] for i in range(self.nbEvSplaMax)]
        d['pSpoEv'] = [self.pSpoEv[i] for i in range(self.nbEvSplaMax)]
        return d

    def getEclinSoil(self, TTem, EclaSol):
        """
        Returns infectious droplets emmited by the soil at date TTem
        """
        lc.getEclinSoil(c_float(TTem), c_double(EclaSol), byref(self))

def ParCycle_fromC(ParCycle_C):
    p = ParCycle()
    lc.copyParCycle(byref(ParCycle_C),byref(p))
    return(p)
        
        
class Lesions(cLesions):
    """
    A class for simulating the infectious cycle of foliar pathogen producing lesions
    """

    def __init__(self,par):
        super(Lesions,self).__init__()
        lc.newLes_void(byref(self),byref(par))

    def __del__(self):
        lc.freeLes(byref(self))

    def to_dict(self):
        d = dict([(f,getattr(self,f)) for f,_ in self._fields_])
        d['LesRes'] = [v for v in self.LesRec]
        d['NewUdin'] = dict([(f,getattr(self.NewUdin,f)) for f,_ in self.NewUdin._fields_])
        d['PopUdin'] = PopMyc_fromC(self.PopUdin).to_dict()
        d['PopInc'] = PopMyc_fromC(self.PopInc).to_dict()
        d['PopHlat'] = PopMyc_fromC(self.PopHlat).to_dict()
        d['PopCroi'] = PopMyc_fromC(self.PopCroi).to_dict()
        d['PopMat'] = PopMyc_fromC(self.PopMat).to_dict()
        d['PopChlo'] = PopMyc_fromC(self.PopChlo).to_dict()
        d['PopSpo'] = PopMyc_fromC(self.PopSpo).to_dict()
        d['par'] = ParCycle_fromC(self.par).to_dict()
        return d
        
    def resetLesrec(self):
        """
        Reset records in LesRec
        """
        lc.resetLesRec(byref(self))

    def ResetLes(self):
        """
        Reset Lesion to initial state
        """
        lc.ResetLes(byref(self))

    def addUdin(self, N=1, Q=1, age=0):
        """ajout de nouvelles Udins dans NewUdin, lors de l'interception.age = delai(jour) entre le debut du pas de temps et le depot, N = Nombres d'udin, Q = nombres de spores
        """
        lc.addUdin(byref(self),c_float(N), c_float(Q), c_float(age))
        
    def addCtoChlo(self, Cmyc):
        """
        Add a cohorts of incubating lesions to Chlo.
        age = agedd since transition to Chlo
        N = number of lesions
        Q = total surface of the cohort
        """
        lc.addCtoChlo(byref(self),Cmyc)

    def AfficheLesion(self):
        """
        Print Lesion detailed state
        """
        lc.AfficheLesion(byref(self))

    def AfficheLesRec(self):
        """
        Print Current records of the lesion
        """
        lc.AfficheLesRec(byref(self))

    def devLes(self, Svert, dt, T, PPFD, Rh, efficacy = {'protectant' : 0, 'eradicant' : 0}, debug = 0):
        """
        Simulate Lesion developpement  on a Green area Svert during dt
        Return the updated  green Surface
        """
        assert len(T) >= dt and len(PPFD) >= dt and len(Rh) >= dt
        tab = c_float * dt
        Sv = c_float(Svert)
        lc.devLes_dp(byref(self),byref(Sv),int(dt),tab(*T[0:dt]),tab(*PPFD[0:dt]),tab(*Rh[0:dt]),c_float(efficacy['protectant']), c_float(efficacy['eradicant']), int(debug))
        return Sv.value
        
    def senLes(self, StoRec, Svert, Sdsen):
        """
        Update Lesions during senescence. StoRec is the Lesion surface to be replaced by senescent area.
        return update of Green area and natural apical senescence area.
        Seenescene transforms (kills/replaces) chlorotic and incubating area into natural senescence
        Necrotic area 'absorb' part of StoRec (as if they are recovered) but remain necrotic (simulates same rate of progression of apical senescence on leaves as in controled conditions)
        """
        Sv = c_float(Svert)
        Sd = c_float(Sdsen)
        lc.senLes(c_float(StoRec), byref(self), byref(Sv), byref(Sd))
        return Sv.value, Sd.value
