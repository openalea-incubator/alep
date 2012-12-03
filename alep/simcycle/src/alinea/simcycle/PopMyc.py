#  Python access to septo3D/PopMyc.h structure and functions


from ctypes import *
from alinea.simcycle.libsimcycle import *
lc = lib_simcycle


# Python classes for syntaxic sugar


class ClassMyc(cClassMyc):
    """
    A class of N lesion of age age whose surface is Q
    """

    def __repr__(self):
        s = c_char_p()
        s = lc.strClass(byref(self))
        return s

    def ResetCmyc(self):
        lc.ResetCmyc(byref(self))

    def expandClass(self,mult,Slmax,txsat):
        """ expand Q to Q * mult, controling for Q < Slmax and delta(Q) = txsat * Q * mult """
        return lc.expandClass(byref(self),c_float(mult),c_float(Slmax), c_float(txsat))
    
    def to_dict(self):
        d = dict([(f,getattr(self,f)) for f,_ in self._fields_])
        return d
        
class StackCmyc(cStackCmyc):
    """
    Priority stack of Cmyc Classes
    """

    def __init__(self,size=10):
        super(StackCmyc,self).__init__()
        lc.newStackCmyc_void(byref(self),int(size))

    def __del__(self):
        lc.freeStackCmyc(byref(self))

    def Ntot(self):
        """ Return the sum of N of all classes within the stack"""
        return lc.Ntot_sc(byref(self))

    def Qdisp(self):
        """ Sum of Q for classes of index >= ilim """
        return lc.Qdisp_sc(byref(self))

    def killInStack(self, S2kill):
        """ Try to kill classes with index > ilim up to S2kill. returns suface effectively killed """
        return lc.killInStack(byref(self),c_float(S2kill))

class StackBpop(cStackBpop):
    """
    Priority stack of Bpop Classes
    """

    def __init__(self,size=10):
        super(StackBpop,self).__init__()
        lc.newStackBpop_void(byref(self),int(size))

    def __del__(self):
        lc.freeStackBpop(byref(self))

    def Qdisp(self):
        """ Sum of Q for classes of index >= ilim """
        return lc.Qdisp_sbp(byref(self))

    def killInStack(self, S2kill):
        """ Try to kill classes with index > ilim up to S2kill. returns suface effectively killed """
        return lc.killInStackBp(byref(self),c_float(S2kill))


class PopMyc(cPopMyc):
    """
    A class for populations of Lesion
    """

    def __init__(self,agemin = 0, agemax = 100, largClass = 10,popType = 'Les'):

        if popType is 'Udin':
            pt = 1
        elif popType is 'Spo':
            pt = 2
        elif popType is 'Hlat':
            pt = 3
        else:
            pt = 0
        super(PopMyc,self).__init__()
        lc.newPopMyc_void(byref(self),c_float(agemin), c_float(agemax), c_float(largClass), c_int(pt))
        #pointer(self)[0] = lc.newPopMyc(c_float(agemin), c_float(agemax), c_float(largClass), c_int(pt))

    def __del__(self):
        """ free c pointers Q,N and decT before deletion"""
        lc.freePop(byref(self))

    def __repr__(self):
        s = c_char_p()
        s = lc.strPop(byref(self))
        return s

    def ResetPop(self):
        """ Reset Population Content """
        lc.ResetPop(byref(self))

    def ResetClass(self,numC):
        """ Reset Contents of class numC """
        lc.ResetClass(int(numC),byref(self))

    def ResetBpop(self,pos):
        """ Reset Contents of class in position pos """
        lc.ResetBpop(int(pos),byref(self))

    def posC(self, numC):
        """ return index of Class numC"""
        return lc.posC(int(numC),byref(self))

    def numCAge(self, age):
        """ return Class number containing lesions of age age """
        return lc.numCAge(c_float(age),byref(self))
    
    def borneInf(self,numC):
        """ returns inferior boundary of class numC  """
        return lc.borneInf(int(numC),byref(self))

    def decTAge(self,age):
        """ return Class decT containing lesions of a given age """
        return lc.decTAge(c_float(age),byref(self))

    def ageC(self, numC):
        """ return age of Class numC """
        return lc.ageC(int(numC),byref(self))

    def bornemax(self):
        """ returns superior boundary of the pop """
        return lc.bornemax(byref(self))

    def getCmyc(self, numC):
        """ return ClassMyc numC  """
        c = ClassMyc()
        lc.getCmyc_void(pointer(c),numC,byref(self))
        return c

    def getBarPop(self, numC):
        """ return BarPop numC """
        bp = BarPop()
        lc.getBarPop_void(pointer(bp),numC,byref(self))
        return bpc

    def PutBarPop(self, bp):
        """ Put BarPop bp in pop """
        lc.PutBarPop(bp,byref(self))

    def Qtot(self):
        """ returns sum of Q  """
        return lc.Qtot(byref(self))

    def Ntot(self):
        """ returns sum of N  """
        return lc.Ntot(byref(self))

    def NQtot(self):
        """ returns sum of  NQ products of pop """
        return lc.NQtot(byref(self))

    def NFilledC(self):
        """ Count of non-empty classes within the pop """
        return lc.NFilledC(byref(self))
        
    def getQ(self):
        """ return the list of Q of all cohorts, sorted by class index """
        Q = (c_float * self.nbClass)()
        lc.getQ(byref(self),pointer(Q))
        return [val for val in Q]

    def getN(self):
        """ return the list of N of all cohorts, sorted by class index """
        N = (c_float * self.nbClass)()
        lc.getN(byref(self),pointer(N))
        return [val for val in N]
        
    def getdecT(self):
        """ return the list of decT of all cohorts, sorted by class index """
        decT = (c_float * self.nbClass)()
        lc.getdecT(byref(self),pointer(decT))
        return [val for val in decT]

    def NtotEntre(self,agemin,agemax):
        """ returns total N inluded beteween agemin and agemeax  """
        return lc.NtotEntre(byref(self),c_float(agemin),c_float(agemax))

    def QtotEntre(self,agemin,agemax):
        """ returns total Q inluded beteween agemin and agemeax """
        return lc.QtotEntre(byref(self),c_float(agemin),c_float(agemax))

    def NQtotEntre(self,agemin,agemax):
        """ returns total NQ inluded beteween agemin and agemeax """
        return lc.NQtotEntre(byref(self),c_float(agemin),c_float(agemax))

    def PopAsStack(self):
        """ return a priority stack (older first) of classes within the pop """
        s = StackCmyc()
        lc.PopAsStack_void(pointer(s),byref(self))
        return s

    def PopAsStack_bp(self):
        """ return a priority stack (older first) of classes within the pop """
        s = cStackBpop()
        lc.PopAsStack_bp_void(pointer(s),byref(self))
        return s
        
    def getType(self):
        """ return the type of the Pop """
        t = (c_char * 20)()
        lc.getType(pointer(t),byref(self))
        return t.value

    def PrintPop(self):
        """ Print details on pop contents """
        lc.PrintPop(self)

    def AddCtoPop(self,Cmyc):
        """ Add Class Cmyc to the pop, returns the class number """
        lc.AddCtoPop(Cmyc, byref(self))
                      
    def RemoveQtoPop(self,Q):
        """ Decreases pop surface by Q up to zero, returns effective surface removed. Class surface reduction is proportional to class Q / Qtot"""
        lc.RemoveQtoPop(c_float(Q), byref(self))

    def RemoveNtoPop(self,S):
        """ Decreases pop surface by S up to zero, by decreasing lesions number in each class so that surface decrease by class is the same. Returns effective surface removed """
        lc.RemoveNtoPop(c_float(Q), byref(self))

    def TireCbyN(self, debug = 0):
        """ Return a random Class number, the propability for a class to be retruned beeing proportional to its lesion ammount (N) """
        lc.TireCbyNd(byref(self), int(debug))

    def RemoveQbyN(self, Q, pKillLast = 0, debug = 0):
        """ Decreases pop surface by Q up to zero, returns effective surface removed. Class surface reduction is proportional to class N / Ntot.pKillLast is the probability of killing the lesion that makes the rounding """
        lc.RemoveQbyNd(c_float(Q), pointer(self), c_float(pKillLast), int(debug))

    def ClassEqQ(self):
        """ returns a ClassMyc with population averaged age, Total Q and N ammount """
        Cmyc = ClassMyc()
        lc.ClassEqQ_void(pointer(Cmyc), byref(self))
        return Cmyc

    def AgePop(self, dt):
        """ Update population so that age = age + dt. Returns the population of classes leaving the pop """
        popout = cPopMyc()
        lc.AgePop_void(pointer(popout), pointer(self), c_float(dt))
        return popout

    def calcdemC(self, mult, Slmax):
        """ Returns the surface needed by a pop if the surface of all its class are multiplied by mult, but limited to Slmax by lesion """
        return lc.calcdemC(byref(self), c_float(mult), c_float(Slmax))

    def expandPop(self, multmax, Slmax, txsat):
        """ Expand the surface of the pop by multiplying the surface of all its classes by multmax and adjusting to txsat. The surface of individual lesions is limited to Slmax; Returns the newly created surface """
        lc.expandPop(byref(self),c_float(multmax), c_float(Slmax), c_float(txsat))
    
    def to_dict(self):
        d = dict([(f,getattr(self,f)) for f,_ in self._fields_])
        d['Q'] = [self.Q[i] for i in range(self.nbClass)]
        d['N'] = [self.N[i] for i in range(self.nbClass)]
        d['decT'] = [self.decT[i] for i in range(self.nbClass)]
        return d


# Functions outside classes

def AddCtoPops(C,Nbis,Qbis,pop,popbis):
    """ Add a class in pop and in popbis with modified contents """
    return lc.AddCtoPops(C,c_float(Nbis),c_float(Qbis),byref(pop),byref(popbis))

def propEntiere(N,p):
    """ Rounds N*p with probability """
    return lc.propEntiere(c_float(N),c_float(p))

# Utilites for interactive python wrapping of cPopMyc objects

def PopMyc_fromC(PopMyc_C):
    """
    Return PopMyc python class instance from cPopMyc object
    """
    p = PopMyc(PopMyc_C.agemin,PopMyc_C.agemax,PopMyc_C.largClass,PopMyc_C.type)
    lc.copyPopMyc(byref(PopMyc_C),byref(p))
    return p

##def newPopMyc(agemin = 0, agemax = 100, largClass = 10.):
##    pop = PopMyc()
##    n = int(ceil(float(agemax - agemin) / largClass))
##    if n < 1:
##        n = 1
##    pop.nbClass = n 
##    pop.agemin = agemin
##    pop.agemax = agemax
##    # astuce pour les pointer sur des tableaux
##    pop.Q = cast((c_float * n)(), POINTER(c_float))
##    pop.N = cast((c_float * n)(), POINTER(c_float))
##    pop.decT = cast((c_float * n)(), POINTER(c_float))
##    pop.largClass = largClass
##    pop.pC0 = 0
##    pop.TimeBuffer = 0
##    return pop



# functions
                      
# Factories : doesn't work on windows where, contrary to linux, restype is used for checking (error in calling convention): does that mean that only simple ctypes are allowed for restype ?
##lc.newCmyc_p.restype = cClassMyc
##lc.newPopMyc.restype = PopMyc
##lc.getCmyc.restype = ClassMyc
##lc.getBarPop.restype = BarPop
##lc.ClassEqQ.restype = ClassMyc
##lc.AgePop.restype = PopMyc

