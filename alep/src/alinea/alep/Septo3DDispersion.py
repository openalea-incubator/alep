"""
 fake class for septo3D dispersion model
"""




class Septoria(object):
    
    def emission(leaf):
       
        """change state of rings"""
        return NSpores
        
        
class SplashDispersal(object):
    
    from alinea.Ldisp.Dispersor import *
    
    def __init__(RainItensity,RainDuration)
        self.Ip = RainItensity
        self.Dp=RainDuration
        self.parD = parDisp()
        disp = newDispersor()
        pass
        
    
        
    def disperse(g, dispersal_units):
        
        ilayers = [newILayerts(leaf) for leaf in g]
        
        self.disp.newDispersion(g)
        
        self.disp.DispInterception(self.Ip)
        self.disp.DispEclas(self.Ip,parD)
        
        for d in dispersal_units:
            leaf = g.node(d.vid)
            les = leaf.l
            
            eclas = self.disp.getQL(Sleaf,ilayers[d.vid]);
            eclasSpo = eclas * les.Sspo# to be in lesion class 
            eclins = eclas * (1 -peV-PnoSpo)
            Qspores = dispersal_unit.NSpores / eclins
            self.disp.addEclin(elcin, QSpores)
            DispSplash(Ip,Dp)
            Udin = getQL(Sleaf[j], &(Disp.Udin),&(Disp.LAI), ISect[j]);
            #here compute relatiove position
            deposits[d.vid].NSpores = udin.NSpores