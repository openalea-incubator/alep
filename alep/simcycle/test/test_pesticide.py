from alinea.simcycle.pesticide import global_efficacy
from alinea.simcycle.Lesions import *
dose = {'Epoxiconazole':0.0125}

def test_call():
    eff = global_efficacy(dose)

def test_lesion():
    p=ParCycle()
    l=Lesions(p)
    l.addUdin()
    nh = 24
    T,PPFD, Rh = [20] * nh,[0]*24, [100]*24 
    sb= 1
    sa=l.devLes(sb,24,T,PPFD,Rh)
    print('Svert avant : %f, Svert apres: %f\n'%(sb,sa))
    
def test_lesion_with_pesticide():
    p=ParCycle()
    l=Lesions(p)
    l.addUdin()
    nh = 24
    T,PPFD, Rh = [20] * nh,[0]*24, [100]*24 
    sb= 1
    eff = global_efficacy(dose)
    sa=l.devLes(sb,24,T,PPFD,Rh,efficacy = eff)
    print('Svert avant : %f, Svert apres: %f\n'%(sb,sa))
    l.AfficheLesion()