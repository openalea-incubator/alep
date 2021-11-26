from alinea.alep.septoria import *
import random as rd

def variable_septoria(distri_chlorosis = {'mu':200., 'sigma':30.}):
    assert (isinstance(distri_chlorosis, dict) and 
            list(distri_chlorosis.keys()) == ['mu', 'sigma']), ("distri_chlorosis must be in the form {'mu':200., 'sigma':30.}")
    rnd = rd.Random(1)
    septoria_parameters['mu'] = distri_chlorosis['mu']
    septoria_parameters['sigma'] = distri_chlorosis['sigma']

    class VariableSeptoria(SeptoriaLesion):
        def __init__(self, mutable = True):
            super(VariableSeptoria, self).__init__(mutable = mutable)
            self.fungus.degree_days_to_chlorosis = rnd.gauss(mu=self.fungus.mu,
                                                             sigma=self.fungus.sigma)

    class VariableSeptoriaFungus(Fungus):
        def __init__(self, Lesion=VariableSeptoria, 
                        DispersalUnit=SeptoriaDU, parameters=septoria_parameters):
            super(VariableSeptoriaFungus, self).__init__(Lesion=Lesion,
                                                         DispersalUnit = DispersalUnit, 
                                                         parameters=parameters)
                                                         
    return VariableSeptoriaFungus()
