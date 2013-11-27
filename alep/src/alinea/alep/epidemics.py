""" Definitions of assemby of modules/processes that define completely an epidemics interacting with a plant
"""

class SeptoriaSepto3D(object):
    from alinea.alep.septoria import plugin_septoria
    from alinea.alep.inoculation import RandomInoculation
    from alinea.alep.growth_control import NoPriorityGrowthControl
    from alinea.alep.infection_control import BiotrophDUPositionModel
    from alinea.alep.senescence import WheatSeptoriaPositionedSenescence
    from alinea.alep.dispersal_emission import SeptoriaRainEmission
    from alinea.septo3d.dispersion.alep_interfaces import Septo3DTransport
    from alinea.alep.washing import RapillyWashing

   
    def __init__(g,label='LeafElement', domain, domain_area, convUnit):
    #a priori si on considere ces attributs comme generique, plutot viser un ininit de ces champs sans argument, et faire le setting ensuite. Ceci afin d'etre compatible avec plugin (qui ne remplace que l'etape import)
        self.disease = plugin_septoria()
        self.inoculator = RandomInoculation()
        self.growth_controler = NoPriorityGrowthControl()
        self.infection_controler = BiotrophDUPositionModel()
        self.sen_model = WheatSeptoriaPositionedSenescence(g, label='LeafElement')
        self.emitter = SeptoriaRainEmission()
        self.transporter = Septo3DTransport(domain=domain, domain_area=domain_area, convUnit = convUnit)
        self.washor = RapillyWashing()
    