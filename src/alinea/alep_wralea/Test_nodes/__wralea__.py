
# This file has been generated at Wed Feb 20 15:02:32 2013

from openalea.core import *


__name__ = 'Alep.Test_nodes'

__editable__ = True
__description__ = ''
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '0.8.0'
__authors__ = ''
__institutes__ = None
__icon__ = ''


__all__ = []

test_protocol_nodes_wheat = Factory(name='wheat',
                                    nodemodule='alinea.alep_wralea.Test_nodes.test_protocol_nodes',
                                    nodeclass='wheat',
                                    )
__all__.append('test_protocol_nodes_wheat')

test_protocol_nodes_microclimate = Factory(name='microclimate',
                                    nodemodule='alinea.alep_wralea.Test_nodes.test_protocol_nodes',
                                    nodeclass='microclimate',
                                    )
__all__.append('test_protocol_nodes_microclimate')

test_protocol_nodes_weather_reader = Factory(name='weather_reader',
                                    nodemodule='alinea.alep_wralea.Test_nodes.test_protocol_nodes',
                                    nodeclass='weather_reader',
                                    )
__all__.append('test_protocol_nodes_weather_reader')

test_protocol_nodes_scene_from_g = Factory(name='scene_from_g',
                                           nodemodule='alinea.alep_wralea.Test_nodes.test_protocol_nodes',
                                           nodeclass='scene_from_g',
                                           )
__all__.append('test_protocol_nodes_scene_from_g')

test_protocol_nodes_set_properties_node = Factory(name='set_properties_node',
                                             nodemodule='alinea.alep_wralea.Test_nodes.test_protocol_nodes',
                                             nodeclass='set_properties_node',
                                             )
__all__.append('test_protocol_nodes_set_properties_node')

test_protocol_nodes_distribute_dispersal_units = Factory(name='distribute_dispersal_units',
                                                         nodemodule='alinea.alep_wralea.Test_nodes.test_protocol_nodes',
                                                         nodeclass='distribute_dispersal_units',
                                                         )
__all__.append('test_protocol_nodes_distribute_dispersal_units')

test_protocol_nodes_distribute_lesions = Factory(name='distribute_lesions',
                                                 nodemodule='alinea.alep_wralea.Test_nodes.test_protocol_nodes',
                                                 nodeclass='distribute_lesions',
                                                 )
__all__.append('test_protocol_nodes_distribute_lesions')




