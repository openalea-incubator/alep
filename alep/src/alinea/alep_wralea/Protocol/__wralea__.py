
# This file has been generated at Wed Feb 20 15:02:32 2013

from openalea.core import *


__name__ = 'Alep.Protocol'

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
        
protocol_nodes_initiate = Factory(name='initiate',
                nodemodule='protocol_nodes',
                nodeclass='initiate',
               )
__all__.append('protocol_nodes_initiate')

protocol_nodes_infect = Factory(name='infect',
                nodemodule='protocol_nodes',
                nodeclass='infect',
               )
__all__.append('protocol_nodes_infect')

protocol_nodes_update = Factory(name='update',
                nodemodule='protocol_nodes',
                nodeclass='update',
               )
__all__.append('protocol_nodes_update')

protocol_nodes_disperse = Factory(name='disperse',
                nodemodule='protocol_nodes',
                nodeclass='disperse',
               )
__all__.append('protocol_nodes_disperse')

protocol_nodes_wash = Factory(name='wash',
                nodemodule='protocol_nodes',
                nodeclass='wash',
               )
__all__.append('protocol_nodes_wash')

protocol_nodes_growth_control = Factory(name='growth_control',
                nodemodule='protocol_nodes',
                nodeclass='growth_control',
               )
__all__.append('protocol_nodes_growth_control')


