
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
                nodemodule='alinea.alep.protocol',
                nodeclass='initiate',
               )
__all__.append('protocol_nodes_initiate')

protocol_nodes_infect = Factory(name='infect',
                nodemodule='alinea.alep.protocol',
                nodeclass='infect',
               )
__all__.append('protocol_nodes_infect')

protocol_nodes_update = Factory(name='update',
                nodemodule='alinea.alep.protocol',
                nodeclass='update',
               )
__all__.append('protocol_nodes_update')

protocol_nodes_disperse = Factory(name='disperse',
                nodemodule='alinea.alep.protocol',
                nodeclass='disperse',
                inputs=[{'interface': None, 'name': 'g', 'value': None},
                        {'interface': None, 'name': 'emission model', 'value': None, 'desc': 'None if emission is handled by lesions'},
                        {'interface': None, 'name': 'transport model', 'value': None, 'desc': ''},
                        {'interface': IStr, 'name': 'fungus_name', 'value': 'septoria'},
                        {'interface': IStr, 'name': 'label', 'value': 'LeafElement'},
                        {'interface': IBool, 'name': 'activate', 'value': True},
                        ],
                outputs=[{'interface': None, 'name': 'g', 'desc': ''},
                         {'interface': IFloat, 'name': 'nb_dus', 'desc': ''},
                ],
               )
__all__.append('protocol_nodes_disperse')

protocol_nodes_wash = Factory(name='wash',
                nodemodule='alinea.alep.protocol',
                nodeclass='wash',
               )
__all__.append('protocol_nodes_wash')

protocol_nodes_contamination = Factory(name='external contamination',
                nodemodule='alinea.alep.protocol',
                nodeclass='external_contamination',
               )
__all__.append('protocol_nodes_contamination')
