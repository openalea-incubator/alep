# -*- python -*-
#
#       container.mesh: node for container.topomesh
#
#       Copyright 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#

__doc__="""
node definition for container.mesh package
"""

__license__= "Cecill-C"
__revision__=" $Id: __wralea__.py 7865 2010-02-08 18:04:39Z cokelaer $ "

from openalea.core import Factory
from openalea.core.interface import *

__name__ = "openalea.container.mesh"
__alias__ = ['container.mesh']
__version__ = '0.0.3'
__license__ = "Cecill-C"
__authors__ = 'Jerome Chopard'
__institutes__ = 'INRIA/CIRAD'
__description__ = 'container.mesh Node library.'
__url__ = 'http://openalea.gforge.inria.fr'

__all__ = []

#########################################
#
#       read
#
#########################################
read = Factory( name= "read",
                                description= "",
                                category = "",
                                nodemodule = "mesh",
                                nodeclass = "read",
                                inputs=(dict(name="filename", interface=IFileStr,),),
                                outputs=(dict(name="mesh", interface=None,),
                                                dict(name="mesh_prop", interface=None,),),
                        )

__all__.append('read')

#########################################
#
#       edition
#
#########################################
remove_wisp = Factory( name= "remove_wisp",
                                description= "",
                                category = "",
                                nodemodule = "mesh",
                                nodeclass = "remove_wisp",
                                inputs=(dict(name="mesh", interface=None,),
                                                dict(name="scale", interface=IInt,),
                                                dict(name="wid", interface=IInt),),
                                outputs=(dict(name="mesh", interface=None,),),
                        )

__all__.append('remove_wisp')

add_wisp = Factory( name= "add_wisp",
                                description= "",
                                category = "",
                                nodemodule = "mesh",
                                nodeclass = "add_wisp",
                                inputs=(dict(name="mesh", interface=None,),
                                                dict(name="scale", interface=IInt,),
                                                dict(name="wid", interface=IInt, value=None),),
                                outputs=(dict(name="mesh", interface=None,),
                                                dict(name="wid", interface=IInt,),),
                        )

__all__.append('add_wisp')

#########################################
#
#       cleaning
#
#########################################
clean_geometry = Factory( name= "clean_geometry",
                                description= "",
                                category = "",
                                nodemodule = "mesh",
                                nodeclass = "clean_geometry",
                                inputs=(dict(name="mesh", interface=None,),),
                                outputs=(dict(name="mesh", interface=None,),),
                        )

__all__.append('clean_geometry')

clean_orphans = Factory( name= "clean_orphans",
                                description= "",
                                category = "",
                                nodemodule = "mesh",
                                nodeclass = "clean_orphans",
                                inputs=(dict(name="mesh", interface=None,),),
                                outputs=(dict(name="mesh", interface=None,),),
                        )

__all__.append('clean_orphans')
