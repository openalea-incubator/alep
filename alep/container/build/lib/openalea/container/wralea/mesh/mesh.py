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
__revision__=" $Id: mesh.py 7865 2010-02-08 18:04:39Z cokelaer $ "

from openalea.vmanalysis.serial.mesh_read import read
from openalea.container import topomesh_algo as algo

#################################################
#
#       read write
#
#################################################
def read_mesh (filename) :#TODO find a way to define properties more conveniently
    """
    read a mesh from a file and return
    a tuple mesh,mesh_prop
    """
    print "read mesh"
    return read(filename)

#################################################
#
#       basic edition
#
#################################################
def remove_wisp (mesh, scale, wid) :
    """
    remove a wisp from the mesh
    and return it
    scale: scale of the wisp
    wid: id of the wisp to remove
    """
    mesh.remove_wisp(scale,wid)
    return mesh,

def add_wisp (mesh, scale, wid) :
    """
    add a wisp in the mesh
    and return the mesh and the id of the added wisp
    scale: scale of the wisp
    wid: id of the wisp to add (may be None)
    """
    wid = mesh.add_wisp(scale,wid)
    return mesh,wid
#################################################
#
#       cleaning edition
#
#################################################
def clean_geometry (mesh) :
    """
    remove wisps not geometrically defined
    i.e. with no borders
    """
    algo.clean_geometry(mesh)
    return mesh,

def clean_orphans (mesh) :
    """
    remove wisps with no regions
    """
    algo.clean_orphans(mesh)
    return mesh,
