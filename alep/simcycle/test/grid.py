""" Embed a set of triangles into a 1D or 3D grid.

"""
from openalea.plantgl.all import *
import numpy as np

def surface(pts):
    a, b, c = [Vector3(*x) for x in pts]
    return norm((b-a)^(c-a))/2.
    
def triangles():
    """ Generate a set of triangles for test purpose.
    """
    return [((0.5,0,0), (0,0,1), (1,0,1))]

def grid1d(h=1., n=2):
    return np.linspace(0,h,n+1,endpoint=False)

def fill():
    tris = triangles()
    z = grid1d()
    surfaces = np.zeros_like(z)
    for (p1, p2, p3) in tris:
        # check on which each points are
        pts = [p1, p2, p3]
        zs = p1[-1], p2[-1], p3[-1]
        # Compute where each z are in the grid (which z layer)
        indices = z.searchsorted(zs)-1

        # Reorder the three points to have always the same configuration
        # i.e. pt with z min, pt with z max, and middle point.
        m, M = indices.min(), indices.max()
        indices = list(indices)
        im = indices.index(m); del indices[im]
        iM = indices.index(M); del indices[iM]
        ii = indices[0]
        pm = pts[im]
        pM = pts[iM]
        pi = pts[ii]

        # discretize all the segments of the triangles
        # cases
        # 1. all points in the same layer
        if im == ii == iM:
            surfaces[im] += surface(pts)
        else:
            #S12: 
            dz= pi[-1]-pm[-1]
            Smi = [np.poly1d([(pm[i]+pi[i])/dz, -(pm[i]*pi[-1]+pm[-1]*pi[i] )/dz]) for i in range(3)]

            dz= pM[-1]-pi[-1]
            SiM = [np.poly1d([(pi[i]+pM[i])/dz, -(pi[i]*pM[-1]+pi[-1]*pM[i] )/dz]) for i in range(3)]

            dz= pM[-1]-pm[-1]
            SmM = [np.poly1d([(pm[i]+pM[i])/dz, -(pm[i]*pM[-1]+pm[-1]*pM[i] )/dz]) for i in range(3)]

            zs = z[im+1:ii]
            ptsmi = []
            if im == ii:
                # cut segments pmpM and pipM
                # add quad then tri
                pass

# p1 + t*(p2-p1)
# but test with t = z
# or compute t where Z1(z)+ts*(Z2-Z1) = Zs



        
