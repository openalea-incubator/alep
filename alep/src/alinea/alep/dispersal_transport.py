""" Gather different strategies for modeling dispersal of fungus propagules.

"""

# Imports #########################################################################
import random
import collections

# Useful function #################################################################
def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)

# Random dispersal ################################################################
class RandomDispersal:
    """ Template class for a dispersal model that complies with the guidelines of Alep.
    
    A class for a model of dispersal must include a method 'disperse'. In this example,
    dispersal units are randomly distributed.
    
    """
    def disperse(self, g, dispersal_units, time_control = None):
        """ Example method for dispersal with random distribution.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        dispersal_units : dict
            Dispersal units emitted by the lesions on leaves
            
        Returns
        -------
        deposits : dict
            Dispersal units deposited on new position on leaves
        """
        # vids = scene.todict().keys()
        try:
            dt = time_control.dt
        except:
            dt = 1
        
        vids = [id for id,v in g.property('geometry').iteritems()]
        n = len(vids)
        deposits = {}
        if dt > 0:
            for vid, dlist in dispersal_units.iteritems():
                for d in dlist:
                    if random.random() < 0.1:
                        if n>=1:
                            idx = random.randint(0,n-1)
                            v = vids[idx]
                            d.set_position([0, 0])
                            deposits.setdefault(v,[]).append(d)
        
        return deposits
        
# Septoria rain dispersal #########################################################
class SeptoriaRainDispersal:
    """ Template class for a model of dispersal by rain that complies with the 
    guidelines of Alep.
    
    A class for a model of dispersal must include a method 'disperse'. In this example,
    dispersal units are first distributed updward in a semi-sphere normal to source leaf.
    After, those which are left are deposited downward on leaves comprised in a cylinder
    whose dimension is calculated according to the size of the semi-sphere.    
    """
    
    def __init__(self, k=0.148, precision=0.01, label='LeafElement'):
        """ Initialize the model with fixed parameters.
        
        Parameters
        ----------
        k: float
            Shape parameter of the exponential decreasing function
        precision: float
            Value defining distance max with the exponential decreasing function 
        label: str
            Label of the part of the MTG concerned by the calculation
        """
        from math import log
        self.label = label
        self.k = k
        self.distance_max = -log(precision)/k
   
    def disperse(self, g, dispersal_units, time_control = None):
        """ Compute distribution of dispersal units by rain splash.
        
        1. Upward dispersal:
        For each source of dispersal units, create a semi-sphere of dispersal 
        normal to the surface of source leaf. In the semi-sphere, target 
        leaves are sorted according to the distance from the source.                
        Then distribute dispersal units from the closer target to the more far.
        The number of dispersal units by target leaf is computed as in Robert et al. 2008
        
        2. Downward dispersal:
        Get leaves in a cylinder whose dimensions are related to the dimesions
        of the semi-sphere in step 1. Then distribute dispersal units from top to bottom.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        dispersal_units : dict
            Dispersal units emitted by the lesions on leaves
            
        Returns
        -------
        deposits : dict
            Dispersal units deposited on new position on leaves
        """
        try:
            dt = time_control.dt
        except:
            dt = 1
        
        deposits = {}
        if dt>0:
            from alinea.astk.plantgl_utils import get_area_and_normal
            from alinea.alep.architecture import get_leaves
            from openalea.plantgl import all as pgl
            from collections import OrderedDict
            from math import exp, pi, cos, sin, tan
            from random import shuffle
            import numpy as np
            from copy import copy
            
            dmax = self.distance_max
            tesselator = pgl.Tesselator()
            bbc = pgl.BBoxComputer(tesselator)
            leaves = get_leaves(g, label=self.label)
            centroids = g.property('centroid')
            geometries = g.property('geometry')
            _, norm = get_area_and_normal(geometries)
            areas = g.property('area')          
            
            def centroid(vid):
                if is_iterable(geometries[vid]):
                    bbc.process(pgl.Scene(geometries[vid]))
                else:
                    bbc.process(pgl.Scene([pgl.Shape(geometries[vid])]))
                center = bbc.result.getCenter()
                centroids[vid] = center
            
            for source, dus in dispersal_units.iteritems():
                nb_tri = len(norm[source])
                borders = np.linspace(0,1,num=nb_tri)
                
                dus_by_tri = {k:filter(lambda x: borders[k]<x.position[0]<=borders[k+1], dus) 
                                for k in range(nb_tri-1)
                                if len(filter(lambda x: borders[k]<x.position[0]<=borders[k+1], dus))>0.}
                
                for k,v in dus_by_tri.iteritems():
                    source_normal = norm[source][k]
                    
                    ## UPWARD ##
                    # All leaves except the source are potential targets
                    targets = list(leaf for leaf in leaves if leaf in geometries.iterkeys())
                    targets.remove(source)
                    
                    # Compute centroids
                    centroid(source)
                    for vid in targets:
                        centroid(vid)
                    
                    # Sort the vids based on the direction 
                    # TODO : modify source angle
                    Origin = centroids[source]
                    vects = {vid:(centroids[vid]-Origin) for vid in targets 
                            if (centroids[vid]-Origin)*source_normal >= 0}
                    
                    # Sort the vids based on the distance
                    distances = {vid:pgl.norm(vects[vid]) for vid in vects if pgl.norm(vects[vid])<dmax}
                    distances = OrderedDict(sorted(distances.iteritems(), key=lambda x: x[1]))
                    
                    # Distribute the dispersal units
                    if len(distances.values())>0:
                        shuffle(v)
                        n = len(v)
                        sphere_area = 2*pi*distances.values()[-1]**2
                        for leaf_id in distances:
                            area_factor = areas[leaf_id]/sphere_area
                            distance_factor = exp(-self.k * distances[leaf_id])
                            qc = min(n, (n * area_factor * distance_factor))
                            
                            deposits[leaf_id] = v[:int(qc)]
                            del v[:int(qc)]
                            # if len(dus) < 1 or len(deposits[leafid]) < 1:
                            if len(v) < 1:
                                for d in v:
                                    d.disable()
                                # break
                    
                    ## DOWNWARD ##
                    vects2 = {vid:(centroids[vid]-Origin) for vid in targets if not vid in vects}
                    projection = {}

                    alpha = pgl.angle(source_normal, (1,0,0))
                    if alpha>=pi/2. or (alpha<pi/2. and source_normal[2]>=0):
                        alpha+=pi/2.
                    beta = pgl.angle(source_normal, (0,0,1))
                    a = dmax
                    b = dmax*cos(beta)
                    
                    for leaf in vects2:
                        if (centroids[leaf]-Origin)*(source_normal[0], source_normal[1], 0) >= 0:
                            # Big side of the projection semi circle
                            copy_centroid = copy(centroids[leaf])
                            copy_origin = copy(Origin)
                            copy_centroid[2] = 0.
                            copy_origin[2] = 0.
                            if pgl.norm(copy_centroid-copy_origin) < dmax:
                                projection[leaf] = vects2[leaf]
                        else:
                            # Small side of the projection semi ellipse
                            x = vects2[leaf][0]
                            y = vects2[leaf][1]
                            x2 = x*cos(alpha)+y*sin(alpha)
                            y2 = -x*sin(alpha)+y*cos(alpha)
                            if (x2**2)/(a**2) + (y2**2)/(b**2) < 1 :
                                projection[leaf] = vects2[leaf]
                    projection = OrderedDict(sorted(projection.items(), key=lambda x:x[1][2], reverse=True))
                    
                    if len(projection)>0:
                        shuffle(v)
                        n = len(v)
                        n_big = int(n*(beta+pi/2.)/pi)
                        n_small = n - n_big
                        for leaf in projection:
                            copy_centroid = copy(centroids[leaf])
                            copy_origin = copy(Origin)
                            copy_centroid[2] = 0.
                            copy_origin[2] = 0.
                            if (centroids[leaf]-Origin)*(source_normal[0],source_normal[1],0) >= 0:
                                area_factor = areas[leaf]/(pi*dmax**2/2.)
                                dist = pgl.norm(copy_centroid-copy_origin)
                                distance_factor = exp(-self.k * dist)
                                qc = min(n_big, (n_big * area_factor * distance_factor))
                                g.node(leaf).color = (0, 180, 0)
                            else:
                                area_factor = areas[leaf]/(pi*a*b/2.)
                                dist = pgl.norm(copy_centroid-copy_origin)/abs(cos(pgl.angle(source_normal, (1,0,0))+pi/2.))
                                distance_factor = exp(-self.k * dist)
                                qc = min(n_small, (n_small * area_factor * distance_factor))
                                g.node(leaf).color = (0, 0, 180)
                                # import pdb
                                # pdb.set_trace()
                            deposits[leaf] = v[:int(qc)]
                            del v[:int(qc)]
                            
                    for leaf in distances:
                        g.node(leaf).color = (180, 0, 0)
                    
                    # Temp
                    # from alinea.adel.mtg_interpreter import plot3d
                    # from openalea.plantgl.all import Viewer
                    # g.node(source).color=(230, 62, 218)
                    # scene = plot3d(g)
                    # Viewer.display(scene)
                    # import pdb
                    # pdb.set_trace()
        return deposits
        
# Powdery mildew wind dispersal ###################################################
class PowderyMildewWindDispersal:
    """ Template class for a model of dispersal by wind that complies with the 
    guidelines of Alep.
    
    A class for a model of dispersal must include a method 'disperse'. In this example,
    dispersal units are distributed in a cone of dispersal in the direction of the wind.
    This model is adapted from the work of Calonnec et al., 2008 on powdery mildew.    
    """
    
    def __init__(self, cid=0.04, a0=45., reduction=100., k_beer=0.5, label='lf'):
        """ Initialize the model with fixed parameters.
        
        Parameters
        ----------
        cid: float
            Spore decay with distance
        a0: float
            Angle of the cone of dispersal
        reduction: float
            Reduction parameter to limit the number of spores reaching a surface.
            Regarding the equation, this parameter must be homogeneous to a surface.
        k_beer: float
            Value of k in the Beer-Lambert law for simplification of interception 
            of spores in the canopy.
        label: str
            Label of the part of the MTG concerned by the calculation
        """
        self.cid = cid
        self.a0 = a0
        self.reduction = reduction
        self.k_beer = k_beer
        self.label = label
        print('')
        print('Be careful conversion l.228 dispersal')
        print('Be careful Beer law commented l.229 dispersal')
    
    def disperse(self, g, dispersal_units, time_control = None):
        """ Compute dispersal of spores of powdery mildew by wind in a cone.
        
        For each source of dispersal units, create a cone of dispersal 
        in which target leaves are sorted:
        1. according to the wind direction
        2. according to the angle a0
        3. according to the distance from the source
        
        Then distribute dispersal units from the closer target to the more far.
        The number of dispersal units by target leaf is computed as in
        Calonnec et al. 2008.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        dispersal_units : dict
            Dispersal units emitted by the lesions on leaves
            
        Returns
        -------
        deposits : dict
            Dispersal units deposited on new position on leaves
        """
        try:
            dt = time_control.dt
        except:
            dt = 1
        
        deposits = {}
        if dt > 0:
            from alinea.alep.architecture import get_leaves
            from openalea.plantgl import all as pgl
            from random import seed, choice, shuffle
            from math import degrees, exp, tan, pi, radians
            from collections import OrderedDict
            geometries = g.property('geometry')
            lengths = g.property('length')
            centroids = g.property('centroid')
            areas = g.property('area')
            wind_directions = g.property('wind_direction')
            tesselator = pgl.Tesselator()
            bbc = pgl.BBoxComputer(tesselator)
        
            leaves = get_leaves(g, label=self.label)

            def centroid(vid):
                if is_iterable(geometries[vid]):
                    bbc.process(pgl.Scene(geometries[vid]))
                else:
                    bbc.process(pgl.Scene([pgl.Shape(geometries[vid])]))
                center = bbc.result.getCenter()
                centroids[vid] = center
            
            def area(vid):
                # areas[vid] = pgl.surface(geometries[vid][0])*1000
                areas[vid] = pgl.surface(geometries[vid][0])

            for source, dus in dispersal_units.iteritems():
                # TODO: Special computation for interception by source leaf
                
                # All other leaves are potential targets
                targets = list(leaf for leaf in leaves if leaf in geometries.iterkeys())
                targets.remove(source)
                
                # Compute centroids
                centroid(source)
                for vid in targets:
                    centroid(vid)
                    # surface(vid)

                # Sort the vids based on the direction 
                Origin = centroids[source]
                vects = {vid:(centroids[vid]-Origin) for vid in targets 
                        if (centroids[vid]-Origin)*wind_directions[source] >= 0}
                
                # Sort the vids based on the angle                
                angles = {vid:degrees(pgl.angle(vect, wind_directions[source])) 
                          for vid, vect in vects.iteritems()
                          if degrees(pgl.angle(vect, wind_directions[source]))<= self.a0}
                
                # Sort the vids based on the distance
                distances = {vid:pgl.norm(vects[vid]) for vid in angles}
                distances = OrderedDict(sorted(distances.iteritems(), key=lambda x: x[1]))
                               
                # Beer law inside cone to take into account leaf coverage
                shuffle(dus)
                n = len(dus)

                if len(distances.values())>0:
                    for leaf in distances:
                        # qc = min(n, (n * (areas[leaf]/self.reduction) * 
                             # exp(-self.cid * distances[leaf]) * 
                             # (self.a0 - angles[leaf])/self.a0))
                        surf_base_cone = pi*(tan(radians(self.a0))*distances[leaf])**2
                        area_factor = min(1, areas[leaf]/surf_base_cone)
                        # import pdb
                        # pdb.set_trace()
                        qc = min(n, (n * area_factor * 
                             exp(-self.cid * distances[leaf]) * 
                             (self.a0 - angles[leaf])/self.a0))
                        
                        # if qc < 1:
                            # for d in dus:
                                # d.disable()
                            # break
                                            
                        deposits[leaf] = dus[:int(qc)]
                        del dus[:int(qc)]
                        # if len(dus) < 1 or len(deposits[leaf]) < 1:
                        if len(dus) < 1:
                            for d in dus:
                                d.disable()
                            break
                        
        return deposits
