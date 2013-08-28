""" Gather different strategies for modeling dispersal of fungus propagules.

"""

# Imports #########################################################################
import random

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

# Powdery mildew wind dispersal ###################################################
class PowderyMildewWindDispersal:
    """ Template class for a model of dispersal by wind that complies with the 
    guidelines of Alep.
    
    A class for a model of dispersal must include a method 'disperse'. In this example,
    dispersal units are distributed in a cone of dispersal in the direction of the wind.
    This model is adapted from the work of Calonnec et al., 2008 on powdery mildew.    
    """
    
    def __init__(self, cid=0.04, a0=45., reduction=100., label='lf'):
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
        label: str
            Label of the part of the MTG concerned by the calculation
        """
        self.cid = cid
        self.a0 = a0
        self.reduction = reduction
        self.label = label
    
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
            from math import degrees, exp
            from collections import OrderedDict
            geometries = g.property('geometry')
            centroids = g.property('centroid')
            areas = g.property('area')
            wind_directions = g.property('wind_direction')
            tesselator = pgl.Tesselator()
            bbc = pgl.BBoxComputer(tesselator)
        
            leaves = get_leaves(g, label=self.label)

            def centroid(vid):
                bbc.process(pgl.Scene(geometries[vid]))
                center = bbc.result.getCenter()
                centroids[vid] = center
            
            def area(vid):
                areas[vid] = pgl.surface(geometries[vid][0])*1000
            
            for source, dus in dispersal_units.iteritems():
                # TODO: Special computation for interception by source leaf
                
                # All other leaves are potential targets
                targets = list(leaves)
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
                
                # Distribute the dispersal units                
                shuffle(dus)
                n = len(dus)
                for leafid in distances:
                    qc = min(n, (n * (areas[leafid]/self.reduction) * 
                         exp(-self.cid * distances[leafid]) * 
                         (self.a0 - angles[leafid])/self.a0))
                    if qc < 1:
                        for d in dus:
                            d.disable()
                        break
                                        
                    deposits[leafid] = dus[:int(qc)]
                    del dus[:int(qc)]
                    if len(dus) < 1 or len(deposits[leafid]) < 1:
                        for d in dus:
                            d.disable()
                        break
                        
        return deposits