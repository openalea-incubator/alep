""" Gather different strategies for modeling dispersal of fungus propagules.

"""

# Imports #########################################################################
import random
import numpy
import collections
from alinea.alep.fungus import Fungus
from alinea.alep.architecture import get_leaves

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
    
    def __init__(self, k=0.148, precision=0.01, label='LeafElement', fungus=None):
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
        if fungus is None:
            self.fungus = Fungus()
        else:
            self.fungus = fungus
   
    def disperse(self, g, dispersal_units, time_control = None, **kwds):
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
            from openalea.plantgl import all as pgl
            from collections import OrderedDict
            from math import exp, pi, cos, sin, tan
            from random import shuffle
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
                borders = numpy.linspace(0,1,num=nb_tri)
                
                dus_by_tri = {k: int(dus / nb_tri) for k in range(nb_tri)}
                dus_by_tri[nb_tri-1] += dus % nb_tri
                
                for k,n in dus_by_tri.iteritems():
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
                        sphere_area = 2*pi*distances.values()[-1]**2
                        for leaf_id in distances:
                            if n >= 1:
                                area_factor = areas[leaf_id]/sphere_area
                                distance_factor = exp(-self.k * distances[leaf_id])
                                qc = min(n, (n * area_factor * distance_factor))
                                deposits[leaf_id] = qc
                                n -= qc

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
                        n_big = int(n*(beta+pi/2.)/pi)
                        n_small = n - n_big
                        for leaf in projection:
                            if n >= 1:
                                copy_centroid = copy(centroids[leaf])
                                copy_origin = copy(Origin)
                                copy_centroid[2] = 0.
                                copy_origin[2] = 0.
                                if (centroids[leaf]-Origin)*(source_normal[0],source_normal[1],0) >= 0:
                                    area_factor = areas[leaf]/(pi*dmax**2/2.)
                                    dist = pgl.norm(copy_centroid-copy_origin)
                                    distance_factor = exp(-self.k * dist)
                                    qc = min(n_big, (n_big * area_factor * distance_factor))
                                else:
                                    area_factor = areas[leaf]/(pi*a*b/2.)
                                    dist = pgl.norm(copy_centroid-copy_origin)/abs(cos(pgl.angle(source_normal, (1,0,0))+pi/2.))
                                    distance_factor = exp(-self.k * dist)
                                    qc = min(n_small, (n_small * area_factor * distance_factor))
                                    # import pdb
                                    # pdb.set_trace()
                                qc = min(qc, n)
                                deposits[leaf] = qc
                                n-=qc
                    
        for vid, dep in deposits.iteritems():
            du = self.fungus.dispersal_unit()
            du.set_nb_dispersal_units(nb_dispersal_units=dep)
            deposits[vid] = [du]
        
        return deposits
        
    def plot_distri_3d(self, g):
        from alinea.alep.architecture import set_property_on_each_id
        from alinea.alep.alep_color import alep_colormap, green_yellow_red
        from alinea.alep.disease_outputs import plot3d_transparency
        from openalea.plantgl.all import Viewer
        import matplotlib.pyplot as plt
        # Compute severity by leaf
        dus = g.property("dispersal_units")
        deposits = {k:sum([du.nb_dispersal_units for du in v]) for k,v in dus.iteritems()}
        set_property_on_each_id(g, 'nb_dispersal_units', deposits)
    
        # Visualization
        vmax = 2*max(deposits.values())/3
        g = alep_colormap(g, 'nb_dispersal_units', cmap=green_yellow_red(), 
                          lognorm=False, zero_to_one=False, 
                          vmax = vmax)
        d = [numpy.arange(vmax)]
        fig, ax = plt.subplots()
        ax.imshow(d, cmap = green_yellow_red())

        for id in g:
            if not id in deposits:
                g.node(id).color = (0,0,0)
                g.node(id).transparency = 0.7
            elif deposits[id]==0.:
                g.node(id).color = (0,0,0)
                g.node(id).transparency = 0.7
            else:
                g.node(id).transparency = 0.                      
        scene = plot3d_transparency(g)
        Viewer.display(scene)
        
# Powdery mildew wind dispersal ###################################################
class PowderyMildewWindDispersal:
    """ Template class for a model of dispersal by wind that complies with the 
    guidelines of Alep.
    
    A class for a model of dispersal must include a method 'disperse'. In this example,
    dispersal units are distributed in a cone of dispersal in the direction of the wind.
    This model is adapted from the work of Calonnec et al., 2008 on powdery mildew.    
    """
    
    def __init__(self, cid=0.04, a0=45., reduction=100., k_beer=0.5, label='lf', wind_direction=(1, 0, 0), fungus=None):
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
        wind_direction (3-tuple) the wind direction vector
        """
        self.cid = cid
        self.a0 = a0
        self.reduction = reduction
        self.k_beer = k_beer
        self.label = label
        self.wind_direction = wind_direction
        print('')
        print('Be careful conversion l.228 dispersal')
        print('Be careful Beer law commented l.229 dispersal')
        
        if fungus is None:
            self.fungus = Fungus()
        else:
            self.fungus = fungus
        
    
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
            from openalea.plantgl import all as pgl
            from random import shuffle
            from math import degrees, exp, tan, pi, radians
            from collections import OrderedDict
            geometries = g.property('geometry')
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
                        if (centroids[vid]-Origin)*wind_directions.get(source, self.wind_direction) >= 0}
                
                # Sort the vids based on the angle                
                angles = {vid:degrees(pgl.angle(vect, wind_directions.get(source, self.wind_direction))) 
                          for vid, vect in vects.iteritems()
                          if degrees(pgl.angle(vect, wind_directions.get(source, self.wind_direction)))<= self.a0}
                
                # Sort the vids based on the distance
                distances = {vid:pgl.norm(vects[vid]) for vid in angles}
                distances = OrderedDict(sorted(distances.iteritems(), key=lambda x: x[1]))
                               
                # Beer law inside cone to take into account leaf coverage
                n = dus

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
                                            
                        deposits[leaf] = int(qc)
                        n -= int(qc)
                        # if len(dus) < 1 or len(deposits[leaf]) < 1:
                        if n < 1:
                            break
                            
        for vid, dep in deposits.iteritems():
            du = self.fungus.dispersal_unit()
            du.set_nb_dispersal_units(nb_dispersal_units=dep)
            deposits[vid] = [du]
                        
        return deposits
        
        
    def plot_distri_3d(self, g):
        from alinea.alep.architecture import set_property_on_each_id
        from alinea.alep.alep_color import alep_colormap, green_yellow_red
        from alinea.alep.disease_outputs import plot3d_transparency
        from openalea.plantgl.all import Viewer
        import matplotlib.pyplot as plt
        # Compute severity by leaf
        dus = g.property("dispersal_units")
        deposits = {k:sum([du.nb_dispersal_units for du in v]) for k,v in dus.iteritems()}
        set_property_on_each_id(g, 'nb_dispersal_units', deposits)
    
        # Visualization
        vmax = 2*max(deposits.values())/3
        g = alep_colormap(g, 'nb_dispersal_units', cmap=green_yellow_red(), 
                          lognorm=False, zero_to_one=False, 
                          vmax = vmax)
        d = [numpy.arange(vmax)]
        fig, ax = plt.subplots()
        ax.imshow(d, cmap = green_yellow_red())

        for id in g:
            if not id in deposits:
                g.node(id).color = (0,0,0)
                g.node(id).transparency = 0.7
            elif deposits[id]==0.:
                g.node(id).color = (0,0,0)
                g.node(id).transparency = 0.7
            else:
                g.node(id).transparency = 0.                      
        scene = plot3d_transparency(g)
        Viewer.display(scene)
        
# Brown Rust wind dispersal (horizontal layers) ###############################

def compute_overlaying(nb_du, area, impact_surface):
    # true_area_impacted = areas[vid] * (1 - min(1., numpy.exp(-nb_du * impact_surface / areas[vid])))
    # new_nb_du = int(true_area_impacted / diameter)
    return max(1, int((1 - min(1., numpy.exp(-nb_du * impact_surface / area)))*(area/impact_surface)))

class BrownRustDispersal:
    import sys
    sys.setrecursionlimit(10000)

    """ Calculate distribution of dispersal units in horizontal layers """
    def __init__(self, fungus = None,
                 group_dus = False,
                 domain_area = 1.,
                 convUnit = 0.01,
                 layer_thickness = 1.,
                 k_dispersal = 0.07,
                 k_beer = 0.5):             
        if fungus is not None:
            self.fungus = fungus
        else:
            self.fungus = Fungus()
        self.group_dus = group_dus
        self.domain_area = domain_area
        self.convUnit = convUnit
        self.layer_thickness = layer_thickness
        self.k_dispersal = k_dispersal
        self.k_beer = k_beer

    def leaves_in_grid(self, g, label = 'LeafElement'):
        from openalea.plantgl import all as pgl
        from alinea.alep.fungus import Fungus
        geometries = g.property('geometry')
        centroids = g.property('centroid')
        areas = g.property('area')
        tesselator = pgl.Tesselator()
        bbc = pgl.BBoxComputer(tesselator)      
        leaves = get_leaves(g, label=label)
        leaves = [l for l in leaves if l in geometries]
        
        # Get centroids        
        def centroid(vid):
            if is_iterable(geometries[vid]):
                bbc.process(pgl.Scene(geometries[vid]))
            else:
                bbc.process(pgl.Scene([pgl.Shape(geometries[vid])]))
            center = bbc.result.getCenter()
            centroids[vid] = center
        
        if len(leaves)>0:
            for vid in leaves:
                centroid(vid)
                
            # Define grid (horizontal layers)
            zs = [c[2] for c in centroids.itervalues()]
            minz = min(zs)
            maxz = max(zs) + self.layer_thickness
            layers = {l:[] for l in numpy.arange(minz, maxz, self.layer_thickness)}
            
            # Distribute leaves in layers
            for vid, coords in centroids.iteritems():
                z = coords[2]
                ls = layers.keys()
                i_layer = numpy.where(map(lambda x: x<=z<x+self.layer_thickness 
                                        if z!=maxz - self.layer_thickness
                                        else x<=z<=x+self.layer_thickness , ls))[0]
                if len(i_layer) > 0. and areas[vid]>1e-10:
                    layers[ls[i_layer]].append(vid)
            
            self.layers = layers
        else:
            self.layers = {}
    def get_dispersal_units(self, g, 
                            fungus_name = "brown_rust", 
                            label = 'LeafElement',
                            weather_data = None, **kwds):
        lesions = g.property('lesions')
        return {vid:sum([l.emission() for l in les if l.is_sporulating 
                        and l.fungus.name.startswith(fungus_name)]) 
                        for vid, les in lesions.iteritems()}

    def disperse(self, g, dispersal_units = {}, weather_data = None,
                 label='LeafElement', domain_area=None, **kwds):
                     
        def expo_decrease(dist):
            return numpy.exp(-self.k_dispersal*dist)
            
        def left_in_canopy(nb_dus, layer, max_height):
            xmax = numpy.arange(100)
            ymax = expo_decrease(xmax)
            x_up = numpy.arange(max_height-layer)
            y_up = expo_decrease(x_up)
            upward = numpy.trapz(y_up,x_up)/numpy.trapz(ymax,xmax)
            x_down = numpy.arange(max_height)
            y_down = expo_decrease(x_down)
            downward = numpy.trapz(y_down,x_down)/numpy.trapz(ymax,xmax)
            return nb_dus*(upward+downward)/2.
        
        self.leaves_in_grid(g, label=label)
        # Group DUs in layers and reduce global number according to Beer Law
        areas = g.property('area')
        geom = g.property('geometry')
        labels = g.property('label')
        areas = {k:v for k,v in areas.iteritems() 
                if k in geom and labels[k].startswith('LeafElement')}
        total_area = sum(areas.values())
        if domain_area is None:
            domain_area = self.domain_area
        lai = total_area*(self.convUnit**2)/domain_area if domain_area>0. else 0.
        beer_factor = 1-numpy.exp(-self.k_beer * lai)
        dus_by_layer = {layer:sum([dispersal_units[v]
                        for v in vids if v in dispersal_units])*beer_factor
                        for layer, vids in self.layers.iteritems()}
        max_height = max(self.layers.keys())
        dus_by_layer = {layer:left_in_canopy(nb_dus, layer, max_height)
                        for layer, nb_dus 
                        in dus_by_layer.iteritems() if nb_dus>0.}
                
        def sum_nb(nb_leaves, nb_du):
            if nb_leaves == 1:
                return [nb_du]
            elif nb_du == 0:
                return [0] + sum_nb(nb_leaves-1, nb_du)
            else:
                nb_du_avg = float(nb_du/nb_leaves)
                nb_du_sup = 2.*nb_du_avg
                if nb_du_sup >= 1:
                    nb_on_vid = int(round(max(0, min(nb_du, numpy.random.normal(nb_du_avg, nb_du_sup)))))
                else:
                    nb_on_vid = 1 if numpy.random.random()<nb_du_sup else 0
                return [nb_on_vid] + sum_nb(nb_leaves-1, nb_du - nb_on_vid)
                
        deposits = {}
        for source, nb_dus in dus_by_layer.iteritems():
            probas = {}
            for target, vids in self.layers.iteritems():
                nb_vids = len(vids)
                if nb_vids>0:
                    dist = abs(source - target)
                    area_l = (sum([areas[vid] for vid in vids])*self.convUnit**2)/self.domain_area
                    proba_lai = 1-numpy.exp(-self.k_beer * area_l)
                    proba_dist = numpy.exp(-self.k_dispersal*dist)
                    try:
                        probas[target] += proba_lai*proba_dist
                    except:
                        probas[target] = proba_lai*proba_dist
            total_probas = sum([p for p in probas.itervalues()])
            
            for layer, proba in probas.iteritems():
                proba /= total_probas
                nb_depo_layer = numpy.random.binomial(nb_dus, proba)
                vids = self.layers[layer]
                nb_vids = len(vids)
                distribution_by_leaf = sum_nb(nb_vids, nb_depo_layer)
                numpy.random.shuffle(distribution_by_leaf)
                for i_lf, lf in enumerate(vids):
                    if distribution_by_leaf[i_lf] > 0.:
                        depo = distribution_by_leaf[i_lf]
                        depo = compute_overlaying(depo, areas[lf], numpy.pi*0.0015**2)
                        try:
                            deposits[lf] += depo
                        except:
                            deposits[lf] = depo

        for vid, nb_dus in deposits.iteritems():
            if self.group_dus==True:
                du = self.fungus.dispersal_unit()
                du.set_nb_dispersal_units(nb_dispersal_units = nb_dus)
                deposits[vid] = [du]
            else:
                dus = []
                for d in range(nb_dus):
                    du = self.fungus.dispersal_unit(1)
                    dus.append(du)
                deposits[vid] = dus
        return deposits

    def plot_layers(self, g):
        from alinea.alep.architecture import set_property_on_each_id
        from alinea.alep.alep_color import alep_colormap
        from alinea.alep.disease_outputs import plot3d_transparency
        from openalea.plantgl.all import Viewer
        # Compute severity by leaf
        self.leaves_in_grid(g)
        layers = self.layers
        layer_by_leaf = {vid:k for k,v in layers.iteritems() for vid in v}
        set_property_on_each_id(g, 'height', layer_by_leaf)
    
        # Visualization
        g = alep_colormap(g, 'height', cmap='prism', 
                          lognorm=False, zero_to_one=False)

        geometries = g.property('geometry') 
        leaves = get_leaves(g, label='LeafElement')
        leaves = [l for l in leaves if l in geometries] 
        transp = {vid:0. for k,v in layers.iteritems() for vid in v}
        set_property_on_each_id(g, 'transparency', transp)                         
        scene = plot3d_transparency(g)
        Viewer.display(scene)
        
    def view_distri_layers(self, g, nb_dispersal_units = 1e5, vmax = None,
                           position_source = 3./5, return_df=False):
        from alinea.alep.architecture import set_property_on_each_id
        from alinea.alep.alep_color import alep_colormap, green_yellow_red
        from alinea.alep.disease_outputs import plot3d_transparency
        from openalea.plantgl.all import Viewer
        import matplotlib.pyplot as plt
        import pandas as pd
        # Compute severity by leaf
        self.leaves_in_grid(g)
        layers = self.layers.keys()
        layers.sort()        
        layer = layers[int(position_source*len(layers))]
        if len(self.layers[layer])>0:
            leaf = self.layers[layer][0]
        else:
            non_empty = {k:v for k,v in self.layers.iteritems() if len(v)>0}
            leaf = self.layers[min(non_empty.keys(), key=lambda k:abs(k-layer))][0]
        deposits = {k:sum([du.nb_dispersal_units for du in v]) for k,v in 
                    self.disperse(g, dispersal_units = {leaf : nb_dispersal_units}).iteritems()}  
        set_property_on_each_id(g, 'nb_dispersal_units', deposits)
    
        # Visualization
        if vmax is None:
            vmax = 2*max(deposits.values())/3
        g = alep_colormap(g, 'nb_dispersal_units', cmap=green_yellow_red(), 
                          lognorm=False, zero_to_one=False, 
                          vmax = vmax)
        d = [numpy.arange(vmax)]
        fig, ax = plt.subplots()
        ax.imshow(d, cmap = green_yellow_red())

        transp = {vid:0. for k,v in self.layers.iteritems() for vid in v}
        set_property_on_each_id(g, 'transparency', transp)   
        for id in g:
            if not id in deposits:
                g.node(id).color = (0,0,0)
                g.node(id).transparency = 0.7
            elif deposits[id]==0.:
                g.node(id).color = (0,0,0)
                g.node(id).transparency = 0.7
            else:
                g.node(id).transparency = 0.                      
        scene = plot3d_transparency(g)
        Viewer.display(scene)
        
        if return_df==True:        
            depo_layer = {k:sum([deposits[vid] for vid in v if vid in deposits])
                            for k,v in self.layers.iteritems()}
            df = pd.DataFrame([[k,v] for k, v in depo_layer.iteritems()])
            df = df.sort(0)
            df[2] = df[1]/df[1].sum()
            fig, ax = plt.subplots()
            ax.plot(df[2], df[0], 'k')
            ax.set_ylabel('Height', fontsize = 16)
            ax.set_xlabel('Proportion of deposits', fontsize = 16)
            return df
        
    def plot_distri_layers(self, g, nb_dispersal_units = 1000, 
                           position_source = 3./5):
        import pandas as pd
        import matplotlib.pyplot as plt
        # Compute severity by leaf
        self.leaves_in_grid(g)
        layers = self.layers.keys()
        layers.sort()        
        layer = layers[int(position_source*len(layers))]
        if len(self.layers[layer])>0:
            leaf = self.layers[layer][0]
        else:
            non_empty = {k:v for k,v in self.layers.iteritems() if len(v)>0}
            leaf = self.layers[min(non_empty.keys(), key=lambda k:abs(k-layer))][0]
        deposits = {k:sum([du.nb_dispersal_units for du in v]) for k,v in 
                    self.disperse(g, dispersal_units = {leaf : nb_dispersal_units}).iteritems()}
        depo_layer = {k:sum([deposits[vid] for vid in v if vid in deposits])
                        for k,v in self.layers.iteritems()}
        df = pd.DataFrame([[k,v] for k, v in depo_layer.iteritems()])
        df = df.sort(0)
        df[2] = df[1]/df[1].sum()
        fig, ax = plt.subplots()
        ax.plot(df[2], df[0], 'k')
        ax.set_ylabel('Height', fontsize = 16)
        ax.set_xlabel('Proportion of deposits', fontsize = 16)
        
        return d