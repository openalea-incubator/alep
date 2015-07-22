""" Utilities for computing outputs from the disease model.

The aim of this module is to provide all the tools needed to compute
the outputs of the disease models. 
"""


def count_lesions(g):
    """ Count lesions of the mtg.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    nb_lesions: int
        Number of lesions on the MTG
    """
    return len(sum(g.property('lesions').values(), []))

def count_dispersal_units(g):
    """ Count dispersal units of the mtg.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    nb_dispersal_units: int
        Number of dispersal units on the MTG
    """
    return len(sum(g.property('dispersal_units').values(), []))
    
def count_lesions_by_leaf(g):
    """ Count lesions on each leaf of the MTG.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
        
    Returns
    -------
    nb_lesions_by_leaf: dict([id:nb_lesions])
        Number of lesions on each part of the MTG given by the label
    """
    lesions = g.property('lesions')
    return {k:len(v) for k,v in lesions.iteritems()}

def count_dispersal_units_by_leaf(g, label='LeafElement'):
    """ Count dispersal units on each part of the MTG given by the label.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    nb_dispersal_units_by_leaf: dict([id:nb_dispersal_units])
        Number of dispersal units on each part of the MTG given by the label
    """
    dispersal_units = g.property('dispersal_units')
    return {k:len(v) for k,v in dispersal_units.iteritems()}

def plot_lesions(g):
    """ plot the plant with infected elements in red """
    from alinea.adel.mtg_interpreter import plot3d
    from openalea.plantgl.all import Viewer
    
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'lesions' in n.properties():
            n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)

def plot_dispersal_units(g):
    """ plot the plant with infected elements in red """
    from alinea.adel.mtg_interpreter import plot3d
    from openalea.plantgl.all import Viewer
    
    green = (0,180,0)
    red = (180, 0, 0)
    for v in g.vertices(scale=g.max_scale()) : 
        n = g.node(v)
        if 'dispersal_units' in n.properties():
            n.color = red
        else : 
            n.color = green
    
    scene = plot3d(g)
    Viewer.display(scene)
    
def compute_lesion_areas_by_leaf(g, label='LeafElement'):
    """ Compute lesion area on each part of the MTG given by the label.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    lesion_surfaces_by_leaf: dict([id:surface_lesions])
        Surface of the lesions on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    lesions = g.property('lesions')
    return {vid:(sum(l.surface for l in lesions[vid])
            if vid in lesions.keys() else 0.) for vid in vids} 

def compute_green_lesion_areas_by_leaf(g, label='LeafElement'):
    """ Compute lesion areas on each green part of the MTG given by the label.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    green_lesion_area_by_leaf: dict([id:lesion_area])
        Surface of the lesions on each green part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    lesions = g.property('lesions')
    areas = g.property('area')
    green_lengths = g.property('green_length')
    sen_lengths = g.property('senesced_length')
    
    gla = {}
    for vid in vids:
        if vid in lesions.keys():
            les_surf = sum(l.surface_alive for l in lesions[vid])
            ratio_sen = sen_lengths[vid]/(sen_lengths[vid]+green_lengths[vid]) if (sen_lengths[vid]+green_lengths[vid])>0. else 0.
            # /!\ TODO : Can be replaced by green_areas[vid]/senesced_areas[vid]
            if les_surf<=areas[vid]:
                gla[vid]=les_surf*(1-ratio_sen)
            else:
                gla[vid]=les_surf-(areas[vid]*ratio_sen)
        else:
            gla[vid]=0.
    return gla

def compute_healthy_area_by_leaf(g, label='LeafElement'):
    """ Compute healthy area on each part of the MTG given by the label.
    
    Healthy area is green area (without senescence) minus the surface of lesions.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    healthy_by_leaf: dict([id:healthy_area])
        Healthy area on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
    # green_areas = g.property('green_area')

    areas = g.property('area')
    labels = g.property('label')
    # positions_senescence = g.property('position_senescence')
    sen_lengths = g.property('senesced_length')
    green_lengths = g.property('green_length')
    senesced_areas = {k:v*(sen_lengths[k]/(sen_lengths[k]+green_lengths[k]) if (sen_lengths[k]+green_lengths[k])>0. else 0.) for k,v in areas.iteritems() if labels[k].startswith(label)}
    
    # if len(positions_senescence)>0:
        # senesced_areas = {k:v*(1-positions_senescence[k]) for k,v in areas.iteritems() if labels[k].startswith(label)}
    # else:
        # senesced_areas = {k:0. for k,v in areas.iteritems() if labels[k].startswith(label)}
    green_lesion_areas = compute_green_lesion_areas_by_leaf(g, label)
    
    # return {vid:(areas[vid] - (senesced_areas[vid] + green_lesion_areas[vid])
        # if round(areas[vid],10)>round((senesced_areas[vid] + green_lesion_areas[vid]),10) else 0.)
        # for vid in vids}
        
    return {vid:(areas[vid] - (senesced_areas[vid] + green_lesion_areas[vid])) for vid in vids}
    
def compute_severity_by_leaf(g, label='LeafElement'):
    """ Compute severity of the disease on each part of the MTG given by the label.
    
    Severity is the ratio between disease surface and total leaf area (in %).
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    severity_by_leaf: dict([id:severity])
        Severity on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    leaves = get_leaves(g, label=label)
    total_areas = g.property('area')
    lesion_areas = compute_lesion_areas_by_leaf(g, label)
    
    # Calculate by blade
    blades = np.array_split(leaves,np.where(np.diff(leaves)!=1)[0]+1)
    sev={}
    for bl in blades:
        area_bl = np.array([total_areas[lf] for lf in bl])
        if any(area_bl==0.):
            for lf in bl[area_bl==0.]:
                sev[lf]=0.
            bl = np.delete(bl,np.where(area_bl==0.))
            area_bl = np.delete(area_bl,np.where(area_bl==0.))
        les_bl = np.array([lesion_areas[lf] for lf in bl])
        diff = area_bl - les_bl
        if any(diff<0):
            for lf in bl[diff<0]:
                sev[lf]=100.
            to_share = abs(sum(diff[diff<0]))
            bl = np.delete(bl,np.where(diff<0))
            area_bl = np.delete(area_bl,np.where(diff<0))
            diff = np.delete(diff,np.where(diff<0))
            diff*=1-to_share/sum(diff)
        for ind in range(len(bl)):
            if diff[ind]>area_bl[ind]:
                import pdb
                pdb.set_trace()
            sev[bl[ind]] = max(0, min(100, 100.*(1-diff[ind]/area_bl[ind])))
    
    #return {vid:(100*lesion_areas[vid]/float(total_areas[vid]) if total_areas[vid]>0. else 0.) for vid in vids}
    return sev
    
def compute_senescence_by_leaf(g, label='LeafElement'):
    """ Compute senescence on parts of the MTG given by the label.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    senescence_by_leaf: dict([id:senescence_area])
        Senescence on each part of the MTG given by the label    
    """
    labels = g.property('label')
    total_areas = {k:v for k,v in g.property('area').iteritems() if labels[k].startswith(label)}
    pos_sen = g.property('position_senescence')
    sen = {}
    for vid in total_areas.iterkeys():
        sen[vid] = total_areas[vid]*(1-pos_sen[vid])
    return sen
    
def compute_senescence_necrosis_by_leaf(g, label='LeafElement'):
    """ Compute senescence and lesion necrosis on green parts. 
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrosis_senescence_by_leaf: dict([id:necrosis_senescence_area])
        Senescence and lesion necrosis on each part of the MTG given by the label    
    """
    total_areas = g.property('area')
    healthy_areas = g.property('healthy_area')
    lesions = g.property('lesions')
    pos_sen = g.property('position_senescence')
    nec_sen = {}
    for vid in total_areas.iterkeys():
        if vid in lesions.keys():
            # Note G.Garin 16/12/13:
            # Little hack when senescence reaches leaf basis to account 
            # non localized lesion growth with available space. 
            nec_on_green = sum(lesion.necrotic_area for lesion in lesions[vid] if not lesion.is_senescent)
            les_on_green = sum(lesion.surface for lesion in lesions[vid] if not lesion.is_senescent)
            ratio_nec_on_green = nec_on_green/les_on_green if les_on_green>0. else 0.
            nec = min(nec_on_green, (total_areas[vid]*pos_sen[vid] - healthy_areas[vid])*ratio_nec_on_green)
            sen = total_areas[vid]*(1-pos_sen[vid])
            nec_sen[vid] = nec + sen
        else:
            nec_sen[vid] = 0.
    return nec_sen

def compute_necrosis_percentage_by_leaf(g, label='LeafElement'):
    """ Compute necrosis percentage on each part of the MTG given by the label.
    
    Necrosis percentage is the ratio between necrotic area and total leaf area.
    A tissue is necrotic if it is covered by a lesion in one of these states:
        - NECROTIC
        - SPORULATING
        - EMPTY
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrosis_by_leaf: dict([id:necrosis_percentage])
        Necrosis percentage on each part of the MTG given by the label
    """
    from alinea.alep.architecture import get_leaves
    leaves = get_leaves(g, label=label)
    total_areas = g.property('area')
    lesions = g.property('lesions')
    necrotic_areas = {}
    for vid in total_areas.iterkeys():
        if vid in lesions.keys():
            necrotic_areas[vid] = sum([lesion.necrotic_area for lesion in lesions[vid]])
        else:
            necrotic_areas[vid] = 0.
            
    # Calculate by blade
    blades = np.array_split(leaves,np.where(np.diff(leaves)!=1)[0]+1)
    necrosis_by_leaf={}
    for bl in blades:
        area_bl = np.array([total_areas[lf] for lf in bl])
        if any(area_bl==0.):
            for lf in bl[area_bl==0.]:
                necrosis_by_leaf[lf]=0.
            bl = np.delete(bl,np.where(area_bl==0.))
            area_bl = np.delete(area_bl,np.where(area_bl==0.))
        nec_bl = np.array([necrotic_areas[lf] for lf in bl])
        diff = area_bl - nec_bl
        if any(diff<0):
            for lf in bl[diff<0]:
                necrosis_by_leaf[lf]=100.
            to_share = abs(sum(diff[diff<0]))
            bl = np.delete(bl,np.where(diff<0))
            area_bl = np.delete(area_bl,np.where(diff<0))
            diff = np.delete(diff,np.where(diff<0))
            diff*=1-to_share/sum(diff)
        for ind in range(len(bl)):
            if diff[ind]>area_bl[ind]:
                import pdb
                pdb.set_trace()
            necrosis_by_leaf[bl[ind]] = 100.*(1-diff[ind]/area_bl[ind])
    
    #return {vid:(100*necrotic_areas[vid]/float(total_areas[vid]) if total_areas[vid]>0. else 0.) for vid in vids}
    return necrosis_by_leaf
    
def compute_necrotic_area_by_leaf(g, label='LeafElement'):
    """ Compute necrosis percentage on each part of the MTG given by the label.
    
    Necrosis percentage is the ratio between necrotic area and total leaf area.
    A tissue is necrotic if it is covered by a lesion in one of these states:
        - NECROTIC
        - SPORULATING
        - EMPTY
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrotic_area_by_leaf: dict([id:necrotic_area])
        Necrotic area on each part of the MTG given by the label
    """
    total_areas = g.property('area')
    lesions = g.property('lesions')
    necrotic_areas = {}
    for vid in total_areas.iterkeys():
        if vid in lesions.keys():
            necrotic_areas[vid] = sum(lesion.necrotic_area for lesion in lesions[vid])
        else:
            necrotic_areas[vid] = 0.
    return necrotic_areas
    
def compute_total_severity(g, label='LeafElement'):
    """ Compute disease severity on the whole plant.
    
    Severity is the ratio between disease surface and green leaf area (in %).
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    severity: float
        Ratio between disease surface and green leaf area (in %)
    """
    from numpy import mean
    severities = compute_severity_by_leaf(g, label=label)
    return mean(severities.values())
    
def compute_total_necrosis_percentage(g, label='LeafElement'):
    """ Compute necrosis percentage on the whole plant.
    
    Necrosis percentage ratio between necrotic (and sporulating) disease surface and total area of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrosis_percentage: float
        Ratio between necrotic (and sporulating) disease area and total area of leaves (in %)
    """   
    from numpy import mean
    nec = compute_necrosis_percentage_by_leaf(g, label=label)
    return mean(nec.values())

def compute_total_necrotic_area(g, label='LeafElement'):
    """ Compute necrosis percentage on the whole plant.
    
    Necrosis percentage ratio between necrotic (and sporulating) disease surface and total area of leaves.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    necrotic_area: float
        Total area of leaves covered by necrotic surfaces of lesions (in cm2)
    """
    from numpy import mean
    nec = compute_necrotic_area_by_leaf(g, label=label)
    return sum(nec.values())

def compute_normalised_audpc(necrosis, total_area):
    """ Compute the normalised AUDPC as in Robert et al. 2008
    
    "AUDPC is calculated as the area below the curve of pycnidia bearing
    necrotic leaf area. The latter is the leaf area that is or has been
    covered by sporulating lesions and thus represents the sporulation 
    potential of the leaf. The normalised AUDPC is obtained by dividing
    the AUDPC by a theoretical maximum value corresponding to the situation
    where the leaf is infected fully just after its emergence."
    
    Parameters
    ----------
    necrosis: list or array
        Historical values of necrosis percentage
    total_area: list or array
        Historical values of leaf area (same length than necrosis)
        
    Returns
    -------
    normalised_audpc: float
       AUDPC divided by a theoretical maximum value
    """
    import numpy as np
    from scipy.integrate import trapz
    full_necrosis = np.array([100. if total_area[k]>0. else 0. 
                              for k in range(len(total_area))])
    
    audpc = trapz(necrosis, dx=1)
    theo_audpc = trapz(full_necrosis, dx=1)
    return 100*audpc/theo_audpc if theo_audpc>0. else 0.
 
def plot3d_transparency(g, 
               leaf_material = None,
               stem_material = None,
               soil_material = None,
               colors = None,
               transparencies = None):
    """
    Returns a plantgl scene from an mtg.
    """
    from openalea.plantgl import all as pgl
    Material = pgl.Material
    Color3 = pgl.Color3
    Shape = pgl.Shape
    Scene = pgl.Scene
    
    if colors is None:
        if leaf_material is None:
            leaf_material = Material(Color3(0,180,0))
        if stem_material is None:
            stem_material = Material(Color3(0,130,0))
        if soil_material is None:
            soil_material = Material(Color3(170, 85,0))
        colors = g.property('color')

    transparencies = g.property('transparency')
    
    geometries = g.property('geometry')
    greeness = g.property('is_green')
    labels = g.property('label')
    scene = Scene()

    def geom2shape(vid, mesh, scene):
        shape = None
        if isinstance(mesh, list):
            for m in mesh:
                geom2shape(vid, m, scene)
            return
        if mesh is None:
            return
        if isinstance(mesh, Shape):
            shape = mesh
            mesh = mesh.geometry
        label = labels.get(vid)
        is_green = greeness.get(vid)
        if colors:
            if transparencies==None:
                shape = Shape(mesh, Material(Color3(* colors.get(vid, [0,0,0]) )))
            else:
                shape = Shape(mesh, Material(Color3(* colors.get(vid, [0,0,0]) ), transparency=transparencies.get(vid,0)))
        elif not greeness:
            if not shape:
                shape = Shape(mesh)
        elif label.startswith('Stem') and is_green:
            shape = Shape(mesh, stem_material)
        elif label.startswith('Leaf') and is_green:
            shape = Shape(mesh, leaf_material)
        elif not is_green:
            shape = Shape(mesh, soil_material)
        shape.id = vid
        scene.add(shape)

    for vid, mesh in geometries.iteritems():
        geom2shape(vid, mesh, scene)
    return scene

def plot_severity_by_leaf(g, senescence=True, transparency=None, label='LeafElement'):
    """ Display the MTG with colored leaves according to disease severity 
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy
    senescence: bool
        True if senescence must be displayed, False otherwise
    transparency: float[0:1]
        Transparency of the part of the MTG without lesion
    label: str
        Label of the part of the MTG concerned by the calculation
        
    Returns
    -------
    scene:
        Scene containing the MTG attacked by the disease
    
    """
    from alinea.alep.architecture import set_property_on_each_id, get_leaves
    from alinea.alep.disease_outputs import compute_severity_by_leaf
    from alinea.alep.alep_color import alep_colormap, green_yellow_red
    from alinea.adel.mtg_interpreter import plot3d
    from openalea.plantgl.all import Viewer
    # Compute severity by leaf
    severity_by_leaf = compute_severity_by_leaf(g, label=label)
    set_property_on_each_id(g, 'severity', severity_by_leaf, label=label)

    # Visualization
    g = alep_colormap(g, 'severity', cmap=green_yellow_red(levels=100),
                      lognorm=False, zero_to_one=False, vmax=100)

    if senescence==True:
        leaves = get_leaves(g, label=label)
        # pos_sen = g.property('position_senescence')
        sen_lengths = g.property('senesced_length')
        green_lengths = g.property('green_length')
        for leaf in leaves:
            if sen_lengths[leaf]>0. and round(green_lengths[leaf],15)==0.:
                g.node(leaf).color = (157, 72, 7)
    
    if transparency!=None:
        for id in g:
            if not id in severity_by_leaf:
                g.node(id).color = (255,255,255)
                g.node(id).transparency = 0.9
            elif severity_by_leaf[id]==0.:
                g.node(id).color = (255,255,255)
                g.node(id).transparency = transparency
            else:
                g.node(id).transparency = 0.
        
        scene = plot3d_transparency(g)
    else:
        scene = plot3d(g)
    Viewer.display(scene)
    return scene

def plot_severity_vine(g, trunk=True, transparency=None, label='lf'):
    from alinea.alep.architecture import set_property_on_each_id
    from alinea.alep.disease_outputs import compute_severity_by_leaf
    from alinea.alep.alep_color import alep_colormap, green_yellow_red
    from alinea.adel.mtg_interpreter import plot3d
    from openalea.plantgl.all import Viewer
    # Compute severity by leaf
    severity_by_leaf = compute_severity_by_leaf(g, label = label)
    set_property_on_each_id(g, 'severity', severity_by_leaf, label = label)
                       
    # Visualization
    g = alep_colormap(g, 'severity', cmap=green_yellow_red(levels=100),
                      lognorm=False, zero_to_one=False, vmax=100)
    brown = (100,70,30)
    if trunk==True:
        trunk_ids = [n for n in g if g.label(n).startswith('tronc')]
        for id in trunk_ids:
            trunk = g.node(id)
            trunk.color = brown
            
    if transparency!=None:
        for id in g:
            if not id in severity_by_leaf:
                g.node(id).color = (255,255,255)
                g.node(id).transparency = 0.9
            elif severity_by_leaf[id]==0.:
                g.node(id).color = (255,255,255)
                g.node(id).transparency = transparency
            else:
                g.node(id).transparency = 0.
        scene = plot3d_transparency(g)
    else:
        scene = plot3d(g)
    Viewer.display(scene)
    return scene
    
def save_image(scene, image_name='%s/img%04d.%s', directory='.', index=0, ext='png'):
    '''
    Save an image of a scene in a specific directory

    Parameters
    ----------

        - scene: a PlantGL scene
        - image_name: a string template 
            The format of the string is dir/img5.png
        - directory (optional: ".") the directory where the images are written
        - index: the index of the image
        - ext : the image format

    Example
    -------

        - Movie:
            convert *.png movie.mpeg
            convert *.png movie.gif
            mencoder "mf://*.png" -mf type=png:fps=25 -ovc lavc -o output.avi
            mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq "mf://*.png" -mf type=png:fps=18 -of avi  -o output.avi
            
    '''
    from openalea.plantgl.all import Viewer
    import os.path
    if not image_name:
        image_name='{directory}/img{index:0>4d}.{ext}'
    filename = image_name.format(directory=directory, index=index, ext=ext)
    Viewer.frameGL.saveImage(filename)
    return scene,
 
######################################################################
from numpy import mean
import numpy as np
from scipy.integrate import trapz

class VineLeafInspector:
    def __init__(self, leaf_id, label='lf'):
        self.leaf_id = leaf_id
        self.label = label
        # Initialize leaf properties to save
        self.leaf_area = []
        self.leaf_green_area = []  
        self.leaf_healthy_area = []
        self.leaf_disease_area = []
        # Initialize surfaces in state
        self.surface_latent = []
        self.surface_spo = []
        self.surface_empty = []
        # Initialize ratios (surfaces in state compared to leaf area)
        self.ratio_latent = []
        self.ratio_spo = []
        self.ratio_empty = []
        # Initialize total severity
        self.severity = []

    def update_data(self, g):
        leaf = g.node(self.leaf_id)
        area = leaf.area
        if area!=None:
            self.leaf_area.append(area)
            self.leaf_green_area.append(leaf.green_area)
            self.leaf_healthy_area.append(leaf.healthy_area)
            self.leaf_disease_area.append(area - leaf.healthy_area)
            
        else:
            area = 0.
            self.leaf_area.append(0.)
            self.leaf_green_area.append(0.)
            self.leaf_healthy_area.append(0.)
            self.leaf_disease_area.append(0.)
        
        self.severity.append(100.*(1.-leaf.healthy_area/area) if area>0. else 0.)
        
        try:
            lesions = leaf.lesions
        except:
            lesions = []
        surface_latent = 0.
        surface_spo = 0.
        surface_empty = 0.
        if len(lesions)>0.:
            latent_lesions = [l for l in lesions if l.is_latent()]
            if len(latent_lesions)>0.:
                surface_latent = sum([l.surface for l in latent_lesions])
            
            spo_lesions = [l for l in lesions if l.is_sporulating()]
            if len(spo_lesions)>0.:
                surface_spo = sum([l.surface for l in spo_lesions])
            
            empty_lesions = [l for l in lesions if l.is_empty()]
            if len(empty_lesions)>0.:
                surface_empty = sum([l.surface for l in empty_lesions])
                
        self.surface_latent.append(surface_latent)
        self.ratio_latent.append(100.*surface_latent/area if area>0. else 0.)
        self.surface_spo.append(surface_spo)
        self.ratio_spo.append(100.*surface_spo/area if area>0. else 0.)
        self.surface_empty.append(surface_empty)
        self.ratio_empty.append(100.*surface_empty/area if area>0. else 0.)

###############################################################################
class AdelWheatRecorder(object):
    """ Record simulation output on every leaf of main stems in a dataframe during simulation """
    def __init__(self, group_dus = True, 
                 fungus_name = 'template', increment = 1000):
        self.fungus_name = fungus_name
        self.group_dus = group_dus
        self.increment = increment
        columns = ['date', 'degree_days', 'num_plant', 'num_leaf_bottom', 
                   'leaf_area', 'leaf_green_area', 'fnl']
        self.data = pandas.DataFrame(data = [[np.nan for col in columns] 
                                    for i in range(increment)], 
                                    columns = columns)
    
    def get_values_single_leaf(self, g, date, degree_days, id_list):
        dict_lf = {}
        dict_lf['date'] = date
        dict_lf['degree_days'] = degree_days
        
        # Update leaf properties
        areas = g.property('area')
        green_areas = g.property('green_area')
        fnls = g.property('nff')
        a_label_splitted = self.a_labels[id_list[0]].split('_')
        dict_lf['num_plant'] = int(a_label_splitted[0].split('plant')[1])
        dict_lf['num_leaf_bottom'] = int(a_label_splitted[2].split('metamer')[1])
        dict_lf['leaf_area'] = sum([areas[id] for id in id_list])
        dict_lf['leaf_green_area'] = sum([green_areas[id] for id in id_list])
        dict_lf['fnl'] =  fnls[g.complex_at_scale(id_list[0], 2)]
        return dict_lf

    def increment_data(self):
        df = pandas.DataFrame(data = [[np.nan for col in self.data.columns] 
                                for i in range(self.increment)], 
                                columns = self.data.columns)
        self.data = pandas.concat([self.data, df])
        self.data = self.data.reset_index(drop = True)
        
    def add_line_from_dict(self, dict_lf={}):
        indx = self.data[self.data['date'].isnull()].index[0]
        self.data.loc[indx, :] = pandas.Series(dict_lf)
        
    def get_ids_on_blade(self, components, geometries=None, 
                         areas=None, labels=None):
        return [vid for vid in components if geometries.get(vid) is not None 
                                        and areas.get(vid) is not None
                                        and labels.get(vid).startswith('LeafElement')]
                                                
    def record(self, g, date = None, degree_days = None):
        self.a_labels = {vid:lab for vid, lab
                         in adel_labels(g, scale = 5).iteritems() 
                         if 'LeafElement' in lab}
        v_length = g.property('visible_length')
        labels = g.property('label')
        geometries = g.property('geometry')
        areas = g.property('area')
        blades = [id for id,lb in labels.iteritems() if lb.startswith('blade') and v_length[id]>0]
        for blade in blades:
            id_list = self.get_ids_on_blade(g.components(blade),
                                            geometries, areas, labels)
            if len(id_list)>0 and 'MS' in self.a_labels[id_list[0]]:
                dict_lf = self.get_values_single_leaf(g = g, date = date, 
                                                      degree_days = degree_days, 
                                                      id_list = id_list)
                self.add_line_from_dict(dict_lf)
                if len(self.data[self.data['date'].isnull()]) == 0.:
                    self.increment_data()
                    
    def add_leaf_numbers(self):
        self.data['axis'] = 'MS'
        for pl in set(self.data['num_plant']):
            df_ax = self.data[self.data['num_plant'] == pl]
            fnl = df_ax['fnl'].max()
            df_ax.loc[:, 'num_leaf_top'] = fnl - df_ax['num_leaf_bottom'] + 1
            self.data.loc[df_ax.index, 'num_leaf_top'] = fnl - df_ax['num_leaf_bottom'] + 1
            for date in set(df_ax['degree_days']):
                df_date = df_ax[df_ax['degree_days'] == date]
                current_max_bottom = df_date['num_leaf_bottom'].max()
                cur_max_leaf_top = df_date['fnl'] - current_max_bottom + 1
                self.data.loc[df_date.index, 'cur_max_leaf_top'] = cur_max_leaf_top
                self.data.loc[df_date.index, 'cur_num_leaf_top'] = cur_max_leaf_top - df_date['num_leaf_top'] + 1
                
    def add_variety(self, variety = None):
        self.data['variety'] = variety
    
    def post_treatment(self, variety = None):
        self.data = self.data[~pandas.isnull(self.data['date'])]
        self.add_leaf_numbers()
        if variety is not None:
            self.add_variety(variety=variety)
           
    def save(self, filename):
        self.data.to_csv(filename, index = False)
        
######################################################################
import numpy
import pandas
from scipy.integrate import simps
try:
    import cPickle as pickle
except:
    import pickle
from alinea.adel.newmtg import adel_labels

class AdelSeptoRecorder(AdelWheatRecorder):
    """ Record simulation output on every leaf of main stems in a dataframe during simulation """
    def __init__(self, group_dus = True, 
                 fungus_name = 'septoria', 
                 increment = 1000):
        super(AdelSeptoRecorder, self).__init__(group_dus = group_dus, 
                                                fungus_name = fungus_name,
                                                increment = increment)
        columns = ['date', 'degree_days', 'num_plant', 'num_leaf_bottom', 'leaf_area', 
                   'leaf_green_area', 'leaf_length', 'leaf_senesced_length', 'fnl', 
                   'nb_dispersal_units', 'nb_lesions', 'nb_lesions_on_green', 
                   'surface_inc', 'surface_chlo', 'surface_nec', 'surface_nec_on_green', 
                   'surface_spo', 'surface_spo_on_green', 'surface_empty', 
                   'surface_empty_on_green', 'surface_dead']
        self.data = pandas.DataFrame(data = [[np.nan for col in columns] for i in range(self.increment)], 
                                     columns = columns)
    
    def get_values_single_leaf(self, g, date, degree_days, id_list):
        dict_lf = {}
        dict_lf['date'] = date
        dict_lf['degree_days'] = degree_days
        
        # Update leaf properties
        areas = g.property('area')
        green_areas = g.property('green_area')
        lengths = g.property('length')
        senesced_lengths = g.property('senesced_length')
        fnls = g.property('nff')
        a_label_splitted = self.a_labels[id_list[0]].split('_')
        dict_lf['num_plant'] = int(a_label_splitted[0].split('plant')[1])
        dict_lf['num_leaf_bottom'] = int(a_label_splitted[2].split('metamer')[1])
        dict_lf['leaf_area'] = sum([areas[id] for id in id_list])
        dict_lf['leaf_green_area'] = sum([green_areas[id] for id in id_list])
        dict_lf['leaf_length'] = sum([lengths[id] for id in id_list])
        dict_lf['leaf_senesced_length'] = sum([senesced_lengths[id] for id in id_list])
        dict_lf['fnl'] =  fnls[g.complex_at_scale(id_list[0], 2)]

        # Update properties of dispersal units and lesions
        nb_dus = 0
        nb_lesions = 0
        nb_lesions_on_green = 0
        surface_inc = 0.
        surface_chlo = 0.
        surface_nec = 0.
        surface_nec_on_green = 0.
        surface_spo = 0.
        surface_spo_on_green = 0.
        surface_empty = 0.
        surface_empty_on_green = 0.
        surface_dead = 0.
        
        for id in id_list:
            leaf = g.node(id)
            if 'dispersal_units' in leaf.properties():
                for du in leaf.dispersal_units:
                    if du.fungus.name == self.fungus_name:
                        if self.group_dus:
                            nb_dus += du.nb_dispersal_units
                        else:
                            nb_dus += 1
                                
            if 'lesions' in leaf.properties():
                for les in leaf.lesions:
                    if les.fungus.name == self.fungus_name:
                        if self.group_dus:
                            nb_les = les.nb_lesions
                            nb_les_on_green = les.nb_lesions_non_sen
                            nb_lesions += nb_les
                            nb_lesions_on_green += nb_les_on_green
                            ratio_green = float(nb_les_on_green)/nb_les if nb_les>0. else 0.
                            surface_nec_on_green += les.surface_nec *  ratio_green
                            surface_spo_on_green += les.surface_spo * ratio_green
                            surface_empty_on_green += les.surface_empty * ratio_green
                        else:
                            nb_lesions += 1
                            if les.position[0][0]>leaf.senesced_length:
                                nb_lesions_on_green += 1
                                surface_nec_on_green = les.surface_nec_on_green
                                surface_spo_on_green = les.surface_spo
                                surface_empty_on_green = les.surface_empty
                        surface_inc += les.surface_inc
                        surface_chlo += les.surface_chlo
                        surface_nec += les.surface_nec
                        surface_spo += les.surface_spo
                        surface_empty += les.surface_empty
                        surface_dead += les.surface_dead
        
        dict_lf['nb_dispersal_units'] = nb_dus
        dict_lf['nb_lesions'] = nb_lesions
        dict_lf['nb_lesions_on_green'] = nb_lesions_on_green
        dict_lf['surface_inc'] = surface_inc
        dict_lf['surface_chlo'] = surface_chlo
        dict_lf['surface_nec'] = surface_nec
        dict_lf['surface_spo'] = surface_spo
        dict_lf['surface_empty'] = surface_empty
        dict_lf['surface_dead'] = surface_dead
        dict_lf['surface_nec_on_green'] = surface_nec_on_green
        dict_lf['surface_spo_on_green'] = surface_spo_on_green
        dict_lf['surface_empty_on_green'] = surface_empty_on_green
        return dict_lf
                
    def leaf_senesced_area(self):
        self.data['leaf_senesced_area'] = self.data['leaf_area'] - \
                                          self.data['leaf_green_area']
    
    def leaf_disease_area(self):
        self.data['leaf_disease_area'] = self.data['surface_inc'] + \
                                         self.data['surface_chlo'] + \
                                         self.data['surface_nec'] + \
                                         self.data['surface_spo'] + \
                                         self.data['surface_empty']
                                         
    def leaf_necrotic_area(self):
        self.data['leaf_necrotic_area'] = self.data['surface_nec'] + \
                                          self.data['surface_spo'] + \
                                          self.data['surface_empty']
        

    def leaf_necrotic_area_on_green(self):
        self.data['leaf_necrotic_area_on_green'] = self.data['surface_nec_on_green'] + \
                                                   self.data['surface_spo_on_green'] + \
                                                   self.data['surface_empty_on_green']

    def _ratio(self, variable='leaf_necrotic_area', against='leaf_area'):
        r = []
        for ind in self.data.index:
            a = self.data[against][ind]
            if a > 0.:
                r.append(self.data[variable][ind]/a)
            else:
                r.append(0.)
        return r
    
    def severity(self):
        """ Necrotic area of lesions compared to total leaf area """
        if not 'leaf_necrotic_area' in self.data:
            self.leaf_necrotic_area()
        self.data['severity'] = self._ratio(variable='leaf_necrotic_area',
                                            against='leaf_area')
        
    def severity_on_green(self):
        """ Necrotic area of lesions on green compared to green leaf area """
        if not 'leaf_necrotic_area_on_green' in self.data:
            self.leaf_necrotic_area_on_green()
        self.data['severity_on_green'] = self._ratio(variable='leaf_necrotic_area_on_green',
                                                     against='leaf_green_area')
    
    def ratios(self):
        self.data['ratio_inc'] = self._ratio(variable='surface_inc',
                                             against='leaf_area')
        self.data['ratio_chlo'] = self._ratio(variable='surface_chlo',
                                              against='leaf_area')
        self.data['ratio_nec'] = self._ratio(variable='surface_nec',
                                             against='leaf_area')
        self.data['ratio_spo'] = self._ratio(variable='surface_spo',
                                             against='leaf_area')
        self.data['ratio_empty'] = self._ratio(variable='surface_empty',
                                               against='leaf_area')

    def get_audpc(self, variable='severity'):
        for pl in set(self.data['num_plant']):
            df_pl =  self.data[self.data['num_plant'] == pl]
            for lf in set(df_pl['num_leaf_top']):
                df_lf = df_pl[df_pl['num_leaf_top'] == lf]
                ind_data_lf = df_lf.index
                if round(df_lf['leaf_green_area'][pandas.notnull(df_lf['leaf_disease_area'])].iloc[-1],10)==0.:
                    data = df_lf[variable][df_lf['leaf_green_area']>0]
                    ddays = df_lf['degree_days'][df_lf['leaf_green_area']>0]
                    data_ref = numpy.ones(len(data))
                    if len(data[data>0])>0:
                        audpc = simps(data[data>0], ddays[data>0])
                    else:
                        audpc = 0.
                    self.data.loc[ind_data_lf, 'audpc'] = audpc
                    audpc_ref = simps(data_ref[data_ref>0], ddays[data_ref>0])
                    self.data.loc[ind_data_lf, 'normalized_audpc'] = audpc/audpc_ref if audpc_ref>0. else 0.
                else:
                    self.data.loc[ind_data_lf, 'audpc'] = np.nan
                    self.data.loc[ind_data_lf, 'normalized_audpc'] = np.nan    
    
    def post_treatment(self, variety = None):
        self.data = self.data[~pandas.isnull(self.data['date'])]
        self.add_leaf_numbers()
        self.leaf_senesced_area()
        self.leaf_necrotic_area()
        self.leaf_necrotic_area_on_green()
        self.leaf_disease_area()
        self.ratios()
        self.severity()
        self.severity_on_green()
        self.get_audpc()
        if variety is not None:
            self.add_variety(variety=variety)
           
def get_recorder(*filenames):
    recorder = []
    for file in filenames:
        f_rec = open(file)
        recorder.append(pickle.load(f_rec))
        f_rec.close()
    return recorder if len(recorder)>1 else recorder[0]

###############################################################################
class BrownRustRecorder(AdelWheatRecorder):
    """ Record simulation output on every leaf of main stems in a dataframe during simulation """
    def __init__(self, group_dus = True, 
                 fungus_name = 'brown_rust', 
                 increment = 1000):
        super(BrownRustRecorder, self).__init__(group_dus = group_dus, 
                                                fungus_name = fungus_name,
                                                increment = increment)
        columns = ['date', 'degree_days', 'num_plant', 'num_leaf_bottom', 'leaf_area', 
                   'leaf_green_area', 'leaf_length', 'leaf_senesced_length', 'fnl', 
                   'nb_dispersal_units', 'nb_lesions',
                   'surface_sink', 'surface_chlo', 'surface_spo', 
                   'surface_empty', 'surface_dead']
        self.data = pandas.DataFrame(data = [[np.nan for col in columns] 
                                     for i in range(self.increment)],
                                     columns = columns)
    
    def get_values_single_leaf(self, g, date, degree_days, id_list):
        dict_lf = {}
        dict_lf['date'] = date
        dict_lf['degree_days'] = degree_days
        
        # Update leaf properties
        areas = g.property('area')
        green_areas = g.property('green_area')
        lengths = g.property('length')
        senesced_lengths = g.property('senesced_length')
        fnls = g.property('nff')
        a_label_splitted = self.a_labels[id_list[0]].split('_')
        dict_lf['num_plant'] = int(a_label_splitted[0].split('plant')[1])
        dict_lf['num_leaf_bottom'] = int(a_label_splitted[2].split('metamer')[1])
        dict_lf['leaf_area'] = sum([areas[id] for id in id_list])
        dict_lf['leaf_green_area'] = sum([green_areas[id] for id in id_list])
        dict_lf['leaf_length'] = sum([lengths[id] for id in id_list])
        dict_lf['leaf_senesced_length'] = sum([senesced_lengths[id] for id in id_list])
        dict_lf['fnl'] =  fnls[g.complex_at_scale(id_list[0], 2)]

        # Update properties of dispersal units and lesions
        nb_dus = 0
        nb_lesions = 0
        surface_sink = 0.
        surface_chlo = 0.
        surface_spo = 0.
        surface_empty = 0.
        surface_dead = 0.
        
        for id in id_list:
            leaf = g.node(id)
            if 'dispersal_units' in leaf.properties():
                for du in leaf.dispersal_units:
                    if du.fungus.name == self.fungus_name:
                        if self.group_dus:
                            nb_dus += du.nb_dispersal_units
                        else:
                            nb_dus += 1
                                
            if 'lesions' in leaf.properties():
                for les in leaf.lesions:
                    if les.fungus.name == self.fungus_name:
                        if self.group_dus:
                            nb_lesions += les.nb_lesions
                        else:
                            nb_lesions += 1
                        surface_sink += les.surface_sink
                        surface_chlo += les.surface_chlo
                        surface_spo += les.surface_spo
                        surface_empty += les.surface_empty
                        surface_dead += les.surface_dead
        
        dict_lf['nb_dispersal_units'] = nb_dus
        dict_lf['nb_lesions'] = nb_lesions
        dict_lf['surface_sink'] = surface_sink
        dict_lf['surface_chlo'] = surface_chlo
        dict_lf['surface_spo'] = surface_spo
        dict_lf['surface_empty'] = surface_empty
        dict_lf['surface_dead'] = surface_dead
        return dict_lf
    
    def leaf_senesced_area(self):
        self.data['leaf_senesced_area'] = self.data['leaf_area'] - self.data['leaf_green_area']
    
    def leaf_disease_area(self):
        self.data['leaf_disease_area'] = self.data['surface_sink'] + \
                                         self.data['surface_chlo'] + \
                                         self.data['surface_spo'] + \
                                         self.data['surface_empty'] + \
                                         self.data['surface_dead']

    def surface_alive(self):
        self.data['surface_alive'] = self.data['surface_sink'] + \
                                         self.data['surface_chlo'] + \
                                         self.data['surface_spo'] + \
                                         self.data['surface_empty']

    def add_max_leaf_green_area(self):
        for ind in self.data.index:
            df = self.data[(self.data['num_plant']==self.data['num_plant'][ind]) & \
                            (self.data['num_leaf_bottom']==self.data['num_leaf_bottom'][ind])]
            self.data.loc[ind,'max_leaf_green_area'] = df['leaf_green_area'].max()
            
    def get_ratio(self, variable='surface_alive'):
        output = []
        df = self.data
        for ind in df.index:
            if round(df['leaf_area'][ind], 16)>0.:
                output.append(self.data[variable][ind]/self.data['leaf_area'][ind])
            else:
                output.append(np.nan)
        return output

    def ratios(self):
        if not 'max_leaf_green_area' in self.data.columns:
            self.add_max_leaf_green_area()
        self.data['ratio_sink'] = self.get_ratio(variable='surface_sink')
        self.data['ratio_chlo'] = self.get_ratio(variable='surface_chlo')
        self.data['ratio_spo'] = self.get_ratio(variable='surface_spo')
        self.data['ratio_empty'] = self.get_ratio(variable='surface_empty')

    def severity(self):
        if not 'leaf_disease_area' in self.data:
            self.leaf_disease_area()
        if not 'max_leaf_green_area' in self.data.columns:
            self.add_max_leaf_green_area()
        self.data['severity'] = self.get_ratio(variable='surface_alive')

    def add_variety(self, variety = None):
        self.data['variety'] = variety
        
    def get_audpc(self, variable='severity'):
        for pl in set(self.data['num_plant']):
            df_pl =  self.data[self.data['num_plant'] == pl]
            for lf in set(df_pl['num_leaf_top']):
                df_lf = df_pl[df_pl['num_leaf_top'] == lf]
                ind_data_lf = df_lf.index
                if round(df_lf['leaf_green_area'][pandas.notnull(df_lf['leaf_disease_area'])].iloc[-1],10)==0.:
                    data = df_lf[variable][df_lf['leaf_green_area']>0]
                    ddays = df_lf['degree_days'][df_lf['leaf_green_area']>0]
                    data_ref = numpy.ones(len(data))
                    if len(data[data>0])>0:
                        audpc = simps(data[data>0], ddays[data>0])
                    else:
                        audpc = 0.
                    self.data.loc[ind_data_lf, 'audpc'] = audpc
                    audpc_ref = simps(data_ref[data_ref>0], ddays[data_ref>0])
                    self.data.loc[ind_data_lf, 'normalized_audpc'] = audpc/audpc_ref if audpc_ref>0. else 0.
                else:
                    self.data.loc[ind_data_lf, 'audpc'] = np.nan
                    self.data.loc[ind_data_lf, 'normalized_audpc'] = np.nan 
    
    def post_treatment(self, variety = None):
        self.data = self.data[~pandas.isnull(self.data['date'])]
        self.add_leaf_numbers()
        self.leaf_senesced_area()
        self.leaf_disease_area()
        self.surface_alive()
        self.ratios()
        self.severity()
        self.get_audpc()
        if variety is not None:
            self.add_variety(variety=variety)
        self.data['axis'] = 'MS'

def get_data_without_death(data, num_leaf = 'num_leaf_bottom'):
    datas = []
    for lf in set(data[num_leaf]):
        df = data[data[num_leaf] == lf]
        df_count = df.groupby('degree_days').count()
        last_date = df_count[df_count['num_plant'] == max(df_count['num_plant'])].index[-1]
        datas.append(df[df['degree_days']<=last_date])
    return pandas.concat(datas)
    
###############################################################################
class SeptoRustRecorder:
    """ Record simulation output on every leaf of main stems in a dataframe 
        during simulation of septoria and brown rust coupled epidemics """
    def __init__(self, group_dus = True, 
                 increment = 1000):
        super(AdelSeptoRecorder, self).__init__(group_dus = group_dus, 
                                                increment = increment)
        columns = ['date', 'degree_days', 'num_plant', 'num_leaf_bottom', 
                   'fnl', 'leaf_area', 'leaf_green_area', 
                   'surface_septo', 'surface_rust']
        self.data = pandas.DataFrame(data = [[np.nan for col in columns] for i in range(1000)], 
                                     columns = columns)
    
    def get_values_single_leaf(self, g, date, degree_days, id_list):
        dict_lf = {}
        dict_lf['date'] = date
        dict_lf['degree_days'] = degree_days
        
        # Update leaf properties
        areas = g.property('area')
        green_areas = g.property('green_area')
        fnls = g.property('nff')
        a_label_splitted = self.a_labels[id_list[0]].split('_')
        dict_lf['num_plant'] = int(a_label_splitted[0].split('plant')[1])
        dict_lf['num_leaf_bottom'] = int(a_label_splitted[2].split('metamer')[1])
        dict_lf['leaf_area'] = sum([areas[id] for id in id_list])
        dict_lf['leaf_green_area'] = sum([green_areas[id] for id in id_list])
        dict_lf['fnl'] =  fnls[g.complex_at_scale(id_list[0], 2)]

        # Update properties of dispersal units and lesions
        surface_septo = 0.
        surface_rust = 0.
        
        for id in id_list:
            leaf = g.node(id)
            if 'lesions' in leaf.properties():
                for les in leaf.lesions:
                    if les.fungus.name == 'septoria':
                        surface_septo += les.surface
                    elif les.fungus.name == 'brown_rust':
                        surface_rust += les.surface
        
        dict_lf['surface_septo'] = surface_septo
        dict_lf['surface_rust'] = surface_rust
        return dict_lf
    
    def leaf_senesced_area(self):
        self.data['leaf_senesced_area'] = self.data['leaf_area'] - self.data['leaf_green_area']
    
    def leaf_disease_area(self):
        self.data['leaf_disease_area'] = self.data['surface_septo'] + self.data['surface_rust']

    def severity(self):
        if not 'leaf_disease_area' in self.data:
            self.leaf_disease_area()
        self.data['severity'] = [self.data['leaf_disease_area'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['severity_septo'] = [self.data['surface_septo'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['severity_rust'] = [self.data['surface_rust'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]

    def add_variety(self, variety = None):
        self.data['variety'] = variety
    
    def post_treatment(self, variety = None):
        self.data = self.data[~pandas.isnull(self.data['date'])]
        self.add_variety(variety = variety)        
        self.add_leaf_numbers()
        self.leaf_senesced_area()
        self.leaf_disease_area()
        self.severity()
        if variety is not None:
            self.add_variety()
        self.data['axis'] = 'MS'
    
# Plotting ####################################################################
import matplotlib.pyplot as plt
from math import sqrt

def variance(lst):
    """
    Uses standard variance formula (sum of each (data point - mean) squared)
    all divided by number of data points
    """
    mu = numpy.mean(lst)
    return 1.0/(len(lst)-1) * sum([(i-mu)**2 for i in lst])
        
def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean.
    The constant factor can be looked up from a table, for 95pcent confidence
    on a reasonable size sample (>=500) 1.96 is used.
    """
    from scipy.stats import t
    
    n, v = len(lst), variance(lst)
    c = t.interval(perc_conf * 1.0 / 100, n-1)[1]
    
    return sqrt(v/n) * c

def plot_mean(data, variable = 'severity', xaxis = 'degree_days', 
              error_bars = False, error_method = 'confidence_interval', 
              marker = 'd', empty_marker = False, linestyle = '-', color = 'b', 
              alpha = None, title = None, xlabel = None, ylabel = None,
              xlims = None, ylims = None, ax = None, return_ax = False):
    if variable in data.columns:
        if ax == None:
            fig, ax = plt.subplots()
        if empty_marker == False:
            markerfacecolor = color
        else:
            markerfacecolor = 'none'
        if alpha is None:
            alpha = 1
            
        df = data[pandas.notnull(data.loc[:,variable])].loc[:, [xaxis, variable]]
        df_mean = df.groupby(xaxis).mean()
        df['nb_rep'] = map(lambda x: df[xaxis].value_counts()[x], df[xaxis])
        if error_bars == True and len(df['nb_rep'])>0 and min(df['nb_rep'])>1:
            if error_method == 'confidence_interval':
                df_err = df.groupby(xaxis).agg(conf_int)
            elif error_method == 'std_deviation':
                df_err = df.groupby(xaxis).std()
            else:
                raise ValueError("'error_method' unknown: 'try confidence_interval' or 'std_deviation'")
            ax.errorbar(df_mean.index, df_mean[variable], yerr = df_err[variable].values,
                        marker = marker, linestyle = linestyle, color = color,
                        markerfacecolor = markerfacecolor,  markeredgecolor = color)
        else:
            ax.plot(df_mean.index, df_mean[variable],
                    marker = marker, linestyle = linestyle, color = color,
                    markerfacecolor = markerfacecolor,  markeredgecolor = color)
        if title is not None:
            ax.set_title(title, fontsize = 18)
        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize = 18)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize = 18)
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        if return_ax == True:
            return ax
            
def plot_sum(data, variable = 'severity', xaxis = 'degree_days', 
              marker = 'd', linestyle = '-', color = 'b', alpha = None,
              title = None, xlabel = None, ylabel = None,
              xlims = None, ylims = None, ax = None, return_ax = False):
    
    def get_mean_data(variable):
        if variable in data.columns:
            df_sum.loc[:, variable] = df_mean.loc[:, variable]    
    
    if alpha is None:
        alpha = 1
    
    if variable in data.columns:
        if ax == None:
            fig, ax = plt.subplots()
        df = data[pandas.notnull(data.loc[:,variable])]
        df_mean = df.groupby(xaxis).mean()

        df_sum = df.groupby(xaxis).sum()
        for var in ['HS', 'num_leaf_bottom', 'num_leaf_top', 'fnl', 
                    'cur_max_leaf_top', 'cur_num_leaf_top']:
            if var != 'xaxis':
                get_mean_data(var)

        ax.plot(df_sum.index, df_sum[variable], color = color, alpha=alpha,
                    marker = marker, linestyle = linestyle)
        if title is not None:
            ax.set_title(title, fontsize = 18)
        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize = 18)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize = 18)
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        if return_ax == True:
            return ax

def plot_by_leaf(data, variable = 'green_area', xaxis = 'degree_days', 
                  leaves = range(1, 14), from_top = True, plant_axis = ['MS'],
                  error_bars = False, error_method = 'confidence_interval', 
                  marker = '', empty_marker = False, linestyle = '-', fixed_color = None, 
                  alpha = None, title = None, legend = True, xlabel = None, ylabel = None,
                  xlims = None, ylims = None, ax = None, return_ax = False, fig_size = (10,8)):
    df = data.copy()
    if ax == None:
        fig, ax = plt.subplots(figsize = fig_size)
    colors = ax._get_lines.set_color_cycle()
    colors = ax._get_lines.color_cycle
    if alpha is None:
        alpha = 1
        
    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'
    
    proxy = []
    labels = []
    for lf in leaves:
        df_lf = df[(df['axis'].isin(plant_axis)) & (df[num_leaf]==lf)]
        if fixed_color == None:
            color = next(colors)
        else:
            color = fixed_color
        plot_mean(df_lf, variable = variable, xaxis = xaxis, 
                  error_bars = error_bars, error_method = error_method, 
                  marker = marker, empty_marker = empty_marker, linestyle = linestyle, 
                  color = color, alpha = alpha, 
                  title = title, xlabel = xlabel, ylabel = ylabel,
                  xlims = xlims, ylims = ylims, ax = ax)
        proxy += [plt.Line2D((0,1),(0,0), color = color, linestyle ='-')]
        labels += ['L%d' %lf]
    
    if legend == True:
        colors = ax._get_lines.set_color_cycle()
        colors = ax._get_lines.color_cycle
        ax.legend(proxy, labels, title = 'Leaf\nNumber',
                    loc='center left', bbox_to_anchor=(1, 0.5))
    if return_ax == True:
        return ax