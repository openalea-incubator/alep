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
        sev_bl = np.zeros(len(les_bl))
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
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
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
    from alinea.alep.architecture import get_leaves
    vids = get_leaves(g, label=label)
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
        
######################################################################
import numpy
import pandas
from scipy.integrate import trapz, simps
from alinea.astk.plantgl_utils import get_height
from alinea.adel.newmtg import adel_ids
from collections import Iterable
try:
    import cPickle as pickle
except:
    import pickle
from alinea.adel.newmtg import adel_labels

class AdelSeptoRecorder:
    """ Record simulation output on every leaf of main stems in a dataframe during simulation """
    def __init__(self, group_dus = True, fungus_name = 'septoria'):
        self.fungus_name = fungus_name
        self.group_dus = group_dus
        columns = ['date', 'degree_days', 'num_plant', 'num_leaf_bottom', 'leaf_area', 
                   'leaf_green_area', 'leaf_length', 'leaf_senesced_length', 'fnl', 
                   'nb_dispersal_units', 'nb_dus_on_green', 'nb_lesions', 'nb_lesions_on_green', 
                   'surface_inc', 'surface_chlo', 'surface_nec', 'surface_spo', 
                   'surface_spo_on_green', 'surface_empty', 'surface_empty_on_green', 
                   'surface_dead']
        self.data = pandas.DataFrame(data = [[np.nan for col in columns] for i in range(20000)], 
                                     columns = columns)
        # self.data = pandas.DataFrame(columns = columns)
    
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
        nb_dus_on_green = 0
        nb_lesions = 0
        nb_lesions_on_green = 0
        surface_inc = 0.
        surface_chlo = 0.
        surface_nec = 0.
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
                            nb_dus += len(du.position)
                            nb_dus_on_green += len(filter(lambda x: x[0]>leaf.senesced_length, du.position))
                        else:
                            nb_dus += 1
                            if du.position[0][0]>leaf.senesced_length:
                                nb_dus_on_green += 1
                                
            if 'lesions' in leaf.properties():
                for les in leaf.lesions:
                    if les.fungus.name == self.fungus_name:
                        if self.group_dus:
                            nb_les = len(les.position)
                            nb_les_on_green = len(filter(lambda x: x[0]>leaf.senesced_length, les.position))
                            nb_lesions += nb_les
                            nb_lesions_on_green += nb_les_on_green
                            surface_spo_on_green = les.surface_spo * nb_les_on_green/nb_les if nb_les>0. else 0.
                            surface_empty_on_green = les.surface_empty * nb_les_on_green/nb_les if nb_les>0. else 0.
                        else:
                            nb_lesions += 1
                            if les.position[0][0]>leaf.senesced_length:
                                nb_lesions_on_green += 1
                                surface_spo_on_green = les.surface_spo
                                surface_empty_on_green = les.surface_empty
                        surface_inc += les.surface_inc
                        surface_chlo += les.surface_chlo
                        surface_nec += les.surface_nec
                        surface_spo += les.surface_spo
                        surface_empty += les.surface_empty
                        surface_dead += les.surface_dead
        
        dict_lf['nb_dispersal_units'] = nb_dus
        dict_lf['nb_dus_on_green'] = nb_dus_on_green
        dict_lf['nb_lesions'] = nb_lesions
        dict_lf['nb_lesions_on_green'] = nb_lesions_on_green
        dict_lf['surface_inc'] = surface_inc
        dict_lf['surface_chlo'] = surface_chlo
        dict_lf['surface_nec'] = surface_nec
        dict_lf['surface_spo'] = surface_spo
        dict_lf['surface_empty'] = surface_empty
        dict_lf['surface_dead'] = surface_dead
        dict_lf['surface_spo_on_green'] = surface_spo_on_green
        dict_lf['surface_empty_on_green'] = surface_empty_on_green
        return dict_lf
        
    def record(self, g, date = None, degree_days = None):
        self.a_labels = {vid:lab for vid, lab in adel_labels(g, scale = 5).iteritems() 
                            if 'LeafElement' in lab}
        v_length = g.property('visible_length')
        labels = g.property('label')
        geometries = g.property('geometry')
        areas = g.property('area')
        blades = [id for id,lb in labels.iteritems() if lb.startswith('blade') and v_length[id]>0]
        for blade in blades:
            id_list = [id for id in g.components(blade) if geometries.get(id) is not None 
                                                        and areas.get(id) is not None
                                                        and labels[id].startswith('LeafElement')]
            if len(id_list)>0 and 'MS' in self.a_labels[id_list[0]]:
                dict_lf = self.get_values_single_leaf(g = g, date = date, 
                                                      degree_days = degree_days, 
                                                      id_list = id_list)
                indx = self.data[self.data['date'].isnull()].index[0]
                self.data.loc[indx, :] = pandas.Series(dict_lf)
                if len(self.data[self.data['date'].isnull()]) == 0.:
                    df = pandas.DataFrame(data = [[np.nan for col in self.data.columns] 
                                            for i in range(20000)], columns = self.data.columns)
                    self.data = pandas.concat([self.data, df])
                    self.data = self.data.reset_index(drop = True)
                
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
    
    def leaf_senesced_area(self):
        self.data['leaf_senesced_area'] = self.data['leaf_area'] - self.data['leaf_green_area']
    
    def leaf_disease_area(self):
        self.data['leaf_disease_area'] = self.data['surface_inc'] + self.data['surface_chlo'] + self.data['surface_nec'] + self.data['surface_spo'] + self.data['surface_empty'] + self.data['surface_dead']

    def leaf_lesion_area_on_green(self):
        if not 'leaf_disease_area' in self.data:
            self.leaf_disease_area()
        if not 'leaf_senesced_area' in self.data:
            self.leaf_senesced_area()
        self.data['leaf_lesion_area_on_green'] = [self.data['leaf_disease_area'][ind]*(1 - self.data['leaf_senesced_area'][ind]/self.data['leaf_area'][ind]) if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
    
    def leaf_necrotic_area_on_green(self):
        if not 'necrosis' in self.data:
            self.necrosis()
        if not 'leaf_senesced_area' in self.data:
            self.leaf_senesced_area()
        self.data['leaf_necrotic_area_on_green'] = [self.data['necrosis'][ind]*(1 - self.data['leaf_senesced_area'][ind]/self.data['leaf_area'][ind]) if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
    
    def leaf_necrotic_senescent(self):
        if not 'leaf_necrotic_area_on_green' in self.data:
            self.leaf_necrotic_area_on_green()
        if not 'leaf_senesced_area' in self.data:
            self.leaf_senesced_area()
        self.data['leaf_necrotic_senescent'] = self.data['leaf_senesced_area'] + self.data['leaf_necrotic_area_on_green']
        
    def necrotic_senescent_percentage(self):
        if not 'leaf_necrotic_senescent' in self.data:
            self.leaf_necrotic_senescent()
        self.data['leaf_necrotic_senescent'] = [self.data['leaf_necrotic_senescent'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        
    def leaf_healthy_area(self):
        if not 'leaf_lesion_area_on_green' in self.data:
            self.leaf_lesion_area_on_green()
        self.data['leaf_healthy_area'] = self.data['leaf_green_area'] - self.data['leaf_lesion_area_on_green']
    
    def leaf_unhealthy_area(self):
        if not 'leaf_healthy_area' in self.data:
            self.leaf_healthy_area()
        self.data['leaf_unhealthy_area'] = self.data['leaf_area'] - self.data['leaf_healthy_area']

    def surface_alive(self):
        self.data['surface_alive'] = self.data['surface_inc'] + self.data['surface_chlo'] + self.data['surface_nec'] + self.data['surface_spo'] + self.data['surface_empty']

    def ratios(self):
        self.data['ratio_inc'] = [self.data['surface_inc'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_chlo'] = [self.data['surface_chlo'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_nec'] = [self.data['surface_nec'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_spo'] = [self.data['surface_spo'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_empty'] = [self.data['surface_empty'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_spo_on_green'] = [self.data['surface_spo_on_green'][ind]/self.data['leaf_green_area'][ind] if self.data['leaf_green_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_empty_on_green'] = [self.data['surface_empty_on_green'][ind]/self.data['leaf_green_area'][ind] if self.data['leaf_green_area'][ind]>0. else 0. for ind in self.data.index]

    def severity(self):
        if not 'leaf_disease_area' in self.data:
            self.leaf_disease_area()
        self.data['severity'] = [self.data['leaf_disease_area'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
    
    def severity_on_green(self):
        if not 'severity' in self.data:
            self.severity()
        if not 'leaf_lesion_area_on_green' in self.data:
            self.leaf_lesion_area_on_green()
        self.data['severity_on_green'] = [self.data['leaf_lesion_area_on_green'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]

    def necrosis(self):
        self.data['necrosis'] = self.data['surface_nec'] + self.data['surface_spo'] + self.data['surface_empty']
    
    def necrosis_percentage(self):
        if not 'necrosis' in self.data:
            self.necrosis()
        self.data['necrosis_percentage'] = [self.data['necrosis'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        
    def pycnidia_coverage(self):
        if not 'ratio_spo' in self.data and not 'ratio_empty' in self.data:
            self.ratios()
        self.data['pycnidia_coverage'] = self.data['ratio_spo'] + self.data['ratio_empty']
        
    def pycnidia_coverage_on_green(self):
        if not 'ratio_spo_on_green' in self.data and not 'ratio_empty_on_green' in self.data:
            self.ratios()
        self.data['pycnidia_coverage_on_green'] = self.data['ratio_spo_on_green'] + self.data['ratio_empty_on_green']

    def get_audpc(self, variable='pycnidia_coverage'):
        if not variable in self.data:
            exec('self.'+variable+'()')
        for pl in set(self.data['num_plant']):
            df_pl =  self.data[self.data['num_plant'] == pl]
            for lf in set(df_pl['num_leaf_top']):
                ind_data_lf = (self.data['num_plant'] == pl) & (self.data['num_leaf_top'] == lf)
                df_lf = self.data[ind_data_lf]
                if df_lf['leaf_green_area'][pandas.notnull(df_lf['leaf_green_area'])].iloc[-1]==0.:
                    data = df_lf[variable][df_lf['leaf_green_area']>0]
                    ddays = df_lf['degree_days'][df_lf['leaf_green_area']>0]
                    data_ref = numpy.ones(len(data))
                    audpc = simps(data, ddays)
                    self.data.loc[ind_data_lf, 'audpc_'+variable] = audpc
                    audpc_ref = simps(data_ref, ddays)
                    self.data.loc[ind_data_lf, 'normalized_audpc_'+variable] = audpc/audpc_ref if audpc_ref>0. else 0.
                else:
                    self.data.loc[ind_data_lf, 'audpc_'+variable] = 'audpc not available: leaf has not reached senescence or premature leaf death'
                    self.data.loc[ind_data_lf, 'normalized_audpc_'+variable] = 'audpc not available: leaf has not reached senescence or premature leaf death'
    
    def add_variety(self, variety = None):
        self.data['variety'] = variety
    
    def post_treatment(self, variety = None):
        self.add_leaf_numbers()
        self.leaf_senesced_area()
        self.leaf_disease_area()
        self.leaf_lesion_area_on_green()
        self.leaf_healthy_area()
        self.leaf_unhealthy_area()
        self.surface_alive()
        self.ratios()
        self.severity()
        self.severity_on_green()
        self.necrosis()
        self.necrosis_percentage()
        self.pycnidia_coverage()
        self.pycnidia_coverage_on_green()
        self.get_audpc('necrosis_percentage')
        self.get_audpc('pycnidia_coverage')
        if variety is not None:
            self.add_variety()
        self.data['axis'] = 'MS'
           
    def save(self, filename):
        self.post_treatment()
        self.data.to_csv(filename, index = False)
        
def initiate_all_adel_septo_recorders(g, nsect=5, date_sequence = None, fungus_name = 'septoria'):
    """ Used in the case of recording all blades of the main stem of each plant.
    
    Returns
    -------
    leaf_labels: dict('P1'=dict('F1'= recorder init with list of labels for leaf sectors,
                                ...,
                                'Fi'= recorder),
                      ...,
                      'Pn'=dict('F1' = recorder,
                                ...,
                                'Fi' = recorder))
    
    """
    vids = adel_ids(g)
    labels = g.property('label')
    stems = [id for id,lb in labels.iteritems() if lb.startswith('MS')]
    recorders = {}
    ind_plant = 0
    for st in stems:
        ind_plant += 1
        recorders['P%d' % ind_plant] = {}
        nff = int(g.node(st).properties()['nff'])
        ind_lf = nff+1
        for leaf in range(1, nff+1):
            ind_lf -= 1
            lf_labels = ['plant%d_MS_metamer%d_blade_LeafElement%d' % (ind_plant, leaf, sect) for sect in range(1, nsect+1)]
            recorders['P%d' % ind_plant]['F%d' % ind_lf] = AdelSeptoRecorder(adel_labels = lf_labels, 
                                                                             date_sequence = date_sequence,
                                                                             fungus_name = fungus_name)
            recorders['P%d' % ind_plant]['F%d' % ind_lf].update_vids_with_labels(vids)
    return recorders
        
def num_leaf_to_str(num_leaves=range(1,5)):
    return ['F%d' % lf for lf in num_leaves]
    
def get_recorder(*filenames):
    recorder = []
    for file in filenames:
        f_rec = open(file)
        recorder.append(pickle.load(f_rec))
        f_rec.close()
    return recorder if len(recorder)>1 else recorder[0]

def recorder_to_dataframe(recorder, weather = None, adel = None, skipna = True):
    """ Translate recorder object in DataFrame with the same format as disease notations """
    from alinea.astk.Weather import Weather
    from alinea.echap.disease.septo_data_reader import table_count_notations
    recos = []
    for pl, rec_pl in recorder.iteritems():
        for lf, rec_lf in rec_pl.iteritems():
            rec_lf.data['plant'] = int(pl[1:])
            rec_lf.data['num_leaf_top'] = int(lf[1:])
            rec_lf.data['num_leaf_bottom'] = len(rec_pl) - int(lf[1:]) + 1
            rec_lf.data['variety'] = 'tremie'
            rec_lf.data = rec_lf.data.rename(columns = {'date_sequence':'datetime'})
            # Temp : to move in disease outputs
            rec_lf.data['pycnidia_coverage'] = (rec_lf.data['ratio_spo'] + rec_lf.data['ratio_empty'])*100
            rec_lf.data['pycnidia_coverage_on_green'] = (rec_lf.data['ratio_spo'] + rec_lf.data['ratio_empty'])*100
            # rec_lf.get_normalized_audpc()
            # rec_lf.data['audpc'] = rec_lf.normalized_audpc
            recos.append(rec_lf.data)
    data_sim = pandas.concat(recos)
    data_sim.index.name = 'datetime'
    
    # Convert ratios in percentage
    list_of_ratios = [var for var in data_sim.columns if var.startswith('ratio')]+['severity', 'necrosis_percentage']
    data_sim[list_of_ratios] = data_sim[list_of_ratios].apply(lambda x: x*100.)
    
    # Add leaf dates
    if weather is not None:
        filename = find_dates_filename(weather)
        data_sim = add_leaf_dates_to_data(data_sim, adel, filename = filename)
    else:
        weather = Weather()
        weather.data = pandas.DataFrame(data = recorder['P1']['F1'].data.degree_days,
                                        index = recorder['P1']['F1'].data.index)
    
    # Ignore data from dates with dead leaves in each layer
    if skipna == True:
        df_count = table_count_notations(data_sim, weather, variable = 'severity', add_ddays = True)
        for lf in df_count.columns:
            df_lf = df_count[lf][map(lambda x: isinstance(x, (numpy.int64, int, float)), 
                                    df_count[lf])]
            nan_dates = df_lf[df_lf<df_lf.max()].reset_index().loc[:,'Date']
            if len(nan_dates)>0:
                for variable in data_sim.columns:
                    if variable not in ['datetime', 'degree_days', 'date_death',
                                        'variety', 'plant', 'num_leaf_top', 'num_leaf_bottom']:
                        data_sim[variable][(data_sim['num_leaf_top']==lf) & 
                                (data_sim.index.isin([d for d in nan_dates]))] = np.nan
    return data_sim
    
def split_recorder_by_fnl(recorder):
    fnls = set(len(v) for v in recorder.itervalues())
    recorders = {}
    for fnl in fnls:
        recorders[fnl] = {k:v for k,v in recorder.iteritems() if len(v)==fnl}
    return recorders
    
def renumber_recorder_from_bottom(recorder):
    bottomed_reco = {}
    for k, v in recorder.iteritems():
        fnl = len(v)
        bottomed_reco[k] = {}
        for kk, vv in v.iteritems():
            num_lf_top = int(''.join(x for x in kk if x.isdigit()))
            num_lf_bottom = fnl - num_lf_top + 1
            bottomed_reco[k]['F%d' % num_lf_bottom] = vv
    return bottomed_reco
    
def mean_by_leaf(recorder, variable='necrosis_percentage', skipna = False):
    ddays = max(recorder.values()[0].values(), key= lambda x: len(x.data.degree_days)).data.degree_days
    leaves = ['F%d' % leaf for leaf in range(1, max(len(v) for v in recorder.itervalues())+1)]
    df_mean_by_leaf = pandas.DataFrame(data={lf:[numpy.nan for i in range(len(ddays))] for lf in leaves}, 
                                        index = ddays, columns = leaves)
    dfs = []
    for lf in leaves:
        df_leaf = pandas.concat([v[lf].data[variable] for v in recorder.itervalues() if lf in v], axis=1)
        dfs.append(df_leaf)
        df_mean_by_leaf[:ddays[len(df_leaf)-1]][lf] = df_leaf.mean(axis=1, skipna = skipna).values
    return df_mean_by_leaf

def mean_audpc_by_leaf(recorder, variable = 'necrosis_percentage', normalized=True):
    def try_get(x, lf):
        try:
            if normalized==True:
                try:
                    return float(recorder[x][lf].normalized_audpc)
                except:
                    return float(recorder[x][lf].get_normalized_audpc(variable = variable))
            else:
                try:
                    return float(recorder[x][lf].audpc)
                except:
                    return float(recorder[x][lf].get_audpc(variable = variable))
        except:
            return np.nan
            
    df_mean_by_leaf = pandas.DataFrame()
    plants = recorder.keys()
    for leaf in range(1, max(len(v) for v in recorder.itervalues())+1):
        lf = 'F%d' % leaf
        df_mean_by_leaf[lf] = map(lambda x: try_get(x, lf), plants)
    return df_mean_by_leaf.mean()
    
def glue_df_means(df_means, nb_rep=5):
    glued = pandas.concat(df_means, axis=1, keys=range(nb_rep))
    glued.swaplevel(0, 1, axis=1).sortlevel(axis=1)
    return glued.groupby(level=1, axis=1).mean()
    
def get_mean_by_leaf(variable='necrosis_percentage', *recorders):
    if len(recorders)==1:
        return mean_by_leaf(recorders[0], variable=variable)
    else:
        df_means = [mean_by_leaf(reco, variable=variable) for reco in recorders]
        return glue_df_means(df_means, len(recorders))