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

class LeafInspector:
    def __init__(self, g, blade_id=None, label='LeafElement'):
        """ Find the ids of the leaf elements on the chosen blade.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        blade_id: int
            Index of the blade to be inspected
        label: str
            Label of the part of the MTG concerned by the calculation
        """
        from numpy import mean
        labels = g.property('label')
        self.label = label
        # Find leaf elements on the blade
        self.ids = [id for id in g.components(blade_id) if labels[id].startswith(label)]
        # Initialize leaf properties to save
        self.leaf_area = []
        self.leaf_green_area = []  
        self.leaf_healthy_area = []
        self.leaf_disease_area = []
        self.leaf_position_senescence = []
        self.leaf_senescence = []
        self.leaf_ratio_senescence = []
        self.leaf_green_lesions_area = []
        # Initialize variables relative to DUs
        self.nb_dus = []
        self.nb_dus_on_green = []
        self.nb_dus_on_healthy = []
        # Initialize variables relative to number of infections
        self.previous_nb_lesions = 0.
        self.nb_infections = []
        # Initialize variables relative to number of lesions
        self.nb_lesions = []
        self.nb_lesions_after_chlo = []
        # Initialize surfaces in state
        self.surface_inc = []
        self.surface_chlo = []
        self.surface_nec = []
        self.surface_spo = []
        self.surface_empty = []
        self.surface_total_nec = []
        self.surface_sen_nec = []
        # Initialize ratios (surfaces in state compared to leaf area)
        self.ratio_inc = []
        self.ratio_chlo = []
        self.ratio_nec = []
        self.ratio_spo = []
        self.ratio_empty = []
        self.ratio_total_nec = []
        self.ratio_sen_nec = []
        # Initialize total severity
        self.severity = []
        # Initialize necrosis percentage
        self.necrosis = []
        # Initialize number of spores in stock
        self.stock_spores = []
        self.nb_spores_emitted = []
                
        # Temp
        self.wetness = []
    
    def update_du_variables(self, g):
        """ Save counts of dispersal units.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy 
        """
        dus = 0.
        dus_on_green = 0.
        dus_on_healthy = 0.
        for id in self.ids:
            leaf = g.node(id)
            try:
                dus += len(leaf.dispersal_units)
                dus_on_green = len([du for du in leaf.dispersal_units if du.position[0] < leaf.position_senescence])
                dus_on_healthy = len([du for du in leaf.dispersal_units if du.can_infect_at_position])
                inactive_dus = len([du for du in leaf.dispersal_units if du.is_active == False])
            except:
                pass
            
            # Temp
            self.wetness.append(leaf.wetness)
            
        self.nb_dus.append(dus)
        self.nb_dus_on_green.append(dus_on_green)
        self.nb_dus_on_healthy.append(dus_on_healthy)
        
    def update_variables(self, g):
        """ Save leaf properties and disease properties.
        
        Save leaf properties and update the computation of severity, 
        necrosis percentage and ratios.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy   
        """
        area = 0.
        green_area = 0.
        disease_area = 0.
        lesion_list = []
        after_chlo = []
        for id in self.ids:
            leaf = g.node(id)
            try:
                lesion_list += leaf.lesions
                after_chlo += [les for les in leaf.lesions if les.status>=les.fungus.CHLOROTIC]
            except:
                pass

        self.nb_lesions.append(len(lesion_list))
        self.nb_lesions_after_chlo.append(len(after_chlo))
        
        if len(self.ids)==1:
            self.leaf_position_senescence.append(g.node(id).position_senescence)
            self.leaf_ratio_senescence.append(1-g.node(id).position_senescence)
       
        self.update_area(g)
        self.update_green_area(g)
        self.update_healthy_area(g)
        self.update_disease_area(g)
        self.update_green_lesions_area(g)
        self.compute_senescence(g)
        
        self.compute_nb_infections(g)
        self.compute_ratios(g)
        self.compute_severity(g)
        self.compute_necrosis(g)
        self.compute_senescence_necrosis(g)
        
        # self.compute_stock(g)
        # self.compute_nb_spores_emitted(g)
    
    def update_area(self,g):
        areas = g.property('area')
        self.leaf_area.append(sum([areas[id] for id in self.ids]))
    
    def update_green_area(self, g):
        green_areas = g.property('green_area')
        self.leaf_green_area.append(sum([green_areas[id] for id in self.ids]))
        
    def update_healthy_area(self, g):
        healthy_areas = g.property('healthy_area')
        self.leaf_healthy_area.append(sum([healthy_areas[id] for id in self.ids]))
    
    def compute_senescence(self, g):
        sen = compute_senescence_by_leaf(g, label=self.label)
        self.leaf_senescence.append(sum([sen[id] for id in self.ids]))
        
    def update_disease_area(self, g):
        disease_area = compute_lesion_areas_by_leaf(g, label=self.label)
        self.leaf_disease_area.append(sum([disease_area[id] for id in self.ids]))
    
    def update_green_lesions_area(self, g):
        green_lesion_areas = compute_green_lesion_areas_by_leaf(g, label=self.label)
        self.leaf_green_lesions_area.append(sum([green_lesion_areas[id] for id in self.ids]))
       
    def compute_nb_infections(self, g):
        """ Compute the number of infections during time step.
        
        Number of infections corresponds to number of new lesions created.        
        """
        lesions = g.property('lesions')
        nb_infections = 0.
        if len(lesions)>0:
            lesions = sum([lesions[id] for id in self.ids if id in lesions.keys()], [])
            if lesions!=0:
                nb_lesions = len(lesions)
                nb_infections = max(0., nb_lesions - self.previous_nb_lesions)
                self.previous_nb_lesions = nb_lesions

        self.nb_infections.append(nb_infections)
            
    def compute_severity(self, g):
        """ Compute severity on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        severities = compute_severity_by_leaf(g, label=self.label)
        self.severity.append(mean([severities[id] for id in self.ids]))
    
    def compute_necrosis(self, g):
        """ Compute necrosis percentage on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        nec = compute_necrosis_percentage_by_leaf(g, label=self.label)
        self.necrosis.append(mean([nec[id] for id in self.ids]))
    
    def compute_senescence_necrosis(self, g):
        nec_sen = compute_senescence_necrosis_by_leaf(g, label=self.label)
        ns = sum([nec_sen[id] for id in self.ids])
        self.surface_sen_nec.append(ns)
        
        total_area = sum([g.node(id).area for id in self.ids])
        self.ratio_sen_nec.append(100. * ns / total_area if total_area>0. else 0.)
    
    def compute_ratios(self, g):
        """ Compute surface of lesions in chosen state on a blade of the MTG.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy    
        """
        total_area = sum([g.node(id).area for id in self.ids])
        lesion_list = []
        for id in self.ids:
            leaf = g.node(id)
            try:
                lesion_list += leaf.lesions
            except:
                pass
        
        surface_inc = 0.
        surface_chlo = 0.
        surface_nec = 0.
        surface_spo = 0.
        
        # Temporary 
        surface_empty = 0.
        surface_total_nec = 0.
        
        if len(lesion_list)>0:
            for l in lesion_list:
                surface_inc += l.surface_inc
                surface_chlo += l.surface_chlo
                surface_nec += l.surface_nec
                surface_spo += l.surface_spo
                # Temporary 
                surface_empty += l.surface_empty
                surface_total_nec += l.necrotic_area

        self.surface_inc.append(surface_inc)
        self.surface_chlo.append(surface_chlo)
        self.surface_nec.append(surface_nec)
        self.surface_spo.append(surface_spo)
        self.surface_total_nec.append(surface_total_nec)
        self.ratio_inc.append(100. * surface_inc / total_area if total_area>0. else 0.)
        self.ratio_chlo.append(100. * surface_chlo / total_area if total_area>0. else 0.)
        self.ratio_nec.append(100. * surface_nec / total_area if total_area>0. else 0.)
        self.ratio_spo.append(100. * surface_spo / total_area if total_area>0. else 0.)
        self.ratio_total_nec.append(100. * (surface_total_nec) / total_area if total_area>0. else 0.)
        # Temporary
        self.surface_empty.append(surface_empty)
        self.ratio_empty.append(100. * surface_empty / total_area if total_area>0. else 0.)
            
    def count_dispersal_units(self, g):
        """ count DU of the leaf.
   
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
           
        Returns
        -------
            Number of dus on the leaf
        """
        dispersal_units = g.property('dispersal_units')
        return sum(1 for vid in self.leaf_elements for du in dispersal_units[vid] if du.is_active)
        
    def count_DU_on_healthy(self, g, nb_unwashed):
        """ Count DU of the leaf that are on healthy tissue.
        
        Same calculation as in BiotrophProbaModel. 
        Might not be the actual number of DUs on healthy tissue.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        nb_unwashed: int
            Number of DUs staying after washing
           
        Returns
        -------
            Number of dus on healthy tissue on the leaf
        """
        severity = self.compute_severity(g)/100
        return round(severity * nb_unwashed)

    # def compute_stock(self, g):
        # lesion_list = []
        # for id in self.ids:
            # leaf = g.node(id)
            # try:
                # lesion_list += leaf.lesions
            # except:
                # pass
        # if len(lesion_list)>0.:
            # self.stock_spores.append(sum([(l.stock_spores if l.stock_spores!=None else 0.) for l in lesion_list]))
        # else:
            # self.stock_spores.append(0.)
            
    # def compute_nb_spores_emitted(self, g):
        # lesion_list = []
        # for id in self.ids:
            # leaf = g.node(id)
            # try:
                # lesion_list += leaf.lesions
            # except:
                # pass
        # if len(lesion_list)>0.:
            # self.nb_spores_emitted.append(sum([(l.nb_spores_emitted if l.nb_spores_emitted!=None else 0.) for l in lesion_list]))
        # else:
            # self.nb_spores_emitted.append(0.)
            
    def update_audpc(self):
        try:
            ga = np.array(self.leaf_green_area)
        except:
            raise "Calculation of leaf green area with LeafInspector required"
        try:
            sev = np.array(self.severity)
        except:
            raise "Calculation of severity with LeafInspector required"
        if ga[-1]==0.:
            au = np.zeros(len(sev))
            au[ga>0]=sev[ga>0]
            tdelta = np.arange(len(sev)/24)
            self.audpc = trapz(au[::24], tdelta)
        else:
            self.audpc = 'audpc not available: leaf has not reached senescence'
        
#################################################################################
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
from scipy.integrate import trapz
from alinea.astk.plantgl_utils import get_height
from alinea.adel.newmtg import adel_ids
from collections import Iterable
try:
    import cPickle as pickle
except:
    import pickle

class AdelSeptoRecorder:
    def __init__(self, vids=None, adel_labels=None, group_dus=True):
        """ vids are ids of leaf sectors on the same blade. """
        self.vids = vids
        self.adel_labels = adel_labels
        self.group_dus = group_dus
        self.date_sequence = []
        self.degree_days = []
        # Initialize leaf properties to save
        self.leaf_area = []
        self.leaf_green_area = []
        self.leaf_senesced_length = []
        self.leaf_length = []
        # Initialize variables relative to DUs
        self.nb_dispersal_units = []
        self.nb_dus_on_green = []
        # Initialize variables relative to number of lesions
        self.nb_lesions = []
        self.nb_lesions_on_green = []
        # Initialize surfaces in state
        self.surface_inc = []
        self.surface_chlo = []
        self.surface_nec = []
        self.surface_spo = []
        self.surface_empty = []
        self.surface_dead = []
        # Initialize the index spotting the death of the leaf_element
        self.date_death = None

    def update_vids_with_labels(self, adel_ids):
        self.vids = [adel_ids[lb] for lb in self.adel_labels]
        
    def update_vids(self, vids=[]):
        self.vids = vids
    
    def record(self, g, date, degree_days=None):
        if self.date_death == None:
            vids = [id for id in self.vids if g.node(id).geometry is not None and g.node(id).area is not None]
        
            # Spot disappearing of leaf_element to stop recording
            if sum(self.leaf_area)>0 and len(vids)==0.:
                self.date_death = date
                return
        
            self.date_sequence.append(date)
            self.degree_days.append(degree_days)
        
            # Update leaf properties
            self.leaf_area.append(sum([g.node(id).area for id in vids]))
            self.leaf_green_area.append(sum([g.node(id).green_area for id in vids]))
            self.leaf_length.append(sum([g.node(id).length for id in vids]))
            self.leaf_senesced_length.append(sum([g.node(id).senesced_length for id in vids]))
                
            # Update properties of dispersal units and lesions
            nb_dus = 0
            nb_dus_on_green = 0
            nb_lesions = 0
            nb_lesions_on_green = 0
            surface_inc = 0.
            surface_chlo = 0.
            surface_nec = 0.
            surface_spo = 0.
            surface_empty = 0.
            surface_dead = 0.
            
            for id in vids:
                leaf = g.node(id)
                if 'dispersal_units' in leaf.properties():
                    leaf = g.node(id)
                    if self.group_dus:
                        for du in leaf.dispersal_units:
                            nb_dus += len(du.position)
                            nb_dus_on_green += len(filter(lambda x: x[0]>leaf.senesced_length, du.position))
                    else:
                        nb_dus += len(leaf.dispersal_units)
                        try:
                            nb_dus_on_green += len([du for du in leaf.dispersal_units if du.position[0]>leaf.senesced_length])
                        except:
                            nb_dus_on_green += len([du for du in leaf.dispersal_units if du.position[0][0]>leaf.senesced_length])
                
                if 'lesions' in leaf.properties():
                    if self.group_dus:
                        for les in leaf.lesions:
                            nb_lesions += len(les.position)
                            nb_lesions_on_green += len(filter(lambda x: x[0]>leaf.senesced_length, les.position))
                            surface_inc += les.surface_inc
                            surface_chlo += les.surface_chlo
                            surface_nec += les.surface_nec
                            surface_spo += les.surface_spo
                            surface_empty += les.surface_empty
                            surface_dead += les.surface_dead
                    else:
                        nb_lesions += len(leaf.lesions)
                        for les in leaf.lesions:
                            try:
                                if les.position[0]>leaf.senesced_length:
                                    nb_lesions_on_green += 1
                            except:
                                if les.position[0][0]>leaf.senesced_length:
                                    nb_lesions_on_green += 1
                              
                            surface_inc += les.surface_inc
                            surface_chlo += les.surface_chlo
                            surface_nec += les.surface_nec
                            surface_spo += les.surface_spo
                            surface_empty += les.surface_empty
                            surface_dead += les.surface_dead
            
            self.nb_dispersal_units.append(nb_dus)
            self.nb_dus_on_green.append(nb_dus_on_green)                
            self.nb_lesions.append(nb_lesions)
            self.nb_lesions_on_green.append(nb_lesions_on_green)
            self.surface_inc.append(surface_inc)
            self.surface_chlo.append(surface_chlo)
            self.surface_nec.append(surface_nec)
            self.surface_spo.append(surface_spo)
            self.surface_empty.append(surface_empty)
            self.surface_dead.append(surface_dead)
        
    def create_dataframe(self):
        exclude = ['vids', 'init_vids', 'adel_labels', 'group_dus']
        dico = {k:v for k,v in self.__dict__.iteritems() if k not in exclude}
        self.data = pandas.DataFrame(data=dico, index=self.date_sequence)
    
    def leaf_senesced_area(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        self.data['leaf_senesced_area'] = self.data['leaf_area'] - self.data['leaf_green_area']
    
    def leaf_disease_area(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        self.data['leaf_disease_area'] = self.data['surface_inc'] + self.data['surface_chlo'] + self.data['surface_nec'] + self.data['surface_spo'] + self.data['surface_empty'] + self.data['surface_dead']

    def leaf_lesion_area_on_green(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'leaf_disease_area' in self.data:
            self.leaf_disease_area()
        if not 'leaf_senesced_area' in self.data:
            self.leaf_senesced_area()
        self.data['leaf_lesion_area_on_green'] = [self.data['leaf_disease_area'][ind]*(1 - self.data['leaf_senesced_area'][ind]/self.data['leaf_area'][ind]) if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
    
    def leaf_necrotic_area_on_green(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'necrosis' in self.data:
            self.necrosis()
        if not 'leaf_senesced_area' in self.data:
            self.leaf_senesced_area()
        self.data['leaf_necrotic_area_on_green'] = [self.data['necrosis'][ind]*(1 - self.data['leaf_senesced_area'][ind]/self.data['leaf_area'][ind]) if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
    
    def leaf_necrotic_senescent(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'leaf_necrotic_area_on_green' in self.data:
            self.leaf_necrotic_area_on_green()
        if not 'leaf_senesced_area' in self.data:
            self.leaf_senesced_area()
        self.data['leaf_necrotic_senescent'] = self.data['leaf_senesced_area'] + self.data['leaf_necrotic_area_on_green']
        
    def necrotic_senescent_percentage(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'leaf_necrotic_senescent' in self.data:
            self.leaf_necrotic_senescent()
        self.data['leaf_necrotic_senescent'] = [self.data['leaf_necrotic_senescent'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        
    def leaf_healthy_area(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'leaf_lesion_area_on_green' in self.data:
            self.leaf_lesion_area_on_green()
        self.data['leaf_healthy_area'] = self.data['leaf_green_area'] - self.data['leaf_lesion_area_on_green']
    
    def leaf_unhealthy_area(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'leaf_healthy_area' in self.data:
            self.leaf_healthy_area()
        self.data['leaf_unhealthy_area'] = self.data['leaf_area'] - self.data['leaf_healthy_area']

    def surface_alive(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        self.data['surface_alive'] = self.data['surface_inc'] + self.data['surface_chlo'] + self.data['surface_nec'] + self.data['surface_spo'] + self.data['surface_empty']

    def ratios(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        self.data['ratio_inc'] = [self.data['surface_inc'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_chlo'] = [self.data['surface_chlo'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_nec'] = [self.data['surface_nec'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_spo'] = [self.data['surface_spo'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
        self.data['ratio_empty'] = [self.data['surface_empty'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]

    def severity(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'leaf_disease_area' in self.data:
            self.leaf_disease_area()
        self.data['severity'] = [self.data['leaf_disease_area'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
    
    def severity_on_green(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'severity' in self.data:
            self.severity()
        if not 'leaf_lesion_area_on_green' in self.data:
            self.leaf_lesion_area_on_green()
        self.data['severity_on_green'] = [self.data['leaf_lesion_area_on_green'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]

    def necrosis(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        self.data['necrosis'] = self.data['surface_nec'] + self.data['surface_spo'] + self.data['surface_empty']
    
    def necrosis_percentage(self):
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not 'necrosis' in self.data:
            self.necrosis()
        self.data['necrosis_percentage'] = [self.data['necrosis'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
                
    def get_audpc(self, variable='necrosis_percentage', max_ddays=None):
        """ Variable can be 'severity' or 'necrosis_percentage'. """
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not variable in self.data:
            exec('self.'+variable+'()')
        if self.date_death != None and len(self.leaf_green_area)>0:
            if self.data['leaf_green_area'][-1]==0.:
                ind = self.data.index[self.data['leaf_green_area']>0]
                if max_ddays != None:
                    df = self.data['degree_days'][ind]-self.data['degree_days'][ind][0]<max_ddays
                    ind = df.index
                data = self.data[variable][ind]
                self.audpc = sum([trapz(data[i:j], dx=(j-i).days*24 + (j-i).seconds/3600) for i,j in zip(data.index[:-1], data.index[1:])])
            else:
                self.audpc = 'audpc not available: leaf has not reached senescence'
        else:
            self.audpc = 'audpc not available: premature leaf death'
        return self.audpc
        
    def get_normalized_audpc(self, variable='necrosis_percentage'):
        """ Variable can be 'severity' or 'necrosis_percentage'. """
        if not 'data' in self.__dict__.keys():
            self.create_dataframe()
        if not variable in self.data:
            exec('self.'+variable+'()')
        if self.date_death != None and len(self.leaf_green_area)>0:
            if self.data['leaf_green_area'][-1]==0.:
                ind = self.data.index[self.data['leaf_green_area']>0]
                data = self.data[variable][ind]
                df_ref = pandas.DataFrame(data=numpy.ones(len(ind)), index=ind)
                ref_audpc = sum([trapz(df_ref[i:j].values.flatten(), dx=(j-i).days*24 + (j-i).seconds/3600) for i,j in zip(df_ref.index[:-1], df_ref.index[1:])])
                self.normalized_audpc = sum([trapz(data[i:j], dx=(j-i).days*24 + (j-i).seconds/3600) for i,j in zip(data.index[:-1], data.index[1:])])/ref_audpc if ref_audpc>0. else 0.
            else:
                self.normalized_audpc = 'audpc not available: leaf has not reached senescence'
        else:
            self.normalized_audpc = 'audpc not available: premature leaf death'
        return self.normalized_audpc
                
    def get_complete_dataframe(self):
        self.create_dataframe()
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

    def record_only_leaf_data(self, g, date, degree_days=None):
        vids = [id for id in self.vids if g.node(id).geometry is not None and g.node(id).area is not None]
    
        self.date_sequence.append(date)
        self.degree_days.append(degree_days)
    
        # Update leaf properties
        self.leaf_area.append(sum([g.node(id).area for id in vids]))
        self.leaf_green_area.append(sum([g.node(id).green_area for id in vids]))
        self.leaf_length.append(sum([g.node(id).length for id in vids]))
        self.leaf_senesced_length.append(sum([g.node(id).senesced_length for id in vids]))
        heights = [numpy.mean(get_height({vid:g.node(vid).geometry}).values()) for vid in vids]
        try:
            self.senesced_area.append(sum([g.node(id).senesced_area for id in vids]))
            self.leaf_green_length.append(sum([g.node(id).green_length for id in vids]))
            self.leaf_heights.append(heights if len(heights)>0. else [0.])
        except:
            self.senesced_area = [sum([g.node(id).senesced_area for id in vids])]
            self.leaf_green_length = [sum([g.node(id).green_length for id in vids])]
            self.leaf_heights = [heights if len(heights)>0. else [0.]]
                
    def create_dataframe_only_leaf_data(self):
        exclude = ['group_dus', 'vids', 'adel_labels', 'init_vids', 'surface_inc', 'surface_chlo', 
                   'surface_nec', 'surface_spo', 'surface_dead', 'surface_empty', 
                   'nb_lesions', 'nb_lesions_on_green', 'nb_dispersal_units', 
                   'nb_dus_on_green', 'senesced_length', 'leaf_heights']
        self.leaf_min_height = [min(heights) for heights in self.leaf_heights]
        self.leaf_max_height = [max(heights) for heights in self.leaf_heights]
        self.leaf_mean_height = [numpy.mean(heights) for heights in self.leaf_heights]
        dico = {k:v for k,v in self.__dict__.iteritems() if k not in exclude}
        self.data = pandas.DataFrame(data=dico, index=self.date_sequence)
        
    def inactivate(self, date=None):
        labels = self.adel_labels
        if sum(self.leaf_green_area)==0:
            self.__init__()
        self.adel_labels = labels
        self.date_death = date
        
def initiate_all_adel_septo_recorders(g, nsect=5):
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
            recorders['P%d' % ind_plant]['F%d' % ind_lf] = AdelSeptoRecorder(adel_labels = lf_labels)
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
    
def mean_by_leaf(recorder, variable='necrosis_percentage', skipna = False):
    ddays = recorder.values()[0].values()[0].degree_days
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
    
