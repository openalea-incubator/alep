""" Run a simulation to compare strategies of dispersal. """

# Imports for visualization
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer

# Imports for wheat
from alinea.alep.wheat import initialize_stand
from alinea.alep.architecture import set_properties,set_property_on_each_id, get_leaves

# Imports for disease
from alinea.alep.fungal_objects import DispersalUnit, Lesion
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.alep.dispersal import PowderyMildewWindDispersal
from alinea.alep.disease_outputs import count_dispersal_units_by_leaf
from alinea.alep.alep_color import alep_colormap, green_yellow_red
from alinea.alep.protocol import disperse

from alinea.alep.disease_outputs import plot_dispersal_units

# Useful functions ############################################################
def update_plot(g):
    # Compute severity by leaf
    nb_dus_by_leaf = count_dispersal_units_by_leaf(g, label = 'LeafElement')
    set_property_on_each_id(g, 'nb_dus', nb_dus_by_leaf, label = 'LeafElement')

    # Visualization
    g = alep_colormap(g, 'nb_dus', cmap=green_yellow_red(levels=100),
                      lognorm=False, zero_to_one=False, vmax=100)

    leaves = get_leaves(g, label='LeafElement')
    pos_sen = g.property('position_senescence')
    for leaf in leaves:
        if pos_sen[leaf]==0.:
            g.node(leaf).color = (157, 72, 7)

    scene = plot3d(g)
    Viewer.display(scene)
    return scene

# Dispersal ################################################################### 
# Initialize a wheat canopy
g, wheat, domain_area = initialize_stand(age=700., length=1,
                                        width=1, sowing_density=150,
                                        plant_density=150, inter_row=0.12)
g2, wheat2, domain_area2 = initialize_stand(age=700., length=1,
                                        width=1, sowing_density=150,
                                        plant_density=150, inter_row=0.12)
                                       
# Create a undetermined lesion emmitting a stock of dispersal units
class DummyLesion(Lesion):
    """ Undetermined lesion with only a method emission."""
    def __init__(self, nb_spores=None, position=None):
        super(DummyLesion, self).__init__(nb_spores=nb_spores, position=position)
        class params():
            def __init__(self, name="dummy"):
                self.name = name
        self.fungus = params()
        # Stock of spores needs to be positive
        # TODO : change this condition to a method 'can_emit'
        self.stock_spores = 1
        
    def emission(self, leaf=None):
        return [DispersalUnit(nb_spores=self.nb_spores, status='emitted') for i in range(10000)]

        
# leaf_id = 66798
leaf_id = 74
leaf = g.node(leaf_id)
leaf2 = g2.node(leaf_id)
# leaf.color = (0,0,180)
# leaf2.color = (0,0,180)
leaf.lesions = [DummyLesion()]

dispersor = Septo3DSplash(reference_surface=domain_area)

# Define the wind direction
wind_direction = (1,0,0)
set_properties(g,label = 'LeafElement', wind_direction=wind_direction)
dispersor2 = PowderyMildewWindDispersal()
disperse(g, dispersor, "dummy", label='LeafElement')
disperse(g2, dispersor2, "dummy", label='LeafElement')

scene = plot3d(g)
Viewer.display(scene)

plot_dispersal_units(g)
plot_dispersal_units(g2)

