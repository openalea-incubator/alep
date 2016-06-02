""" Test of the leaf inspector from alep.disease_outputs. """

# General imports
from random import sample, randint, seed

# Imports for wheat
from alinea.alep.wheat import initialize_stand, find_blade_id, find_leaf_ids

# Imports for lesions of septoria
from openalea.vpltk import plugin
from alinea.alep.disease_operation import generate_stock_lesions
from alinea.alep.disease_outputs import LeafInspector

# Useful functions ###################################################
def generate_lesions_with_surfaces(nb_lesions = 100, surface_inc = 0.1, surface_chlo = 0.2, 
                                   surface_nec = 0.3, surface_spo = 0.4):

    def necrotic_area_new(lesion):
        return lesion.surface_nec + lesion.surface_spo
    
    diseases=plugin.discover('alep.disease')
    septoria = diseases['septoria_exchanging_rings'].load()
    LesionKlass = septoria.lesion()
    LesionKlass.necrotic_area = necrotic_area_new
    lesion_stock = [LesionKlass(nb_spores=1) for i in range(nb_lesions)]
    for lesion in lesion_stock:
        lesion.surface_alive = surface_inc + surface_chlo + surface_nec + surface_spo
        lesion.surface_inc = surface_inc
        lesion.surface_chlo = surface_chlo
        lesion.surface_nec = surface_nec
        lesion.surface_spo = surface_spo
    return lesion_stock 
    
def distribute_lesions(g, leaf_ids, lesions):
    seed(7)   
    for le in leaf_ids:
        if len(lesions)>0:
            if le != leaf_ids[-1]:
                selection = sample(lesions, randint(0, len(lesions)-1))
                g.node(le).lesions = selection
                for les in selection:
                    lesions.remove(les)
            else:
                g.node(le).lesions= lesions
                lesions = []

def compute_ratios_new(inspector, g):
    total_area = 0.
    lesion_list = []
    for id in inspector.ids:
        leaf = g.node(id)
        total_area += leaf.area
        if leaf.lesions!=None:
            lesion_list += g.node(id).lesions
    
    surface_inc = 0.
    surface_chlo = 0.
    surface_nec = 0.
    surface_spo = 0.       
    if len(lesion_list)>0:
        for l in lesion_list:
            surface_inc += l.surface_inc
            surface_chlo += l.surface_chlo
            surface_nec += l.surface_nec
            surface_spo += l.surface_spo

    inspector.surface_inc.append(surface_inc)
    inspector.surface_chlo.append(surface_chlo)
    inspector.surface_nec.append(surface_nec)
    inspector.surface_spo.append(surface_spo)
    inspector.ratio_inc.append(100. * surface_inc / total_area if total_area>0. else 0.)
    inspector.ratio_chlo.append(100. * surface_chlo / total_area if total_area>0. else 0.)
    inspector.ratio_nec.append(100. * surface_nec / total_area if total_area>0. else 0.)
    inspector.ratio_spo.append(100. * surface_spo / total_area if total_area>0. else 0.)

# Test ###############################################################
# Initialize a wheat canopy
g, wheat, domain_area = initialize_stand(age=1000., length=0.1, 
                                        width=0.2, sowing_density=150,
                                        plant_density=150, inter_row=0.12)

# Generate a pool of lesions with given surface to be distributed on target leaves
nb_lesions = 50; surface_inc = 0.1; surface_chlo = 0.2; surface_nec = 0.3; surface_spo = 0.4
lesions = generate_lesions_with_surfaces(nb_lesions, surface_inc, surface_chlo, surface_nec, surface_spo)

# Select target blades (top 3)
blade_1 = find_blade_id(g, leaf_rank = 1)
blade_2 = find_blade_id(g, leaf_rank = 2)
blade_3 = find_blade_id(g, leaf_rank = 3)

# Distribute the lesions
target_leaves = find_leaf_ids(g, blade_1) + find_leaf_ids(g, blade_2) + find_leaf_ids(g, blade_3)
distribute_lesions(g, target_leaves, lesions)

# Call leaf inspectors
LeafInspector.compute_ratios = compute_ratios_new
inspector_1 = LeafInspector(g, blade_id = blade_1)
inspector_2 = LeafInspector(g, blade_id = blade_2)
inspector_3 = LeafInspector(g, blade_id = blade_3)

# Compute disease coverage
inspector_1.update_variables(g)
inspector_2.update_variables(g)
inspector_3.update_variables(g)

# Test
assert (round(inspector_1.surface_inc[0] + inspector_2.surface_inc[0] +
        inspector_3.surface_inc[0], 6) == round(nb_lesions*surface_inc, 6))
assert (round(inspector_1.surface_chlo[0] + inspector_2.surface_chlo[0] +
        inspector_3.surface_chlo[0], 6) == round(nb_lesions*surface_chlo, 6))
assert (round(inspector_1.surface_nec[0] + inspector_2.surface_nec[0] +
        inspector_3.surface_nec[0], 6) == round(nb_lesions*surface_nec, 6))
assert (round(inspector_1.surface_spo[0] + inspector_2.surface_spo[0] +
        inspector_3.surface_spo[0], 6) == round(nb_lesions*surface_spo, 6))