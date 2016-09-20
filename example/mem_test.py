from alinea.adel.data_samples import adel_two_metamers_stand
from memory_profiler import profile
from alinea.alep.septoria_age_physio import *
from alinea.alep.simulation_tools.variable_septoria import *

def create_g(leaf_sectors = 1, density = 350, interleaf = 10, 
            leaf_length = 20, leaf_width = 1, Einc = 0):
    g, domain_area, domain, convunit = adel_two_metamers_stand(leaf_sectors = leaf_sectors,
                                                                density = density, 
                                                                interleaf = interleaf, 
                                                                leaf_length = leaf_length,
                                                                leaf_width = leaf_length,
                                                                Einc = Einc)
    return g, domain_area, domain, convunit

def del_dus(g):
    dus = g.property('dispersal_units')
    for vid, du_list in dus.iteritems():
        for du in du_list:
            du.disable()
        dus[vid] = [du for du in du_list if du.is_active]    
    
@profile
def my_func():
    g, domain_area, domain, convunit = create_g()
    septo = SeptoriaFungus()
    g.node(12).dispersal_units = [septo.dispersal_unit(mutable = False) for i in range(100000)]
    fungus = variable_septoria(distri_chlorosis = {'mu':150, 'sigma':30})
    g.node(12).dispersal_units = [fungus.dispersal_unit(mutable = True) for i in range(100000)]
    del_dus(g)
    
if __name__ == '__main__':
    my_func()