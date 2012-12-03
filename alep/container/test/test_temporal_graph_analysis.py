# -*- python -*-
#
#       OpenAlea.Container
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s):  Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                        Frederic Boudon <frederic.boudon@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################
from openalea.container.temporal_graph_analysis import laplacian, mean_abs_dev, temporal_change, relative_temporal_change
from temporal_property_graph_input import create_TPG_VertexProperty


def test_vid_type():
    """ Test the relabeling of the graph vertices """
    g = create_TPG_VertexProperty()
    assert g.vertex_property('old_label') == {0: 0, 1: 1, 2: 0, 3: 1, 4: 2, 5: 3, 6: 4, 7: 0, 8: 1, 9: 2, 10: 3, 11: 4, 12: 5, 13: 6, 14: 7, 15: 8, 16: 9, 17: 10, 18: 11, 19: 12}

def test_laplacian():
    """
    Test the laplacian function in 'temporal_graph_analysis':
        laplacian(i) = [vertex_property(i) - sum(vertex_property(neighbors(i)))] / nb_neighbors
    """
    g = create_TPG_VertexProperty()
    assert laplacian(g,'property', range(7)) == { 0 : -1.,  1 : 1., 2 : -1., 3 : (3-11/3.), 4 : 0., 5 : (5-13/3.), 6 : 1. }
    
def test_mean_abs_dev():
    """
    Test the mean_abs_dev function in 'temporal_graph_analysis':
        mean_abs_dev(i) = [vertex_property(i) - sum( |vertex_property(neighbors(i))| )] / nb_neighbors
    """
    g = create_TPG_VertexProperty()
    assert mean_abs_dev(g,'property', range(7)) == { 0 : 1,  1 : 1, 2 : 1., 3 : 4/3., 4 : 1, 5 : 4/3., 6 : 1 }
    
def test_temporal_change():
    g = create_TPG_VertexProperty()
    res_n = { 0 : 9,  1 : 10, 2 : 13, 3 : 6, 4 : 42, 5 : 24, 6 : 64 }
    res_n1 = {2: 9.0, 3: 9.0, 4: 9.0, 5: 10.0, 6: 10.0, 7: 13.0, 8: 13.0, 9: 6.0, 10: 42.0, 11: 42.0, 12: 42.0, 13: 42.0, 14: 24.0, 15: 24.0, 16: 64.0, 17: 64.0, 18: 64.0, 19: 64.0}
    assert temporal_change(g,'property', range(7)) == res_n
    # If no descendants for a vid, no result returned.
    assert temporal_change(g,'property', range(8)) == res_n
    # Check relabelling for ids @ t_n+1. If no ancestors for a vid, no result returned.
    assert temporal_change(g,'property', range(7), labels_at_t_n = False) == res_n1
    # If no descendants for a vid, no result returned.
    assert temporal_change(g,'property', range(8), check_full_lineage=False ) == res_n

def test_relative_temporal_change():
    g = create_TPG_VertexProperty()
    res = { 1 : 10, 2 : 13/2., 3 : 2, 4 : 10.5, 5 : 4.8, 6 : 64/6. }
    assert relative_temporal_change(g,'property', range(1,7)) == res
    # If no descendants for a vid, no result returned.
    assert relative_temporal_change(g,'property', range(1,8)) == res
    # If no descendants for a vid, no result returned.
    assert relative_temporal_change(g,'property', range(1,8), check_full_lineage=False) == res
