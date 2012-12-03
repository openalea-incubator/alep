
__license__ = "Cecill-C"
__revision__ = " $Id: test_property_graph.py 7865 2010-02-08 18:04:39Z cokelaer $ "

# Test node module
from openalea.container import PropertyGraph
from openalea.container import TemporalPropertyGraph
from temporal_property_graph_input import create_TemporalGraph

def test_temporalPropertyGraph():
    """create a graph"""
    g = create_TemporalGraph()

    # ~ begin tests
    
    # ~ Neighbors
    assert g.in_neighbors(13)==set([4, 9, 10, 11, 12])
    assert g.in_neighbors(13, 't')==set([4])
    assert g.in_neighbors(13, 's')==set([9, 10, 11, 12])
    assert g.out_neighbors(4)==set([5, 10, 11, 12, 13])
    assert g.out_neighbors(4, 't')==set([10, 11, 12, 13])
    assert g.out_neighbors(4, 's')==set([5])
    
    # ~ Edges
    assert g.out_edges(5)==set([4, 36, 37])
    assert g.out_edges(5, 't')==set([36, 37])
    assert g.out_edges(5, 's')==set([4])
    assert g.in_edges(5)==set([3, 5, 9])
    assert g.in_edges(5, 't')==set([9])
    assert g.in_edges(5, 's')==set([3, 5])
    
    # ~ Sibling
    assert g.sibling(1)==None
    assert g.sibling(16)==set([17, 18, 19])
    
    # ~ Neighborhood
    assert g.neighborhood(2, 2)==set([0, 1, 2, 3, 4, 5, 7, 8, 9, 10])
    assert g.neighborhood(2, 2, 't')==set([0, 2, 3, 4, 7, 8])
    assert g.neighborhood(2, 2, 's')==set([2, 3, 4, 5])

    # ~ Topologic Distance
    assert g.topological_distance(3)=={0: 1,
                                       1: 2,
                                       2: 1,
                                       3: 0,
                                       4: 1,
                                       5: 1,
                                       6: 2,
                                       7: 2,
                                       8: 2,
                                       9: 1,
                                       10: 2,
                                       11: 2,
                                       12: 2,
                                       13: 2,
                                       14: 2,
                                       15: 2,
                                       16: 3,
                                       17: 3,
                                       18: 3,
                                       19: 3}
                          
    assert g.topological_distance(3, edge_type='t')=={0: 1,
                                                      1: float('Inf'),
                                                      2: 2,
                                                      3: 0,
                                                      4: 2,
                                                      5: float('Inf'),
                                                      6: float('Inf'),
                                                      7: 3,
                                                      8: 3,
                                                      9: 1,
                                                      10: 3,
                                                      11: 3,
                                                      12: 3,
                                                      13: 3,
                                                      14: float('Inf'),
                                                      15: float('Inf'),
                                                      16: float('Inf'),
                                                      17: float('Inf'),
                                                      18: float('Inf'),
                                                      19: float('Inf')}
    
    assert g.topological_distance(3, edge_type='s')=={0: float('Inf'),
                                                      1: float('Inf'),
                                                      2: 1,
                                                      3: 0,
                                                      4: 1,
                                                      5: 1,
                                                      6: 2,
                                                      7: float('Inf'),
                                                      8: float('Inf'),
                                                      9: float('Inf'),
                                                      10: float('Inf'),
                                                      11: float('Inf'),
                                                      12: float('Inf'),
                                                      13: float('Inf'),
                                                      14: float('Inf'),
                                                      15: float('Inf'),
                                                      16: float('Inf'),
                                                      17: float('Inf'),
                                                      18: float('Inf'),
                                                      19: float('Inf')}
    
    assert g.topological_distance(3, edge_dist=lambda x,y: 2)=={0: 2,
                                                                1: 4,
                                                                2: 2,
                                                                3: 0,
                                                                4: 2,
                                                                5: 2,
                                                                6: 4,
                                                                7: 4,
                                                                8: 4,
                                                                9: 2,
                                                                10: 4,
                                                                11: 4,
                                                                12: 4,
                                                                13: 4,
                                                                14: 4,
                                                                15: 4,
                                                                16: 6,
                                                                17: 6,
                                                                18: 6,
                                                                19: 6}

def func_regional_bin(graph, vid):
    return (len(graph.neighborhood(vid))-1)==4

def test_regionalisation():
    g=create_TemporalGraph()
    g.add_region(func_regional_bin, "dist4")
    assert g._graph_property=={'dist4' : [0, 2, 9, 11, 15]}
    assert g._vertex_property["regions"] == {0 : ['dist4'], 
                                             2 : ['dist4'], 
                                             9 : ['dist4'], 
                                             11: ['dist4'], 
                                             15: ['dist4']}
    
    g.add_vertex_to_region(1, "dist4")
    assert 1 in g._graph_property["dist4"]
    assert "dist4" in g._vertex_property["regions"][1]

    g.remove_vertex_from_region(2, "dist4")
    assert not 2 in g._graph_property["dist4"]
    assert g._vertex_property["regions"].get(2) == None
    
    g.remove_vertex_from_region([0, 1], "dist4")
    assert not 0 in g._graph_property["dist4"]
    assert g._vertex_property["regions"].get(0) == None
    assert not 1 in g._graph_property["dist4"]
    assert g._vertex_property["regions"].get(1) == None

    assert not g.is_connected_region("dist4")
    
    g.add_vertex_to_region(3, "dist4")
    g.add_vertex_to_region([4, 5], "dist4")
    
    assert g.is_connected_region("dist4")
    assert not g.is_connected_region("dist4", "s")

    g.remove_region("dist4")
    assert g._graph_property.get("dist4") == None
    assert g._vertex_property["regions"] == {}
