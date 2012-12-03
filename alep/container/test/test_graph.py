from nose import with_setup
from openalea.container import Graph

g=Graph()

def setup_func () :
    for i in xrange(10) :
        g.add_vertex(i)
    for i in xrange(9) :
        g.add_edge(i,i+1,i)

def teardown_func () :
    g.clear()

# ##########################################################
#
# Graph concept
#
# ##########################################################
@with_setup(setup_func,teardown_func)
def test_source () :
    for i in xrange(9) :
        assert g.source(i)==i

@with_setup(setup_func,teardown_func)
def test_target () :
    for i in xrange(9) :
        assert g.target(i)==(i+1)

@with_setup(setup_func,teardown_func)
def test_has_vertex () :
    for i in xrange(10) :
        assert g.has_vertex(i)

@with_setup(setup_func,teardown_func)
def test_has_edge () :
    for i in xrange(9) :
        assert g.has_edge(i)

@with_setup(setup_func,teardown_func)
def test_is_valid () :
    assert g.is_valid()

# ##########################################################
#
# Vertex List Graph Concept
#
# ##########################################################
@with_setup(setup_func,teardown_func)
def test_vertices () :
    assert list(g.vertices())==range(10)

@with_setup(setup_func,teardown_func)
def test_nb_vertices () :
    assert g.nb_vertices()==10

@with_setup(setup_func,teardown_func)
def test_in_neighbors () :
    for i in xrange(9) :
        assert list(g.in_neighbors(i+1))==[i]

@with_setup(setup_func,teardown_func)
def test_out_neighbors () :
    for i in xrange(9) :
        assert list(g.out_neighbors(i))==[i+1]

@with_setup(setup_func,teardown_func)
def test_neighbors () :
    for i in xrange(8) :
        neis=list(g.neighbors(i+1))
        assert i in neis
        assert i+2 in neis

@with_setup(setup_func,teardown_func)
def test_nb_in_neighbors () :
    for i in xrange(9) :
        assert g.nb_in_neighbors(i+1)==1

@with_setup(setup_func,teardown_func)
def test_nb_out_neighbors () :
    for i in xrange(9) :
        assert g.nb_out_neighbors(i)==1

@with_setup(setup_func,teardown_func)
def test_nb_neighbors () :
    for i in xrange(8) :
        assert g.nb_neighbors(i+1)==2

@with_setup(setup_func,teardown_func)
def test_edge () :
    assert g.edge(0,1) == 0
    assert g.edge(0,2) == None

# ##########################################################
#
# Edge List Graph Concept
#
# ##########################################################
@with_setup(setup_func,teardown_func)
def test_edges () :
    assert list(g.edges())==range(9)

@with_setup(setup_func,teardown_func)
def test_nb_edges () :
    assert g.nb_edges()==9

@with_setup(setup_func,teardown_func)
def test_in_edges () :
    for i in xrange(9) :
        assert list(g.in_edges(i+1))==[i]

@with_setup(setup_func,teardown_func)
def test_out_edges () :
    for i in xrange(9) :
        assert list(g.out_edges(i))==[i]

@with_setup(setup_func,teardown_func)
def test_vertex_edges () :
    for i in xrange(8) :
        neis=list(g.edges(i+1))
        assert i in neis
        assert i+1 in neis

@with_setup(setup_func,teardown_func)
def test_nb_in_edges () :
    for i in xrange(9) :
        assert g.nb_in_edges(i+1)==1

@with_setup(setup_func,teardown_func)
def test_nb_out_edges () :
    for i in xrange(9) :
        assert g.nb_out_edges(i)==1

@with_setup(setup_func,teardown_func)
def test_nb_edges () :
    for i in xrange(8) :
        assert g.nb_edges(i+1)==2

# ##########################################################
#
# Mutable Vertex Graph concept
#
# ##########################################################
@with_setup(setup_func,teardown_func)
def test_add_vertex () :
    assert g.add_vertex(100)==100
    vid=g.add_vertex()
    assert g.has_vertex(vid)

@with_setup(setup_func,teardown_func)
def test_remove_vertex () :
    g.remove_vertex(5)
    assert not g.has_vertex(5)
    assert not g.has_edge(4)
    assert not g.has_edge(5)
    assert 5 not in list(g.neighbors(6))
    assert 5 not in list(g.neighbors(4))

@with_setup(setup_func,teardown_func)
def test_clear () :
    g.clear()
    assert g.nb_vertices()==0
    assert g.nb_edges()==0

# ##########################################################
#
# Mutable Edge Graph concept
#
# ##########################################################
@with_setup(setup_func,teardown_func)
def test_add_edge () :
    assert g.add_edge(0,9,100)==100
    eid=g.add_edge(2,1)
    assert eid in list(g.in_edges(1))
    assert eid in list(g.out_edges(2))

@with_setup(setup_func,teardown_func)
def test_remove_edge () :
    g.remove_edge(4)
    assert not g.has_edge(4)
    assert 4 not in list(g.neighbors(5))
    assert 5 not in list(g.neighbors(4))

@with_setup(setup_func,teardown_func)
def test_clear_edges () :
    g.clear_edges()
    assert g.nb_vertices()==10
    assert g.nb_edges()==0

# ##########################################################
#
# Extend Graph concept
#
# ##########################################################
@with_setup(setup_func,teardown_func)
def test_extend () :
    trans_vid,trans_eid=g.extend(g)
    assert len(trans_vid)==10
    assert len(trans_eid)==9
