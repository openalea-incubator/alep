/* -*-c++-*- 
 *------------------------------------------------------------------------------
 *                                                                              
 *        openalea.container: utils package                                     
 *                                                                              
 *        Copyright 2006 INRIA - CIRAD - INRA                      
 *                                                                              
 *        File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>         
 *                                                                              
 *        Distributed under the Cecill-C License.                               
 *        See accompanying file LICENSE.txt or copy at                          
 *            http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html       
 *                                                                              
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr                    
 *       
 *                                                                       
 *-----------------------------------------------------------------------------*/
#include <vector>

#include "container/relation.h"
#include "container/graph.h"
#include "container/graph_iterators.h"
#include "container/int_iterator.h"
using namespace container;

#include <boost/smart_ptr.hpp>
#include <boost/python.hpp>
#include "export_iterator.h"
#include "int_range.h"
using namespace boost::python;

typedef PyCustomRange<Graph::in_edge_iterator> graph_local_edge_range;
typedef PyCustomRange<Graph::edge_iterator> graph_edge_range;
typedef PyCustomRange<Graph::vertex_edge_iterator> graph_vertex_edge_range;

graph_local_edge_range export_in_edges (Graph& g, int vid) {
	return graph_local_edge_range(g.in_edges_begin(vid),g.in_edges_end(vid));
}
graph_local_edge_range export_out_edges (Graph& g, int vid) {
	return graph_local_edge_range(g.out_edges_begin(vid),g.out_edges_end(vid));
}
graph_edge_range export_graph_edges (Graph& g) {
	return graph_edge_range(g.edges_begin(),g.edges_end());
}

typedef IntCloneIterator<Graph::edge_iterator> int_edge_iterator;
typedef IntCloneIterator<Graph::vertex_edge_iterator> int_vertex_edge_iterator;


/*PyIntRange export_graph_vertex_edges (Graph& graph, PyObject* arg) {
	if(arg==Py_None) {
		return PyIntRange(int_edge_iterator(&graph.edges_begin()),
				  int_edge_iterator(&graph.edges_end()));
	}
	else {
		int vid=extract<int>(arg);
		return PyIntRange(int_vertex_edge_iterator(&graph.edges_begin(vid)),
				  int_vertex_edge_iterator(&graph.edges_end(vid)));
	}
}*/

graph_vertex_edge_range export_graph_vertex_edges (Graph& g, int vid) {
	return graph_vertex_edge_range(g.edges_begin(vid),g.edges_end(vid));
}

typedef PyCustomRange<Graph::vertex_iterator> graph_vertex_range;
typedef PyCustomRange<Graph::in_neighbor_iterator> graph_in_neighbor_range;
typedef PyCustomRange<Graph::out_neighbor_iterator> graph_out_neighbor_range;
typedef PyCustomRange<Graph::neighbor_iterator> graph_neighbor_range;

graph_vertex_range export_iter_vertices (Graph& g) {
	return graph_vertex_range(g.vertices_begin(),g.vertices_end());
}
graph_in_neighbor_range export_in_neighbors (Graph& g, int vid) {
	return graph_in_neighbor_range(g.in_neighbors_begin(vid),g.in_neighbors_end(vid));
}
graph_out_neighbor_range export_out_neighbors (Graph& g, int vid) {
	return graph_out_neighbor_range(g.out_neighbors_begin(vid),g.out_neighbors_end(vid));
}
graph_neighbor_range export_neighbors (Graph& g, int vid) {
	return graph_neighbor_range(g.neighbors_begin(vid),g.neighbors_end(vid));
}

int export_graph_nb_edges (Graph& graph, PyObject* arg) {
	if(arg==Py_None) {
		return graph.nb_edges();
	}
	else {
		return graph.nb_edges(extract<int>(arg));
	}
}

int export_graph_add_vertex (Graph& graph, PyObject* arg) {
	if(arg==Py_None) {
		return graph.add_vertex();
	}
	else {
		return graph.add_vertex(extract<int>(arg));
	}
}

int export_graph_add_edge (Graph& graph, int sid, int tid, PyObject* arg) {
	if(arg==Py_None) {
		return graph.add_edge(sid,tid);
	}
	else {
		return graph.add_edge(sid,tid,extract<int>(arg));
	}
}

class TotoA {
public:
	virtual void test () const {std::cout << "A" << std::endl;}
	void test_neighbors_iterators () {
		Graph g;
		for(int i=0;i<1000;++i) {
			g.add_vertex(i);
		}
		for(int i=0;i<1000;++i) {
			g.add_edge(i,(i+1)%1000);
		}
		std::vector<int> vid_list;
		for(int i=0;i<1000;++i) {
			Graph::neighbor_iterator it=g.neighbors_begin(i);
			for(it;it!=g.neighbors_end(i);++it) {
				vid_list.push_back(*it);
			}
		}
		std::cout << "end neighbors " << vid_list.size() << std::endl;
	}
	void test_edges_iterators () {
		Graph g;
		for(int i=0;i<1000;++i) {
			g.add_vertex(i);
		}
		for(int i=0;i<1000;++i) {
			g.add_edge(i,(i+1)%1000);
		}
		std::vector<int> eid_list;
		for(int i=0;i<1000;++i) {
			Graph::vertex_edge_iterator it=g.edges_begin(i);
			for(it;it!=g.edges_end(i);++it) {
				eid_list.push_back(*it);
			}
		}
		std::cout << "end edgers " << eid_list.size() << std::endl;
	}
};

void export_graph () {
	export_custom_range<Graph::vertex_edge_iterator>("_PyGraphVertexEdgeRange");
	export_custom_range<Graph::vertex_iterator>("_PyGraphVertexRange");
	export_custom_range<Graph::neighbor_iterator>("_PyGraphNeighborRange");

	class_<TotoA>("TotoA")
		.def("test",&TotoA::test)
		.def("test_neighbor_iterator",&TotoA::test_neighbors_iterators)
		.def("test_edge_iterator",&TotoA::test_edges_iterators);


	class_<Graph, bases<Relation> >("Graph", "graph")
		//graph concept
		.def("has_vertex",&Graph::has_vertex,"test wether a vertex is in the graph")
		.def("__contains__",&Graph::has_vertex,"test wether a vertex is in the graph")
		.def("has_edge",&Graph::has_edge,"test wether an eid is in the graph")
		//edge list
		.def("in_edges",&export_in_edges,"iterator on in edges connected to a vertex")
		.def("nb_in_edges",&Graph::nb_in_edges,"number of in edges connected to a vertex")
		.def("out_edges",&export_out_edges,"iterator on out edges connected to a vertex")
		.def("nb_out_edges",&Graph::nb_out_edges,"number of out edges connected to a vertex")
		.def("edges",&export_graph_edges,"iterator on all edges")
		.def("nb_edges",(int (Graph::*) () const)& Graph::nb_edges,"total number of edges")
		.def("edges",&export_graph_vertex_edges,"iterator on all edges connected to a vertex")
		.def("nb_edges",&export_graph_nb_edges,"number of out edges connected to a vertex")
		//vertex list
		.def("vertices",&export_iter_vertices,"iterator on all vertices")
		.def("nb_vertices",&Graph::nb_vertices,"number of vertices in the graph")
		.def("__len__",&Graph::nb_vertices,"number of vertices in the graph")
		.def("in_neighbors",&export_in_neighbors,"iterator on neighbors of a vertex")
		.def("nb_in_neighbors",&Graph::nb_in_neighbors,"number of neighbors of this vertex")
		.def("out_neighbors",&export_out_neighbors,"iterator on neighbors of a vertex")
		.def("nb_out_neighbors",&Graph::nb_out_neighbors,"number of neighbors of this vertex")
		.def("neighbors",&export_neighbors,"iterator on neighbors of a vertex")
		.def("nb_neighbors",&Graph::nb_neighbors,"number of neighbors of this vertex")
		//mutable
		.def("add_vertex",(int (Graph::*) ())& Graph::add_vertex,"add a new vertex in the graph")
		//.def("add_vertex",(int (Graph::*) (int))& Graph::add_vertex,"add a new vertex in the graph")
		.def("add_vertex",&export_graph_add_vertex,"add a new vertex in the graph")
		.def("remove_vertex",&Graph::remove_vertex,"remove a vertex from the graph")
		.def("add_edge",(int (Graph::*) (int,int))& Graph::add_edge,"add a new edge between two vertices")
		.def("add_link",(int (Graph::*) (int,int))& Graph::add_link,"add a new edge between two vertices")
		//.def("add_edge",(int (Graph::*) (int,int,int))& Graph::add_edge,"add a new edge between two vertices")
		.def("add_edge",&export_graph_add_edge,"add a new edge between two vertices")
		//.def("add_link",(int (Graph::*) (int,int,int))& Graph::add_link,"add a new edge between two vertices")
		.def("add_link",&export_graph_add_edge,"add a new edge between two vertices")
		.def("remove_edge",&Graph::remove_edge,"remove an edge, do not remove corresponding vertices")
		.def("remove_link",&Graph::remove_link,"remove an edge, do not remove corresponding vertices")
		.def("clear_edges",&Graph::clear_edges,"remove all links")
		.def("clear",&Graph::clear,"remove all vertices from the graph")
		//debug
		.def("state",&Graph::state,"debug function");

}
