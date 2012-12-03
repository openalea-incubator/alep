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

#ifndef __CONTAINER_GRAPH_H__
#define __CONTAINER_GRAPH_H__

#include <utility>
#include <map>
#include "container/id_map.h"
#include "container/relation.h"
#include "container/relation_iterators.h"

namespace container {
	class EdgeIterator;
	class VertexIterator;

	class Graph : public Relation {
	private:
		typedef std::pair<std::set<int>,std::set<int> > vpair;
		typedef IdMap<vpair>::iterator v_iterator;
	public:
		typedef IdMap<vpair>::key_iterator vertex_iterator;//iter of int
		typedef LinkTargetIterator out_neighbor_iterator;//iter of int
		typedef LinkSourceIterator in_neighbor_iterator;//iter of int
		typedef VertexIterator neighbor_iterator;//iter of int

		typedef Relation::iterator edge_iterator;//iter of int
		typedef std::set<int>::iterator out_edge_iterator;//iter of int
		typedef std::set<int>::iterator in_edge_iterator;//iter of int
		typedef EdgeIterator vertex_edge_iterator;//iter of int

	public:
		struct InvalidVertex : public IdMap<vpair>::InvalidId {
		//small exception raised when attempting to access an invalid vertex
			InvalidVertex (int l) : IdMap<vpair>::InvalidId(l) {};
			virtual ~InvalidVertex() throw () {};
			InvalidVertex* copy () { return new InvalidVertex(id); }
			virtual void rethrow () { throw InvalidVertex(id); }
		};
		typedef Relation::InvalidLink InvalidEdge;
	private:
		IdMap<vpair> vertex_links;//map of vertex_id:(set(in_links),set(out_links))
	private:
		void test_existing_vertex (int vid) const;
		//throw InvalidVertex if vid not in the graph

		void add_edge_to_vertex (int sid, int tid, int eid);
		//internal function to add eid to vertex_links

	public:
		Graph() : Relation() {}
		~Graph() {}

		bool has_vertex (int vid) const {
		//return true if the vertex is inside the graph
			return vertex_links.find(vid)!=vertex_links.end();
		}

		bool has_edge (int eid) const {
		//return true if the edge is in the graph
			return has_link(eid);
		}

		//
		//
		// edge concept
		//
		//
		in_edge_iterator in_edges_begin (int vid) {
		//iterator on all edges going to a vertex
			test_existing_vertex(vid);
			return vertex_links.find(vid)->second.first.begin();
		}

		in_edge_iterator in_edges_end (int vid) {
		//iterator on all edges going to a vertex
			test_existing_vertex(vid);
			return vertex_links.find(vid)->second.first.end();
		}

		int nb_in_edges (int vid) const {
		//number of edges going to a vertex
			test_existing_vertex(vid);
			return vertex_links.find(vid)->second.first.size();
		}

		out_edge_iterator out_edges_begin (int vid) {
		//iterator on edges going outide a vertex
			test_existing_vertex(vid);
			return vertex_links.find(vid)->second.second.begin();
		}

		out_edge_iterator out_edges_end (int vid) {
		//iterator on edges going outside a vertex
			test_existing_vertex(vid);
			return vertex_links.find(vid)->second.second.end();
		}

		int nb_out_edges (int vid) const {
		//number of edges going outside a vertex
			test_existing_vertex(vid);
			return vertex_links.find(vid)->second.second.size();
		}

		edge_iterator edges_begin () {
		//iterator on all edges in the graph
			return begin();
		}

		edge_iterator edges_end () {
		//iterator on all edges in the graph
			return end();
		}

		int nb_edges () const {
		//number of edges in the graph
			return size();
		}

		vertex_edge_iterator edges_begin (int vid);
		//iterator on edges around a vertex

		vertex_edge_iterator edges_end (int vid);
		//iterator on edges around a vertex

		int nb_edges (int vid) const {
		//number of edges around a vertex
			test_existing_vertex(vid);
			return nb_in_edges(vid)+nb_out_edges(vid);
		}

		//
		//
		// vertex concept
		//
		//
		vertex_iterator vertices_begin () {
		//iterator on all vertices in the graph
			return vertex_links.key_begin();
		}

		vertex_iterator vertices_end () {
		//iterator on all vertices in the graph
			return vertex_links.key_end();
		}

		int nb_vertices () const {
		//total number of vertices in the graph
			return vertex_links.size();
		}

		in_neighbor_iterator in_neighbors_begin (int vid) {
		//iterator on all vertices connected to a vertex
			return LinkSourceIterator(*this,in_edges_begin(vid));
		}

		in_neighbor_iterator in_neighbors_end (int vid) {
		//iterator on all vertices connected to a vertex
			return LinkSourceIterator(*this,in_edges_end(vid));
		}

		int nb_in_neighbors (int vid) const {
		//number of vertices connected to a vertex
			return nb_in_edges(vid);
		}

		out_neighbor_iterator out_neighbors_begin (int vid) {
		//iterator on all vertices connected to a vertex
			return LinkTargetIterator(*this,out_edges_begin(vid));
		}

		out_neighbor_iterator out_neighbors_end (int vid) {
		//iterator on all vertices connected to a vertex
			return LinkTargetIterator(*this,out_edges_end(vid));
		}

		int nb_out_neighbors (int vid) const {
		//number of vertices connected to a vertex
			return nb_out_edges(vid);
		}

		neighbor_iterator neighbors_begin (int vid);
		//iterator on all vertices connected to a vertex

		neighbor_iterator neighbors_end (int vid);
		//iterator on all vertices connected to a vertex

		int nb_neighbors (int vid) const {
		//number of vertices connected to a vertex
			return nb_edges(vid);
		}


		//
		//
		// mutable concept
		//
		//
		int add_vertex () {
		//add a new vertex to the graph
		//return id used for this vertex
			return vertex_links.add(vpair(std::set<int>(),std::set<int>()));
		}

		int add_vertex (int vid) {
		//try to add a vertex to the graph using the provided id
		//if vid already used throw InvalidVertex
			try {
				return vertex_links.add(vpair(std::set<int>(),std::set<int>()),vid);
			}
			catch (IdMap<vpair>::InvalidId) {
				throw InvalidVertex(vid);
			}
		}

		void remove_vertex (int vid);
		//try to remove a vertex form the graph
		//if the vertex does not exist throw InvalidVertex

		int add_edge (int sid, int tid);
		//add a link between vertex sid and tid
		//throw InvalidVertex if sid or tid not in graph
		//return id used for this edge

		int add_link (int sid, int tid) {
		//same as add_edge
			return add_edge(sid,tid);
		}

		int add_edge (int sid, int tid, int eid);
		//try to add a link between vertex sid and tid
		//using eid as edge id
		//throw InvalidVertex if sid or tid not in graph
		//if eid already used, throw InvalidEdge
		
		int add_link (int sid, int tid, int eid) {
		//same as add_edge
			return add_edge(sid,tid,eid);
		}

		void remove_edge (int eid);
		//try to remove a link from the graph
		//throw InvalidEdge if the edge is not in the graph

		void remove_link (int lid) {
		//try to remove a link from the graph
		//throw InvalidEdge if the edge is not in the graph
			remove_edge(lid);
		}

		void clear_edges ();
		//remove all links from this graphs

		void clear ();
		//remove all vertices and links from the graph

		//debug
		void state () {
			Relation::state();
		}
	};
};
#endif //__CONTAINER_GRAPH_H__
