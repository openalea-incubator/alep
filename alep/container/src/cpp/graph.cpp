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

#include <iostream>
#include <vector>
#include "container/graph.h"
#include "container/graph_iterators.h"

namespace container {
	void Graph::test_existing_vertex (int vid) const {
		if(!has_vertex(vid)) {
			throw Graph::InvalidVertex(vid);
		}
	}

	void Graph::add_edge_to_vertex (int sid, int tid, int eid) {
		vertex_links.find(sid)->second.second.insert(eid);
		vertex_links.find(tid)->second.first.insert(eid);
	}

	Graph::vertex_edge_iterator Graph::edges_begin (int vid) {
		return EdgeIterator(in_edges_begin(vid),in_edges_end(vid),out_edges_begin(vid));
	}

	Graph::vertex_edge_iterator Graph::edges_end (int vid) {
		return EdgeIterator(out_edges_end(vid));
	}
	
	Graph::neighbor_iterator Graph::neighbors_begin (int vid) {
		return VertexIterator(*this,in_edges_begin(vid),in_edges_end(vid),out_edges_begin(vid));
	}

	Graph::neighbor_iterator Graph::neighbors_end (int vid) {
		return VertexIterator(*this,out_edges_end(vid));
	}
	
	void Graph::remove_vertex (int vid) {
		std::vector<int> tmp;
		for(Graph::vertex_edge_iterator it=edges_begin(vid);it!=edges_end(vid);++it) {
			tmp.push_back(*it);
		}
		for(std::vector<int>::iterator it=tmp.begin();it!=tmp.end();++it) {
			remove_edge(*it);
		}
		vertex_links.erase(vid);
	}

	int Graph::add_edge (int sid, int tid) {
		test_existing_vertex(sid);
		test_existing_vertex(tid);
		int eid=Relation::add_link(sid,tid);
		add_edge_to_vertex(sid,tid,eid);
		return eid;
	}

	int Graph::add_edge (int sid, int tid, int eid) {
		test_existing_vertex(sid);
		test_existing_vertex(tid);
		add_edge_to_vertex(sid,tid,Relation::add_link(sid,tid,eid));
		return eid;
	}

	void Graph::remove_edge (int eid) {
		vertex_links.find(source(eid))->second.second.erase(eid);
		vertex_links.find(target(eid))->second.first.erase(eid);
		Relation::remove_link(eid);
	}

	void Graph::clear_edges () {
		Relation::clear();
		for(Graph::v_iterator it=vertex_links.begin();it!=vertex_links.end();++it) {
			it->second.first.clear();
			it->second.second.clear();
		}
	}

	void Graph::clear () {
		Relation::clear();
		vertex_links.clear();
	}

};
