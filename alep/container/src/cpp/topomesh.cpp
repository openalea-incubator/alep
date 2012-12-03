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
#include "container/topomesh.h"

namespace container {
	void Topomesh::test_existing_cell (int cid) const {
		if(!has_cell(cid)) {
			throw Topomesh::InvalidCell(cid);
		}
	}

	void Topomesh::test_existing_point (int pid) const {
		if(!has_point(pid)) {
			throw Topomesh::InvalidPoint(pid);
		}
	}

	void Topomesh::add_edge_to_elements (int cid, int pid, int lid) {
		cell_links.find(cid)->second.insert(lid);
		point_links.find(pid)->second.insert(lid);
	}

	Topomesh::neighbor_cell_iterator Topomesh::cell_neighbors_begin (int cid) {
		test_existing_cell(cid);
		//TODO
	}

	Topomesh::neighbor_cell_iterator Topomesh::cell_neighbors_end (int cid) {
		test_existing_cell(cid);
		//TODO
	}

	int Topomesh::nb_cell_neighbors (int cid) {
		test_existing_cell(cid);
		//TODO
	}

	Topomesh::neighbor_point_iterator Topomesh::point_neighbors_begin (int pid) {
		test_existing_point(pid);
		//TODO
	}

	Topomesh::neighbor_point_iterator Topomesh::point_neighbors_end (int pid) {
		test_existing_point(pid);
		//TODO
	}

	int Topomesh::nb_point_neighbors (int pid) {
		test_existing_point(pid);
		//TODO
	}

	void Topomesh::remove_cell (int cid) {
		std::vector<int> tmp;
		for(Topomesh::cell_link_iterator it=cell_links_begin(cid);it!=cell_links_end(cid);++it) {
			tmp.push_back(*it);
		}
		for(std::vector<int>::iterator it=tmp.begin();it!=tmp.end();++it) {
			remove_link(*it);
		}
		cell_links.erase(cid);
	}

	void Topomesh::remove_point (int pid) {
		std::vector<int> tmp;
		for(Topomesh::point_link_iterator it=point_links_begin(pid);it!=point_links_end(pid);++it) {
			tmp.push_back(*it);
		}
		for(std::vector<int>::iterator it=tmp.begin();it!=tmp.end();++it) {
			remove_link(*it);
		}
		point_links.erase(pid);
	}

	void Topomesh::remove_link (int lid) {
		cell_links.find(cell(lid))->second.erase(lid);
		point_links.find(point(lid))->second.erase(lid);
		Relation::remove_link(lid);
	}

	void Topomesh::clear_links () {
		Relation::clear();
		for(Topomesh::cp_iterator it=cell_links.begin();it!=cell_links.end();++it) {
			it->second.clear();
		}
		for(Topomesh::cp_iterator it=point_links.begin();it!=point_links.end();++it) {
			it->second.clear();
		}
	}

	void Topomesh::clear () {
		Relation::clear();
		cell_links.clear();
		point_links.clear();
	}

};
