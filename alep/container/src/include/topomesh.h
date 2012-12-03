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

#ifndef __CONTAINER_TOPOMESH_H__
#define __CONTAINER_TOPOMESH_H__

#include <utility>
#include <map>
#include "container/id_map.h"
#include "container/relation.h"
#include "container/relation_iterators.h"

namespace container {
	class Topomesh : public Relation {
	private:
		typedef IdMap< std::set<int> >::iterator cp_iterator;
	public:
		typedef Relation::iterator link_iterator;//iter of int
		typedef std::set<int>::iterator cell_link_iterator;//iter of int
		typedef std::set<int>::iterator point_link_iterator;//iter of int
		typedef IdMap< std::set<int> >::key_iterator cell_iterator;//iter of int
		typedef LinkSourceIterator point_cell_iterator;//iter of int
		typedef int neighbor_cell_iterator;//iter of int
		typedef IdMap< std::set<int> >::key_iterator point_iterator;//iter of int
		typedef LinkTargetIterator cell_point_iterator;//iter of int
		typedef int neighbor_point_iterator;//iter of int

	public:
		typedef Relation::InvalidLink InvalidLink;
		struct InvalidCell : public IdMap< std::set<int> >::InvalidId {
		//small exception raised when attempting to access an invalid cell
			InvalidCell (int l) : IdMap< std::set<int> >::InvalidId(l) {};
			virtual ~InvalidCell() throw () {};
			InvalidCell* copy () { return new InvalidCell(id); }
			virtual void rethrow () { throw InvalidCell(id); }
		};
		struct InvalidPoint : public IdMap< std::set<int> >::InvalidId {
		//small exception raised when attempting to access an invalid point
			InvalidPoint (int l) : IdMap< std::set<int> >::InvalidId(l) {};
			virtual ~InvalidPoint() throw () {};
			InvalidPoint* copy () { return new InvalidPoint(id); }
			virtual void rethrow () { throw InvalidPoint(id); }
		};
	private:
		IdMap< std::set<int> > cell_links;//map of cell_id:set(out_links)
		IdMap< std::set<int> > point_links;//map of point_id:set(in_links)

	private:
		void test_existing_cell (int cid) const;
		//throw InvalidCell if cid not in the mesh

		void test_existing_point (int pid) const;
		//throw InvalidPoint if pid not in the mesh

		void add_edge_to_elements (int cid, int pid, int lid);
		//internal function to add lid to cell_links and point_links

	public:
		Topomesh() : Relation() {}
		~Topomesh() {}

		bool has_cell (int cid) const {
		//return true if the cell is inside the mesh
			return cell_links.find(cid)!=cell_links.end();
		}

		bool has_point (int pid) const {
		//return true if the point is inside the mesh
			return point_links.find(pid)!=point_links.end();
		}

		//
		//
		// link concept
		//
		//
		cell_link_iterator cell_links_begin (int cid) {
		//iterator on all edges connected to a cell
			test_existing_cell(cid);
			return cell_links.find(cid)->second.begin();
		}

		cell_link_iterator cell_links_end (int cid) {
		//iterator on all edges connected to a cell
			test_existing_cell(cid);
			return cell_links.find(cid)->second.end();
		}

		int nb_cell_links (int cid) {
		//number of edges connected to a cell
			test_existing_cell(cid);
			return cell_links.find(cid)->second.size();
		}

		point_link_iterator point_links_begin (int pid) {
		//iterator on all edges connected to a point
			test_existing_point(pid);
			return point_links.find(pid)->second.begin();
		}

		point_link_iterator point_links_end (int pid) {
		//iterator on all edges connected to a point
			test_existing_point(pid);
			return point_links.find(pid)->second.end();
		}

		int nb_point_links (int pid) {
		//number of edges connected to a point
			test_existing_point(pid);
			return point_links.find(pid)->second.size();
		}

		int cell (int lid) {
		//cell extremity of a  link
			return source(lid);
		}

		int point (int lid) {
		//point extremity of a  link
			return target(lid);
		}

		//
		//
		// cell list concept
		//
		//
		cell_iterator cells_begin () {
		//iterator on all cells
			return cell_links.key_begin();
		}

		cell_iterator cells_end () {
		//iterator on all cells
			return cell_links.key_end();
		}

		int nb_cells () {
		//number of cells in the mesh
			return cell_links.size();
		}

		point_cell_iterator cells_begin (int pid) {
		//iterator on all cells connected to a point
			return LinkSourceIterator(*this,point_links_begin(pid));
		}

		point_cell_iterator cells_end (int pid) {
		//iterator on all cells connected to a point
			return LinkSourceIterator(*this,point_links_end(pid));
		}

		int nb_cells (int pid) {
		//number of cells connected to a point
			return nb_point_links(pid);
		}

		neighbor_cell_iterator cell_neighbors_begin (int cid);
		//iterator on cells around a given cell

		neighbor_cell_iterator cell_neighbors_end (int cid);
		//iterator on cells around a given cell

		int nb_cell_neighbors (int cid);
		//number of cells around a given cell

		//
		//
		// point list concept
		//
		//
		point_iterator points_begin () {
		//iterator on all points
			return point_links.key_begin();
		}

		point_iterator points_end () {
		//iterator on all pointls
			return point_links.key_end();
		}

		int nb_points () {
		//number of points in the mesh
			return point_links.size();
		}

		cell_point_iterator points_begin (int cid) {
		//iterator on all points connected to a cell
			return LinkTargetIterator(*this,cell_links_begin(cid));
		}

		cell_point_iterator points_end (int cid) {
		//iterator on all points connected to a cell
			return LinkTargetIterator(*this,cell_links_end(cid));
		}

		int nb_points (int cid) {
		//number of points connected to a cell
			return nb_cell_links(cid);
		}

		neighbor_point_iterator point_neighbors_begin (int pid);
		//iterator on points around a given point

		neighbor_point_iterator point_neighbors_end (int pid);
		//iterator on points around a given point

		int nb_point_neighbors (int pid);
		//number of points around a given point


		//
		//
		// mutable concept
		//
		//
		int add_cell () {
		//add a new cell to the mesh
		//return id used for this cell
			return cell_links.add(std::set<int>());
		}

		int add_cell (int cid) {
		//try to add a cell to the mesh using the provided id
		//if cid already used throw InvalidCell
			try {
				return cell_links.add(std::set<int>(),cid);
			}
			catch (IdMap< std::set<int> >::InvalidId) {
				throw InvalidCell(cid);
			}
		}

		void remove_cell (int cid);
		//try to remove a cell form the mesh
		//if the cell does not exist throw InvalidCell

		int add_point () {
		//add a new point to the mesh
		//return id used for this point
			return point_links.add(std::set<int>());
		}

		int add_point (int pid) {
		//try to add a point to the mesh using the provided id
		//if pid already used throw InvalidPoint
			try {
				return point_links.add(std::set<int>(),pid);
			}
			catch (IdMap< std::set<int> >::InvalidId) {
				throw InvalidPoint(pid);
			}
		}

		void remove_point (int pid);
		//try to remove a point form the mesh
		//if the point does not exist throw InvalidPoint

		int add_link (int cid, int pid) {
		//add a link between a cell and a point
			test_existing_cell(cid);
			test_existing_point(pid);
			int lid = Relation::add_link(cid,pid);
			add_edge_to_elements(cid,pid,lid);
			return lid;
		}

		int add_link (int cid, int pid, int lid) {
		//try to add a link between a cell and a point
		//using the provided lid
		//throw InvalidLink if the id is already used
			test_existing_cell(cid);
			test_existing_point(pid);
			Relation::add_link(cid,pid,lid);
			add_edge_to_elements(cid,pid,lid);
			return lid;
		}

		void remove_link (int lid);
		//remove a link but not the cell and point

		void clear_links ();
		//remove all links from this mesh

		void clear ();
		//remove all cells, points and links from the mesh

		//debug
		void state () {
			Relation::state();
		}
	};
};
#endif //__CONTAINER_TOPOMESH_H__
