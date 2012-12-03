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

#include "container/relation.h"
#include "container/topomesh.h"
using namespace container;

#include <boost/python.hpp>
#include "export_iterator.h"
using namespace boost::python;

typedef PyCustomRange<Topomesh::cell_link_iterator> mesh_cell_link_range;
typedef PyCustomRange<Topomesh::point_link_iterator> mesh_point_link_range;
typedef PyCustomRange<Topomesh::cell_iterator> mesh_cell_range;
typedef PyCustomRange<Topomesh::point_cell_iterator> mesh_point_cell_range;
typedef PyCustomRange<Topomesh::point_iterator> mesh_point_range;
typedef PyCustomRange<Topomesh::cell_point_iterator> mesh_cell_point_range;

mesh_cell_link_range export_cell_links (Topomesh& m, int cid) {
	return mesh_cell_link_range(m.cell_links_begin(cid),m.cell_links_end(cid));
}
mesh_point_link_range export_point_links (Topomesh& m, int pid) {
	return mesh_point_link_range(m.point_links_begin(pid),m.point_links_end(pid));
}
mesh_cell_range export_cells (Topomesh& m) {
	return mesh_cell_range(m.cells_begin(),m.cells_end());
}
mesh_point_cell_range export_point_cells (Topomesh& m, int pid) {
	return mesh_point_cell_range(m.cells_begin(pid),m.cells_end(pid));
}
mesh_point_range export_points (Topomesh& m) {
	return mesh_point_range(m.points_begin(),m.points_end());
}
mesh_cell_point_range export_cell_points (Topomesh& m, int pid) {
	return mesh_cell_point_range(m.points_begin(pid),m.points_end(pid));
}

int export_mesh_add_cell (Topomesh& mesh, PyObject* arg) {
	if(arg==Py_None) {
		return mesh.add_cell();
	}
	else {
		return mesh.add_cell(extract<int>(arg));
	}
}

int export_mesh_add_point (Topomesh& mesh, PyObject* arg) {
	if(arg==Py_None) {
		return mesh.add_point();
	}
	else {
		return mesh.add_point(extract<int>(arg));
	}
}

int export_mesh_add_link (Topomesh& mesh, int elm1, int elm2, PyObject* arg) {
	if(arg==Py_None) {
		return mesh.add_link(elm1,elm2);
	}
	else {
		return mesh.add_link(elm1,elm2,extract<int>(arg));
	}
}

void export_topomesh () {
	export_custom_range<Topomesh::cell_iterator>("_PyTopomeshCellRange");

	class_<Topomesh, bases<Relation> >("Topomesh", "topomesh")
		//topomesh
		.def("has_cell",&Topomesh::has_cell,"test wether a cell is inside the mesh")
		.def("has_point",&Topomesh::has_point,"test wether a point is inside the mesh")
		//link
		.def("cell_links",&export_cell_links,"links connected to a cell")
		.def("nb_cell_links",&Topomesh::nb_cell_links,"number of links connected to a cell")
		.def("point_links",&export_point_links,"links connected to a point")
		.def("nb_point_links",&Topomesh::nb_point_links,"number of links connected to a point")
		.def("cell",&Topomesh::cell,"cell extremity of a link")
		.def("point",&Topomesh::point,"point extremity of a link")
		//cell list
		.def("cells",&export_cells,"iterator on cells")
		.def("cells",&export_point_cells,"iterator on cells")
		.def("nb_cells",(int (Topomesh::*) ())& Topomesh::nb_cells,"number of cells in the mesh")
		.def("nb_cells",(int (Topomesh::*) (int))& Topomesh::nb_cells,"number of cells in the mesh")
		//point list
		.def("points",&export_points,"iterator on points")
		.def("points",&export_cell_points,"iterator on points")
		.def("nb_points",(int (Topomesh::*) ())& Topomesh::nb_points,"number of points in the mesh")
		.def("nb_points",(int (Topomesh::*) (int))& Topomesh::nb_points,"number of points in the mesh")
		//mutable
		.def("add_cell",(int (Topomesh::*) ())& Topomesh::add_cell,"add a new cell in the mesh")
		//.def("add_cell",(int (Topomesh::*) (int))& Topomesh::add_cell,"add a new cell in the mesh")
		.def("add_cell",&export_mesh_add_cell,"add a new cell in the mesh")
		.def("remove_cell",&Topomesh::remove_cell,"remove a cell from the mesh")
		.def("add_point",(int (Topomesh::*) ())& Topomesh::add_point,"add a new point in the mesh")
		//.def("add_point",(int (Topomesh::*) (int))& Topomesh::add_point,"add a new point in the mesh")
		.def("add_point",&export_mesh_add_point,"add a new point in the mesh")
		.def("remove_point",&Topomesh::remove_point,"remove a point from the mesh")
		.def("add_link",(int (Topomesh::*) (int,int))& Topomesh::add_link,"add a new link between a cell and a point")
		//.def("add_link",(int (Topomesh::*) (int,int,int))& Topomesh::add_link,"add a new link between a cell and a point")
		.def("add_link",&export_mesh_add_link,"add a new link between a cell and a point")
		.def("remove_link",&Topomesh::remove_link,"remove a link, do not remove corresponding cell and point")
		.def("clear_links",&Topomesh::clear_links,"remove all links")
		.def("clear",&Topomesh::clear,"remove all cells ,points and links from the mesh")
		//debug
		.def("state",&Topomesh::state,"debug function");

}
