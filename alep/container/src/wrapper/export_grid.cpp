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

#include "container/grid.h"
using namespace container;

#include <boost/smart_ptr.hpp>
#include <boost/python.hpp>
#include "export_iterator.h"
using namespace boost::python;

typedef PyCustomRange<Grid::iterator> grid_index_range;

static boost::shared_ptr<Grid> export_grid_constructor (const object& shape) {
	tuple tup_shape=tuple(shape);
	Grid::coord_list vec_shape(len(tup_shape));
	for(int i=0;i<len(tup_shape);++i) {
		vec_shape[i]=extract<int>(tup_shape[i]);
	}
	return boost::shared_ptr<Grid>(new Grid(vec_shape));
}

tuple grid_shape (const Grid& g) {
	Grid::coord_list shp=g.shape();
	list l;
	Grid::coord_list::iterator it;
	for(it=shp.begin();it!=shp.end();++it) {
		l.append(*it);
	}
	return tuple(l);
}

void export_grid_set_shape (Grid& g, const object& shape) {
	tuple tup_shape=tuple(shape);
	Grid::coord_list vec_shape(len(tup_shape));
	for(int i=0;i<len(tup_shape);++i) {
		vec_shape[i]=extract<int>(tup_shape[i]);
	}
	g.set_shape(vec_shape);
}

grid_index_range export_grid_iter (const Grid& grid) {
	return grid_index_range(grid.begin(),grid.end());
}

int grid_index (const Grid& g, const object& coords) {
	tuple tup_coords=tuple(coords);
	Grid::coord_list coords_vec(len(tup_coords));
	for(int i=0;i<len(tup_coords);++i) {
		coords_vec[i]=extract<int>(tup_coords[i]);
	}
	return g.index(coords_vec);
}

tuple grid_coordinates (const Grid& g, int index) {
	Grid::coord_list coords=g.coordinates(index);
	list ret;
	for(Grid::coord_list::iterator it=coords.begin();it!=coords.end();++it) {
		ret.append(*it);
	}
	return tuple(ret);
}

void export_grid () {
	export_custom_range<Grid::iterator>("_PyGridIndexRange");

	class_<Grid>("Grid", "multi dimensional regular grid")
		//.def( init<const Grid::tuple&> ())
		.def("__init__",make_constructor(export_grid_constructor),"init")
		.def("dim",&Grid::dim,"dimension of space")
		.def("shape",&grid_shape,"number of cells in each dimension")
		.def("set_shape",&export_grid_set_shape,"set a new shape for this grid")
		.def("size",&Grid::size,"total number of cells in the grid")
		.def("__len__",&Grid::size,"total number of cells in the grid")
		.def("__iter__",&export_grid_iter,"iterate on all cells indexes")
		.def("index",&grid_index,"id of a specific cell in the grid")
		.def("coordinates",&grid_coordinates,"coordinates of a specific cell in the grid")
		.def("state",&Grid::state,"debug function");
}
