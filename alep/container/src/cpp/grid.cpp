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
#include "container/grid.h"

namespace container {
	Grid::Grid (const Grid::coord_list& shape) {
		set_shape(shape);
	}

	void Grid::set_shape (const Grid::coord_list& shape) {
		max_coords.clear();
		offset_values.clear();
		for(Grid::coord_list::const_iterator it=shape.begin();it!=shape.end();++it) {
			max_coords.push_back(*it);
		}
		for(int i=0;i<max_coords.size();++i) {
			offset_values.push_back(1);
		}
		for(int i=0;i<(max_coords.size()-1);++i) {
			offset_values[i+1]=offset_values[i]*max_coords[i];
		}
	}

	int Grid::size () const {
		int s=1;
		Grid::coord_list::const_iterator it;
		for(it=max_coords.begin();it!=max_coords.end();++it) {
			s *= *it;
		}
		return s;
	}

	int Grid::index (const Grid::coord_list& coord) const {
		int ind=0;
		for(int i=0;i<offset_values.size();++i) {
			if(coord[i]<0) {
				throw Grid::OutOfBoundError(coord[i],0);
			}
			else {
				if(coord[i]>max_coords[i]) {
					throw Grid::OutOfBoundError(coord[i],max_coords[i]);
				}
				else {
					ind += coord[i]*offset_values[i];
				}
			}
		}
		return ind;
	}

	Grid::coord_list Grid::coordinates (int ind) const {
		if(ind<0) {
			throw Grid::OutOfBoundError(ind,0);
		}
		else {
			if(ind>size()) {
				throw Grid::OutOfBoundError(ind,size());
			}
			else {
				Grid::coord_list coord(dim());
				for(int i=(dim()-1);i>-1;--i) {
					coord[i]=ind/offset_values[i];
					ind%=offset_values[i];
				}
				return coord;
			}
		}
	}

	void Grid::state () const {
		std::cout << "shape :";
		for(int i=0;i<max_coords.size();++i) {
			std::cout << max_coords[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "size :" << size() << std::endl;
	}
};
