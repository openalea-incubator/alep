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

#include "export_custom_map.h"

void export_iterators();
void export_id_generator ();
void export_id_map ();
void export_grid ();
void export_relation ();
void export_graph ();
void export_topomesh ();
void export_int_range ();
void export_property_map ();

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(_container) {
	export_iterators();
	export_id_generator();
	export_id_map();
	export_grid();
	export_relation();
	export_graph();
	export_topomesh();
	export_int_range();
	export_custom_map<int,int>("Int");
	export_custom_map<int,float>("IntFloat");
	export_property_map();
}
