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

#include "container/id_generator.h"
using namespace container;

#include <boost/python.hpp>
using namespace boost::python;

int export_get_id (IdGenerator& gen, PyObject* arg) {
	if(arg==Py_None) {
		return gen.get_id();
	}
	else {
		return gen.get_id(extract<int>(arg));
	}
}

void export_id_generator () {
	class_<IdGenerator>("IdGenerator", "Generate integer ids")
		.def("get_id",(int (IdGenerator::*) ())& IdGenerator::get_id,"generate a free id")
		//.def("get_id",(int (IdGenerator::*) (int))& IdGenerator::get_id,"generate a free id")
		.def("get_id",&export_get_id,"generate a free id")
		.def("release_id",&IdGenerator::release_id,"free a no longer used id")
		.def("clear",&IdGenerator::clear,"reset id generator")
		.def("state",&IdGenerator::state,"debug function");
}

