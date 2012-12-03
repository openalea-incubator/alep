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

#include <boost/python.hpp>
#include "int_range.h"
using namespace boost::python;

void export_int_range () {
	class_<PyIntRange>("_PyIntRange",init<IntIterator*,IntIterator*>())
		.def("next",&PyIntRange::next)
		.def("__iter__",&PyIntRange::iter,return_internal_reference<1>());
}

