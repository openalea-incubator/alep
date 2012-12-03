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

#include <vector>
#include <map>
#include <set>

#include <boost/python.hpp>
#include "export_iterator.h"
using namespace boost::python;

void export_iterators () {
	export_custom_range<std::map<int,float>::iterator>("_PyMapIntFloatRange");
	export_custom_range<std::map<int,float>::const_iterator>("_PyConstMapIntFloatRange");
	export_custom_range<std::set<int>::iterator>("_PySetIntRange");
	//export_custom_range<std::set<int>::const_iterator>("_PyConstSetIntRange");
	export_custom_range<std::map<int,std::set<int> >::iterator>("_PyMapIntSetIntRange");
	//export_custom_range<std::map<int,std::set<int> >::const_iterator>("_PyConstMapIntSetIntRange");
}
