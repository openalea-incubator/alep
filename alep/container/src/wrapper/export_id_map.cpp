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

#include <set>
#include "export_id_map.h"

typedef std::set<int> set;
#include <boost/python.hpp>
#include "export_iterator.h"
using namespace boost::python;

PyCustomRange<set::iterator> export_set_iter (set& s) {
	return PyCustomRange<set::iterator>(s.begin(),s.end());
}

void export_set_add (set& s, int val) {
	s.insert(val);
}

void export_set_remove (set& s, int val) {
	s.erase(val);
}

bool export_set_contains (const set& s, int val) {
	return s.find(val)!=s.end();
}

void export_id_map () {
	class_<set>("Set","simple export of std::set in python")
		.def("__iter__",&export_set_iter,"iterator on all elements of the set")
		.def("add",&export_set_add,"add a new element in the set")
		.def("remove",&export_set_remove,"remove an element from the set")
		.def("__contains__",&export_set_contains,"test wether an element is in the set")
		.def("clear",&set::clear,"empty the set");

	export_id_map<int>("Int");
	export_id_map<float>("Float");
	export_id_map<set>("SetInt");
}

