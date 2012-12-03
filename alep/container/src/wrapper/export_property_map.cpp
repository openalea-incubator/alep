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

#include "container/property_map.h"
using namespace container;

#include <boost/python.hpp>
#include "export_custom_map.h"
#include "export_iterator.h"
using namespace boost::python;

PropertyMap& export_property_map_iadd (PropertyMap& m, float val) {
	m.iadd(val);
	return m;
}

PropertyMap& export_property_map_isub (PropertyMap& m, float val) {
	m.isub(val);
	return m;
}

PropertyMap& export_property_map_imul (PropertyMap& m, float val) {
	m.imul(val);
	return m;
}

PropertyMap& export_property_map_idiv (PropertyMap& m, float val) {
	m.idiv(val);
	return m;
}

PyCustomRange<FactorMap::key_iterator> export_factor_map_iterkeys (FactorMap& m) {
	return PyCustomRange<FactorMap::key_iterator>(m.key_begin(),m.key_end());
}

#define FactorItem MapItemIterator<FactorMap::iterator>

PyCustomRange< FactorItem > export_factor_map_iteritems (FactorMap& m) {
	return PyCustomRange< FactorItem >(FactorItem(m.begin()),FactorItem(m.end()));
}

void export_property_map () {
	class_<PropertyMap, bases< CustomMap<int,float> > >("PropertyMap", "special dict store a float value with an id")
		.def("__iadd__",&export_property_map_iadd,return_internal_reference<1>(),"add a value to all elements")
		.def("__isub__",&export_property_map_isub,return_internal_reference<1>(),"substract a value to all elements")
		.def("__imul__",&export_property_map_imul,return_internal_reference<1>(),"multiply all elements by a value")
		.def("__idiv__",&export_property_map_idiv,return_internal_reference<1>(),"divide all elements by a value");

	export_custom_range< FactorItem >("_PyFactorMapItemRange");
	export_custom_range< FactorMap::key_iterator >("_PyFactorMapKeyRange");


	class_<FactorMap>("FactorMap",init<FactorMap::associated_map&,float>("property map with units"))
		.def("__getitem__",&FactorMap::getitem,"retrieve a value")
		.def("__setitem__",&FactorMap::setitem,"set a value")
		.def("__iter__",&export_factor_map_iterkeys,"iterator on keys")
		.def("iterkeys",&export_factor_map_iterkeys,"iterator on keys")
		.def("iteritems",&export_factor_map_iteritems,"iterator on tuple (key,val)");
}

#undef FactorItem
