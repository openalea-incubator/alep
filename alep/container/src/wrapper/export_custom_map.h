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

#include "container/custom_map.h"
using namespace container;

#include <boost/python.hpp>
#include <boost/iterator.hpp>
#include "export_iterator.h"
using namespace boost::python;

#define CMap CustomMap<K,D>
#define CMapItem MapItemIterator<typename CustomMap<K,D>::const_iterator >

template <typename K, typename D>
PyCustomRange< CMapItem > export_iteritems (const CMap& m) {
	return PyCustomRange< CMapItem >(CMapItem(m.begin()),CMapItem(m.end()));
}

template <typename K, typename D>
PyCustomRange<typename CMap::const_key_iterator> export_key_iter (const CMap& m) {
	return PyCustomRange<typename CMap::const_key_iterator>(m.key_begin(),m.key_end());
}

template <typename K, typename D>
tuple export_popitem (CMap& m) {
	typename CMap::iterator it=m.popitem();
	return make_tuple(it->first,it->second);
}

template <typename K, typename D>
D export_getitem (CMap& m, const K& key) {
	return m.getitem(key);
}

template <typename K, typename D>
D export_setdefault (CMap& m, const K& key, const D& default_value) {
	return m.setdefault(key,default_value);
}

template <typename K, typename D>
str export_map_str (CMap& m) {
	str ret("{");
	for(typename CMap::iterator it=m.begin();it!=m.end();++it) {
		ret+="(";
		ret+=str(it->first);
		ret+=",";
		ret+=str(it->second);
		ret+=")";
		if(it!=m.end()) {
			ret==",";
		}
	}
	ret+="}";
	return ret;
}

template <typename K, typename D>
void export_custom_map (char* name) {
	std::string key_it_name("_PyConstCustomMap");
	key_it_name+=name;
	key_it_name+="Range";
	export_custom_range<typename CMap::const_key_iterator>(key_it_name.data());
	std::string item_it_name("_PyMapItem");
	item_it_name+=name;
	item_it_name+="Range";
	export_custom_range< CMapItem >(item_it_name.data());

	std::string class_name("Dict");
	class_name+=name;
	class_< CMap >(class_name.data(), "dict with id generator")
		.def("__len__",&CMap::size,"number of elements")
		.def("__contains__",&CMap::has_key,"test wether the key is inside the dict")
		.def("iteritems",&export_iteritems<K,D>,"iterate on all pairs of key,val")
		.def("iterkeys",&export_key_iter<K,D>,"iterator on keys")
		.def("__iter__",&export_key_iter<K,D>,"iterator on keys")
		.def("__getitem__",&export_getitem<K,D>,"get value")
		.def("__setitem__",&CMap::setitem,"set value")
		.def("__delitem__",&CMap::delitem,"remove a key from the dict")
		.def("pop",&CMap::pop,"remove the value and return it")
		.def("popitem",&export_popitem<K,D>,"remove and return a random items from the dict")
		.def("setdefault",&export_setdefault<K,D>,"add an element inside the dictionary if their is nothing at this position")
		.def("copy",&CMap::copy,"return a copy of the dict")
		.def("clear",&CMap::clear,"clear all elements in the map")
		.def("__str__",&export_map_str<K,D>,"string representation of the map");
}

#undef CMapItem
#undef CMap
