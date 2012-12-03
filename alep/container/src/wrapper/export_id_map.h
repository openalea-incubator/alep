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

#include "container/id_map.h"
using namespace container;

#include <boost/python.hpp>
#include <boost/iterator.hpp>
#include "export_iterator.h"
using namespace boost::python;

/*template <typename D>
class MapItemIterator {
private:
	typedef MapItemIterator<D> self;//for internal use
	typedef typename IdMap<D>::const_iterator local_it;
public:
	typedef tuple value_type;
public:
	local_it current_it;
public:
	MapItemIterator(local_it ref) {
		current_it=ref;
	}
	~MapItemIterator() {}
	tuple operator* () const {
		return make_tuple(current_it->first,current_it->second);
	}
	self& operator++ () {
		current_it++;
		return *this;
	}
	self operator++ (int) {
		self ret=*this;
		++current_it;
		return ret;
	}

	bool operator== (const MapItemIterator<D>& other) const {return current_it == other.current_it;}
	bool operator!= (const MapItemIterator<D>& other) const {return current_it != other.current_it;}
	//debug
	void state () const;
};*/

#define ItemIt MapItemIterator<typename IdMap<D>::const_iterator>

template <typename D>
PyCustomRange< ItemIt > export_iteritems (const IdMap<D>& m) {
	return PyCustomRange< ItemIt >(ItemIt(m.begin()),ItemIt(m.end()));
}

template <typename D>
PyCustomRange<typename IdMap<D>::const_key_iterator> export_key_iter (const IdMap<D>& m) {
	return PyCustomRange<typename IdMap<D>::const_key_iterator>(m.key_begin(),m.key_end());
}

template <typename D>
D export_getitem (IdMap<D>& m, const int& key) {
	return m.getitem(key);
}

template <typename D>
D export_pop (IdMap<D>& m, const int& key) {
	D ret=m.getitem(key);
	m.erase(key);
	return ret;
}

template <typename D>
tuple export_popitem (IdMap<D>& m) {
	typename IdMap<D>::iterator it=m.begin();
	if(it==m.end()) {
		throw 0;
	}
	else {
		m.erase(it);
		return make_tuple(it->first,it->second);
	}
}

template <typename D>
IdMap<D> export_copy (const IdMap<D>& m) {
	return IdMap<D>(m);
}

template <typename D>
D export_setdefault (IdMap<D>& m, const int& key, D default_value) {
	if(m.has_key(key)) {
		return m.getitem(key);
	}
	else {
		m.add(default_value,key);
		return default_value;
	}
}

template <typename D>
int export_add (IdMap<D>& m, const D& val, PyObject* arg) {
	if(arg==Py_None) {
		return m.add(val);
	}
	else {
		return m.add(val,extract<int>(arg));
	}
}

template <typename D>
void export_id_map (char* name) {
	std::string key_it_name("_PyConstIdMap");
	key_it_name+=name;
	key_it_name+="Range";
	export_custom_range<typename IdMap<D>::const_key_iterator>(key_it_name.data());
	std::string item_it_name("_PyMapItem");
	item_it_name+=name;
	item_it_name+="Range";
	export_custom_range< ItemIt >(item_it_name.data());

	std::string class_name("IdMap");
	class_name+=name;
	class_< IdMap<D> >(class_name.data(), "dict with id generator")
		.def("__len__",&IdMap<D>::size,"number of elements")
		.def("__contains__",&IdMap<D>::has_key,"test wether the key is inside the dict")
		.def("iteritems",&export_iteritems<D>,"iterate on all pairs of key,val")
		.def("iterkeys",&export_key_iter<D>,"iterator on keys")
		.def("__iter__",&export_key_iter<D>,"iterator on keys")
		.def("__getitem__",&export_getitem<D>,"get value")
		.def("__setitem__",&IdMap<D>::setitem,"set value")
		.def("__delitem__",(void (IdMap<D>::*) (const int &))& IdMap<D>::erase,"remove a key from the dict")
		.def("pop",&export_pop<D>,"remove the value and return it")
		.def("popitem",&export_popitem<D>,"remove and return a random items from the dict")
		.def("setdefault",&export_setdefault<D>,"add an element inside the dictionary if their is nothing at this position")
		.def("copy",&export_copy<D>,"return a copy of the dict")
		.def("add",(int (IdMap<D>::*) (const D&))& IdMap<D>::add,"add a new element")
		.def("add",&export_add<D>,"add a new element")
		//.def("add",(int (IdMap<D>::*) (const D&, const int&))& IdMap<D>::add,"add a new element")
		.def("clear",&IdMap<D>::clear,"clear all elements in the map");
		//.def("state",&IdMap<D>::state,"debug function");
}

