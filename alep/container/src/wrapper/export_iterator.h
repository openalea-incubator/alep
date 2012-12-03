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

#ifndef __CONTAINER_EXPORT_ITERATOR_H__
#define __CONTAINER_EXPORT_ITERATOR_H__

#include <boost/python/errors.hpp>
using namespace boost::python;

template <typename T_iterator>
class MapItemIterator {
private:
	typedef MapItemIterator<T_iterator> self;//for internal use
public:
	typedef tuple value_type;
public:
	T_iterator current_it;
public:
	MapItemIterator(T_iterator ref) : current_it(ref) {}
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

	bool operator== (const self& other) const {return current_it == other.current_it;}
	bool operator!= (const self& other) const {return current_it != other.current_it;}
};

template <typename T_iterator>
class PyCustomRange {
private:
	T_iterator it_begin;
	T_iterator it_end;
public:
	PyCustomRange(T_iterator ref_begin, T_iterator ref_end) : it_begin(ref_begin), it_end(ref_end) {
	}
	~PyCustomRange() {}
	typename T_iterator::value_type next ();
};

template <typename T_iterator>
typename T_iterator::value_type PyCustomRange<T_iterator>::next () {
	if(it_begin==it_end) {
		PyErr_SetString(PyExc_StopIteration,"done iterating");
		throw error_already_set();
	}
	else {
		typename T_iterator::value_type ret=*it_begin;
		++it_begin;
		return ret;
	}
}

template <typename T_iterator>
const PyCustomRange<T_iterator>& iter_func (const PyCustomRange<T_iterator>& it) {
	return it;
}

template <typename T_iterator>
void export_custom_range (const char* class_name) {
	class_<PyCustomRange<T_iterator> >(class_name,init<T_iterator,T_iterator>())
		.def("next",&PyCustomRange<T_iterator>::next)
		.def("__iter__",&iter_func<T_iterator>,return_internal_reference<1>());
}

#endif //__CONTAINER_EXPORT_ITERATOR_H__

