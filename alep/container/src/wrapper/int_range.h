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

#ifndef __CONTAINER_INT_RANGE_H__
#define __CONTAINER_INT_RANGE_H__
#include "container/int_iterator.h"
using namespace container;

#include <boost/python/errors.hpp>
#include <boost/python.hpp>
using namespace boost::python;

class PyIntRange {
private:
	IntIterator* it_begin;
	IntIterator* it_end;
public:
	PyIntRange(IntIterator* ref_begin, IntIterator* ref_end) : it_begin(ref_begin), it_end(ref_end) {
	}
	~PyIntRange() {}
	const PyIntRange& iter () const {
		return *this;
	}
	int next () {
		if(it_begin->operator==(*it_end)) {
			PyErr_SetString(PyExc_StopIteration,"done iterating");
			throw error_already_set();
		}
		else {
			int ret=**it_begin;
			it_begin->operator++();
			return ret;
		}
	}
};
#endif //__CONTAINER_INT_RANGE_H__

