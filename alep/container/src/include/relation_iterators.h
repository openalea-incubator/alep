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

#ifndef __CONTAINER_RELATION_ITERATORS_H__
#define __CONTAINER_RELATION_ITERATORS_H__

#include <set>
#include "container/relation.h"

namespace container {

	class LinkSourceIterator {
	private:
		typedef LinkSourceIterator self;//for internal use
		typedef std::set<int>::iterator local_it;
	public:
		typedef int value_type;
	public:
		local_it current_edge;
		Relation relation;
	public:
		LinkSourceIterator(const Relation& rel, local_it ref) {
			current_edge=ref;
			relation=rel;
		}
		~LinkSourceIterator() {}
		int operator* () const {
			return relation.source(*current_edge);
		}

		self& operator++ () {
			current_edge++;
			return *this;
		}

		self operator++ (int) {
			self ret=*this;
			++current_edge;
			return ret;
		}

		bool operator== (const LinkSourceIterator& other) const {return current_edge == other.current_edge;}
		bool operator!= (const LinkSourceIterator& other) const {return current_edge != other.current_edge;}
		//debug
		void state () const;
	};

	class LinkTargetIterator {
	private:
		typedef LinkTargetIterator self;//for internal use
		typedef std::set<int>::iterator local_it;
	public:
		typedef int value_type;
	public:
		local_it current_edge;
		Relation relation;
	public:
		LinkTargetIterator(const Relation& rel, local_it ref) {
			current_edge=ref;
			relation=rel;
		}
		~LinkTargetIterator() {}
		int operator* () const {
			return relation.target(*current_edge);
		}

		self& operator++ () {
			current_edge++;
			return *this;
		}

		self operator++ (int) {
			self ret=*this;
			++current_edge;
			return ret;
		}

		bool operator== (const LinkTargetIterator& other) const {return current_edge == other.current_edge;}
		bool operator!= (const LinkTargetIterator& other) const {return current_edge != other.current_edge;}
		//debug
		void state () const;
	};

};
#endif //__CONTAINER_RELATION_ITERATORS_H__
