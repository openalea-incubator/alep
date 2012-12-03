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

#ifndef __CONTAINER_GRAPH_ITERATORS_H__
#define __CONTAINER_GRAPH_ITERATORS_H__

#include <set>
#include "container/graph.h"

namespace container {
	class EdgeIterator {
	private:
		typedef EdgeIterator self;//for internal use
	public:
		typedef int value_type;
		typedef std::set<int>::iterator local_it;
	public:
		local_it current_edge;
	private:
		local_it first_end;
		local_it second_begin;
	public:
		EdgeIterator () {}
		EdgeIterator(local_it ref, local_it ref_end=0, local_it ref_begin=0) {
			if(ref==ref_end) {
				current_edge=ref_begin;
				first_end=0;
				second_begin=0;
			}
			else {
				current_edge=ref;
				first_end=ref_end;
				second_begin=ref_begin;
			}
		}
		~EdgeIterator() {}
		int operator* () const {
			return *current_edge;
		}

		self& operator++ () {
			current_edge++;
			if(current_edge==first_end) {
				current_edge=second_begin;
			}
			return *this;
		}

		self operator++ (int) {
			self ret=*this;
			++current_edge;
			if(current_edge==first_end) {
				current_edge=second_begin;
			}
			return ret;
		}

		bool operator== (const EdgeIterator& other) const {return current_edge == other.current_edge;}
		bool operator!= (const EdgeIterator& other) const {return current_edge != other.current_edge;}
		//debug
		void state () const;
	};

	class VertexIterator {
	private:
		typedef VertexIterator self;//for internal use
	public:
		typedef int value_type;
		typedef std::set<int>::iterator local_it;
	public:
		local_it current_edge;
	private:
		local_it first_end;
		local_it second_begin;
		Graph graph;
		bool in_neighbors;
	public:
		VertexIterator(const Graph& g, local_it ref, local_it ref_end=0, local_it ref_begin=0) {
			if(ref==ref_end) {
				current_edge=ref_begin;
				first_end=0;
				second_begin=0;
				in_neighbors=false;
			}
			else {
				current_edge=ref;
				first_end=ref_end;
				second_begin=ref_begin;
				in_neighbors=true;
			}
			graph=g;
		}
		~VertexIterator() {}
		int operator* () const {
			if(in_neighbors) {
				return graph.source(*current_edge);
			}
			else {
				return graph.target(*current_edge);
			}
		}
		self& operator++ () {
			current_edge++;
			if(current_edge==first_end) {
				current_edge=second_begin;
				in_neighbors=false;
			}
			return *this;
		}

		self operator++ (int) {
			self ret=*this;
			++current_edge;
			if(current_edge==first_end) {
				current_edge=second_begin;
				in_neighbors=false;
			}
			return ret;
		}

		bool operator== (const VertexIterator& other) const {return current_edge == other.current_edge;}
		bool operator!= (const VertexIterator& other) const {return current_edge != other.current_edge;}
		//debug
		void state () const;
	};
};
#endif //__CONTAINER_GRAPH_ITERATORS_H__
